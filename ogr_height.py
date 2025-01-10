from pyproj import CRS, Transformer
from osgeo import gdal
import pyproj.exceptions
import math
import multiprocessing
import alive_progress
import pathos.pools
import pandas
import shapely
import geopandas
import json

gdal.UseExceptions()

def get_raster_value(geo_x: float, geo_y: float, ds: gdal.Dataset, band_index: int = 1):
    """Return raster value that corresponds to given coordinates."""
    forward_transform = ds.GetGeoTransform()
    reverse_transform = gdal.InvGeoTransform(forward_transform)
    pixel_coord = gdal.ApplyGeoTransform(reverse_transform, geo_x, geo_y)
    pixel_x = math.floor(pixel_coord[0])
    pixel_y = math.floor(pixel_coord[1])
    band: gdal.Band = ds.GetRasterBand(band_index)
    val_arr = band.ReadAsArray(pixel_x, pixel_y, 1, 1) # Avoid reading the whole raster into memory - read 1x1 array
    return val_arr[0][0]

def getElevation(ds, from_point, to_point, dist, width, method, seg_id):
  #src = "DEM.tif"
  #from_point = [16.424878290636123, 39.38143932639117]
  #to_ponit = [16.525043215589466, 39.386377387105306]
  #in m
  #dist = 1
  #in m
  #width = 1
  ds = gdal.Open(ds)
  # Set up coordinate transform
  #print(from_point)
  #print(to_point)
  proj_str = "+proj=tpeqd +lon_1={} +lat_1={} +lon_2={} +lat_2={}".format(
    from_point[0], from_point[1],
    to_point[0], to_point[1])
  try:
    tpeqd = CRS.from_proj4(proj_str)
  #Seems to be the case for pretty small sections, that float64 explodes
  except pyproj.exceptions.CRSError:
    return {
      "data": [
        get_raster_value(from_point[0], from_point[1], ds),
        get_raster_value(to_point[0], to_point[1], ds)
      ]
    }
  try:
    transformer = Transformer.from_crs(CRS.from_proj4("+proj=latlon"), tpeqd)
  except pyproj.exceptions.ProjError:
    return {
      "data": [
        get_raster_value(from_point[0], from_point[1], ds),
        get_raster_value(to_point[0], to_point[1], ds)
      ]
    }
  # Transfor to tpeqd coordinates
  point_1 = transformer.transform(from_point[0], from_point[1])
  point_2 = transformer.transform(to_point[0], to_point[1])
  # Create an bounding box (minx, miny, maxx, maxy) in tpeqd coordinates
  bbox = (point_1[0], -(width*0.5), point_2[0], (width*0.5))
  # Calculate the number of samples in our profile.
  if ((point_2[0] - point_1[0]) == 0):
    return {
      "data": [
        get_raster_value(from_point[0], from_point[1], ds),
        get_raster_value(to_point[0], to_point[1], ds)
      ]
    }
  num_samples = max(2, int(math.ceil((point_2[0] - point_1[0]) / dist)))
  if ((point_2[0] - point_1[0])/num_samples) == 0:
    return {
      "data": [
        get_raster_value(from_point[0], from_point[1], ds),
        get_raster_value(to_point[0], to_point[1], ds)
      ]
    }
  dx = (point_2[0] - point_1[0])/num_samples
  # Warp it into dataset in tpeqd projection. If args.tif is empty GDAL will
  # interpret it as an in-memory dataset.
  tmp_ds = "/vsimem/test_{}.tif".format(seg_id)
  profile = gdal.Warp(tmp_ds, ds, dstSRS=proj_str, outputBounds=bbox, 
                      height=1, width=num_samples, resampleAlg=method)
  # Extract the pixel values and write to an output file
  data = profile.GetRasterBand(1).ReadAsArray()[0]
  gdal.Unlink(tmp_ds)
  del ds
  del profile
  return {"data": data.tolist(), "dx": dx}

import math
import numpy

def processLine(line, dem, dist, width, method, seg_id):
  ret = {
    "n_segments": len(line.coords) - 1,
    "elevation": []
  }

  for coord in range(len(line.coords) - 1):
    elevation = getElevation(dem, line.coords[coord], line.coords[coord + 1], dist, width, method, "{}_{}".format(seg_id, coord))
    ret["elevation"].append(elevation["data"])

  ret["elevation"] = json.dumps(ret["elevation"])

  return ret

def countSegments(gdf):
  r = 0
  for i in gdf["geometry"]:
    r += len(i.coords) - 1
  return r

def processVectors(vectors, dem, dist, width, method, cores = 1):
  if "id" not in vectors.columns: raise NameError("No Id column on vectors")
  if not vectors["id"].is_unique: raise NameError("The Id column is not unique")
  pool = pathos.pools.ProcessPool(nodes = cores)
  id = 0
  df = {}
  df["id"] = []
  df["elevation"] = []
  with alive_progress.alive_bar(countSegments(vectors)) as bar:
    for data in pool.imap(
      lambda x: processLine(x, dem, dist, width, method, str(id)),
      vectors["geometry"].tolist()
    ):
    #for data in map(lambda x: processLine(x, dem, dist, width, method, str(id)), vectors["geometry"].tolist()):
      df["id"].append(vectors["id"][id])
      df["elevation"].append(data["elevation"])
      bar(data["n_segments"])
      id += 1
  return pandas.DataFrame(df)



def tests():
  #ds = gdal.Open("Float64_Calabria_32633.tif")
  ds = "Float64_Calabria_32633.tif"
  from_point = [16.424878290636123, 39.38143932639117]
  to_point = [16.525043215589466, 39.386377387105306]
  dist = 0.1
  width = 0.1
  method = "bilinear"
  line = shapely.LineString([from_point, to_point])

  #test vectors
  network = geopandas.read_file("sample.gpkg").to_crs(4326)
  ret = processVectors(network, ds, dist, width, method, cores = 10)
  #n = network.merge(ret, on='id', how='left')

  #test elevation
  ret = getElevation(ds, from_point, to_point, dist, width, method = method)

  #test processLine
  ret = processLine(
    shapely.LineString([from_point, to_point]),
    dem,
    1,
    1,
    "bilinear"
  )

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description="Make profile with GDAL")
  parser.add_argument('DEM', metavar='DEM', help='DEM data source')
  parser.add_argument('Network', metavar='network', help='Network data source')
  parser.add_argument('output', metavar='output', help='SQLite db to save')
  parser.add_argument('--width', type=float, default=0.1, help='Profile width (m)')
  parser.add_argument('--dist', type=float, default=0.1, help='Profile sampling distance (m)')
  parser.add_argument('--resample', type=str, default='bilinear', help='Resampling method')
  parser.add_argument('--cores', type=int, default=1, help='Cores for parallel')
  args = parser.parse_args()
  ds = gdal.Open(args.DEM)
  if ds.GetRasterBand(1).DataType != gdal.GDT_Float64:
    raise NameError("The Raster must have Float64 as a type for interpolation")
  del ds
  ret = processVectors(
    geopandas.read_file(args.Network).to_crs(4326),
    args.DEM,
    args.dist,
    args.width,
    args.resample,
    cores = args.cores
  )
  import sqlite3
  conn = sqlite3.connect(args.output)
  ret.to_sql("height", conn, index = False)
  conn.close()