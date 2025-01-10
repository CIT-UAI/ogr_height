# ogr_height

Script que al ingresar un archivo vectorial y un DEM obtienes la altura de cada segmento de las calles.

El DEM debe estar en 4326.
El archivo vectorial debe tener la columna "id"

La salida sera un archivo SQLite, donde la tabla contendra dos columnas.

"id": Id llave del vector.
"elevation": Elevación del vector, es un JSON el cual contiene primero un vector, cada elemento es un segmento, y luego otra lista con los valores interpolados de la altura.

Para ver los parametros:
```
python3 ogr_height.py --help
```
