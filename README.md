# TFGDM17
Programas de Python utilizados en la realización del Trabajo de Fin de Grado de Física durante el curso 2016-2017.


**branchingratios.py**
Calcula las ratios de ramificación (contribuciones significativas a partir de 2 GeV) y la tensión de la brana para un rango de masas dado y las dibuja.

**dm_spectra_f.py**
Este script contiene varias funciones, entre las que se encuentran:
- **getDMspectrum**  Genera una tabla con el modelo de branones requerido. Uno de los parámetros de entrada es finalstate, que hace referencia al canal de desintegración. Si se quiere utilizar el modelo en el que se recogen todos los canales significativos hay que asignarle el valor 'new', de modo que empleará una modificación de **branchingratios** en el cálculo. Los demás parámetros de entrada son la masa del branón, el factor astrofísico, el \textit{boost} y las magnitudes que compondrán la tabla.
- Funciones que obtienen la sensibilidad de distintos telescopios para posteriormente pintarlas en un gráfico y comparar con la densidad espectral de un modelo.

**simuauto.py**
La función que contiene, **analisis**, realiza simulaciones del modelo de branones para una masa y un tiempo de observación dados hasta que la significancia de detección alcanza un valor comprendido entre 24 y 26, y guarda en un archivo de texto la sección eficaz de aniquilación para la que ocurre esto. Si no se alcanzan estos valores de la significancia en un tiempo razonable, de modo que el factor por el que se multiplica o divide se vuelve demasiado pequeño, se toma el valor de la sección eficaz que se haya empleado en la última simulación. Las variables de entrada son la sección eficaz inicial, el conjunto de masas, el tiempo de observación y el número que queremos que aparezca en el nombre de la tabla que se obtiene.

**autoreal.py** 
Este script ejecuta todos los pasos de la simulación y el análisis repetidamente para distintos valores de la masa del branón. Para ello utiliza los dos scripts anteriores. Las variables de entrada son las mismas que las de la función **analisis**.

**plotsigmav.py**
Tras proporcionarle las tablas necesarias, este script dibujará el gráfico en el que se mostrará la sección eficaz en función de la masa.
