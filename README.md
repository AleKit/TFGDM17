# TFGDM17
Programas de Python utilizados en la realización del Trabajo de Fin de Grado de Física durante el curso 2016-2017.


**dm_spectra_f.py**
Este script contiene varias funciones:
- DMSpectrum 
- Sensitivities
- Funciones que dibujan la sensibilidad de distintos telescopios

**simuauto.py**
La función que contiene realiza simulaciones del modelo de branones para una masa y un tiempo de observación dados hasta que la significancia de detección alcanza un valor comprendido entre 24 y 26, e imprime la sección eficaz de aniquilación para la que ocurre esto. Si no se alcanzan estos valores de la significancia en un tiempo razonable, de modo que el factor por el que se multiplica o divide se vuelve demasiado pequeño, se toma el valor de la sección eficaz ...............................

**autoreal.py** 
Este script ejecuta todos los pasos de la simulación y el análisis repetidamente para distintos valores de la masa del branón. Las variables de entrada son el conjunto de masas, ALGO, el tiempo de observación y el número que queremos que aparezca en el nombre de la tabla que se obtiene.
