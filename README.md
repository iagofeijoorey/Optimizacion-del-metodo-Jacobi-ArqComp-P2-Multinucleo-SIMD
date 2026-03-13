# ArqComp_P2: Optimizacion del metodo Jacobi -> Multinúcleo & SIMD
Proyecto desarrollado para la asignatura de Arquitectura de Computadores del grado de Ingeniería Informática en la Universidad de Santiago de Compostela (USC). El equipo de desarrollo está compuesto por Iago Feijóo Rey y Nara Estévez.

## Descripción del Proyecto
Implementación y análisis de rendimiento del método iterativo de Jacobi en C, explorando cuatro niveles de optimización: versión secuencial base, optimización de caché (división de lazos, desenrollamiento, operaciones por bloques), paralelismo con OpenMP y vectorización SIMD con AVX2/FMA. Mediciones realizadas en FinisTerrae III (CESGA) para distintos tamaños de matriz.

## Contenido del Repositorio
* `src/`: Contiene el código fuente de las cuatro versiones del algoritmo.
    - `v1.c`: Versión secuencial base
    - `v2.c`: Versión optimizada con la mejor combinación de las siguientes técnicas de caché
    `v2_options/`: Contiene el código de las 4 variantes individuales de la versión 2 del código.
        - `v2_1.c`: Versión optimizada con reducción del número total de instrucciones a ejecutar.
        - `v2_2.c`: Versión optimizada con división lazos.
        - `v2_3.c`: Versión optimizada con desenrollamiento de lazos.
        - `v2_4.c`: Versión optimizada con realización de operaciones por bloques
    - `v3.c`: Versión paralelizada con OpenMP
    - `v4.c`: Versión vectorizada con AVX2/FMA
* `documentacion/`: Contiene la memoria de la práctica y otros documentos pertinentes.
* `makefile`: Fichero de compilación y ejecución.
* `counter.h`: Cabecera para la medición de ciclos de reloj.
* `README.md`: Este archivo.

## Requisitos del Sistema
* GCC con soporte para OpenMP y AVX2
* Sistema Linux (recomendado: FinisTerrae III, CESGA)

## Uso
```bash
# Compilar todas las versiones
make CFLAGS="-O2"

# Ejecutar todas las versiones con n=1250
make run ARGS="1250"

# Ejecutar con n=1250 y 4 hilos (OpenMP)
make run ARGS="1250 4"
```

## Contacto
Para cualquier duda o sugerencia, puedes ponerte en contacto con el equipo de desarrollo:

**Iago Feijóo Rey**: iagofeijoorey@gmail.com
* [LinkedIn](https://www.linkedin.com/in/iagofeijoorey/)

**Nara Estévez**: nara.estevez@rai.usc.gal

## Languages and Tools used
<p align="left">
<a href="https://www.gnu.org/software/gcc/" target="_blank" rel="noreferrer">
<img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/gcc/gcc-original.svg" alt="gcc" width="40" height="40"/></a>
<a href="https://www.openmp.org/" target="_blank" rel="noreferrer">
<img src="https://upload.wikimedia.org/wikipedia/commons/e/eb/OpenMP_logo.png" alt="openmp" width="40" height="40"/></a>
<a href="https://code.visualstudio.com/" target="_blank" rel="noreferrer">
<img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/vscode/vscode-original.svg" alt="vscode" width="40" height="40"/></a>
</p>
