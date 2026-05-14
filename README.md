# Actividad 3: Análisis de expresión diferencial funcional

Autor: **Javier Luque Serrano**

## Descripción

Este repositorio contiene el código y los resultados de la Actividad 3 del Módulo 4 (Metodología de la Investigación I), centrada en el análisis de expresión diferencial funcional sobre el conjunto de datos de cáncer de mama **TCGA-BRCA** (The Cancer Genome Atlas).

El objetivo es identificar y clasificar los genes cuya expresión difiere significativamente entre tejido sano y tejido tumoral, para comprender mejor los mecanismos moleculares implicados en la proliferación del cáncer de mama.

## Contenido del repositorio

- **`analisis_expresion_diferencial.Rmd`**: Documento R Markdown con todo el flujo de trabajo (preprocesado, análisis de expresión diferencial, enriquecimiento funcional y visualización).
- **`actividad.R`**: Script de R con el código completo para reproducir el análisis.
- **`analisis_expresion_diferencial_JavierLuque.pdf`**: Informe final en PDF generado a partir del Rmd.
- **`clinical_info_TCGA-BRCA.tsv`**: Información clínica de las muestras TCGA-BRCA.
- **`nationwidechildrens.org_clinical_patient_brca.txt`**: Datos clínicos adicionales de pacientes BRCA.
- **`p_values.tsv`**: Tabla con los valores de significancia (p-values) del análisis de expresión diferencial.

> **Nota**: Los archivos de datos de gran tamaño (`BRCA_exp_matrix.tsv` y `data_norm.tsv`) no están incluidos en el repositorio por limitaciones de tamaño de GitHub, pero se pueden regenerar ejecutando el script `actividad.R`.

## Requisitos

Para ejecutar el análisis se necesita **R** (>= 4.0) y las siguientes librerías de Bioconductor/CRAN:

```r
library(edgeR)
library(limma)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(GOSemSim)
library(org.Hs.eg.db)
library(DOSE)
library(clusterProfiler)
library(dplyr)
library(AnnotationDbi)
library(knitr)
library(tinytex)
```

## Cómo usar

1. Clonar el repositorio.
2. Abrir el proyecto en RStudio o establecer el directorio de trabajo en la carpeta del repositorio.
3. Ejecutar el script `actividad.R` paso a paso, o compilar el archivo `analisis_expresion_diferencial.Rmd` para generar el informe completo.

## Licencia

Este trabajo es parte de la formación académica del Máster en Bioinformática y se comparte con fines educativos.
