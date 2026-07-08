# Activity 3: Functional Differential Expression Analysis

Author: **Javier Luque Serrano**
Institution: **European University of Madrid**

## Description

This repository contains the code and results for Activity 3 of Module 4 (Research Methodology I), focused on functional differential expression analysis using the **TCGA-BRCA** breast cancer dataset (The Cancer Genome Atlas).

The goal is to identify and classify genes whose expression differs significantly between healthy and tumor tissue, in order to better understand the molecular mechanisms involved in breast cancer proliferation.

## Repository Contents

- **`analisis_expresion_diferencial.Rmd`**: R Markdown document with the full workflow (preprocessing, differential expression analysis, functional enrichment, and visualization).
- **`actividad.R`**: R script with the complete code to reproduce the analysis.
- **`analisis_expresion_diferencial_JavierLuque.pdf`**: Final PDF report generated from the Rmd.
- **`0JCP001104_UA5_AA1.pdf`**: Assignment guidelines / academic document.
- **`clinical_info_TCGA-BRCA.tsv`**: Clinical information of the TCGA-BRCA samples.
- **`nationwidechildrens.org_clinical_patient_brca.txt`**: Additional clinical data for BRCA patients.
- **`BRCA_exp_matrix.tsv`**: Raw expression matrix for the TCGA-BRCA dataset.
- **`data_norm.tsv`**: Normalized expression matrix (TMM normalization via `edgeR`).
- **`p_values.tsv`**: Table with significance values (p-values) from the differential expression analysis.
- **`last_image.RData`**: Saved R session image with analysis objects.

## Requirements

To run the analysis you need **R** (>= 4.0) and the following Bioconductor/CRAN libraries:

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

## How to Use

1. Clone the repository.
2. Open the project in RStudio or set the working directory to the repository folder.
3. Run the `actividad.R` script step by step, or knit the `analisis_expresion_diferencial.Rmd` file to generate the full report.

## License

This work is part of the Master's in Bioinformatics academic training and is shared for educational purposes.

---

# Actividad 3: Análisis de expresión diferencial funcional

Autor: **Javier Luque Serrano**
Institución: **Universidad Europea de Madrid**

## Descripción

Este repositorio contiene el código y los resultados de la Actividad 3 del Módulo 4 (Metodología de la Investigación I), centrada en el análisis de expresión diferencial funcional sobre el conjunto de datos de cáncer de mama **TCGA-BRCA** (The Cancer Genome Atlas).

El objetivo es identificar y clasificar los genes cuya expresión difiere significativamente entre tejido sano y tejido tumoral, para comprender mejor los mecanismos moleculares implicados en la proliferación del cáncer de mama.

## Contenido del repositorio

- **`analisis_expresion_diferencial.Rmd`**: Documento R Markdown con todo el flujo de trabajo (preprocesado, análisis de expresión diferencial, enriquecimiento funcional y visualización).
- **`actividad.R`**: Script de R con el código completo para reproducir el análisis.
- **`analisis_expresion_diferencial_JavierLuque.pdf`**: Informe final en PDF generado a partir del Rmd.
- **`0JCP001104_UA5_AA1.pdf`**: Documento de instrucciones de la actividad.
- **`clinical_info_TCGA-BRCA.tsv`**: Información clínica de las muestras TCGA-BRCA.
- **`nationwidechildrens.org_clinical_patient_brca.txt`**: Datos clínicos adicionales de pacientes BRCA.
- **`BRCA_exp_matrix.tsv`**: Matriz de expresión en bruto del conjunto de datos TCGA-BRCA.
- **`data_norm.tsv`**: Matriz de expresión normalizada (normalización TMM mediante `edgeR`).
- **`p_values.tsv`**: Tabla con los valores de significancia (p-values) del análisis de expresión diferencial.
- **`last_image.RData`**: Imagen guardada de la sesión de R con los objetos del análisis.

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
