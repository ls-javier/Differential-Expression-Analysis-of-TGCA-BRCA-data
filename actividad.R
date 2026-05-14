### Actividad Unidad 3, Modulo 4
## Analisis de expresion diferencial funcional 
# Nombre y apellidos: JAVIER LUQUE SERRANO

# Reproducible working directory: set this to the folder containing the script
# before running, or open the project so that relative file paths resolve correctly.
# Example: setwd("path/to/this/folder")
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
}

# NormalizaciC3n de la matriz de expresion (incluido en limma)
library(edgeR)
# Analisis de expresion diferencial
library(limma)
# Descargar los datos de TCGA
library(TCGAbiolinks)
library(SummarizedExperiment)
# Hacer grC!ficas
library(ggplot2)
# Complementario a ggplot2: graficas mas detalladas
library(ggrepel)
# Heatmap
library(pheatmap)
# Similitud semC!ntica
library(GOSemSim)
# Base de datos
library(org.Hs.eg.db)
# Libreria Disease Ontology
library(DOSE)
# Enriquecimiento funcional
library(clusterProfiler)
# Manejo de dataframes
library(dplyr)
# Vias metabolicas en KEGG
library(AnnotationDbi)
# Compilar RMarkDown como PDF
library(knitr)
library(tinytex)

# Ver la lista de proyectos que hay en GDC
gdcprojects <- getGDCprojects() # Lista de proyectos disponibles
getProjectSummary('TCGA-BRCA') # Resumen del proyecto TCGA-BRCA

# Construimos una query con la informaciC3n que queremos del proyecto
query_TCGA <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts" # Datos crudos de la expresiC3n de los genes en cada muestra
)

# Ver un resumen del objeto query_TCGA, el cual es una lista
resultados <- getResults(query_TCGA)

# Descargamos los datos
GDCdownload(query_TCGA)

# Preparar los datos en un objeto SummarizedExperiment
BRCA_data <- GDCprepare(query_TCGA, summarizedExperiment = T)

# Acceder a la matriz de expresiC3n (conteos)
counts_matrix <- assay(BRCA_data)

# Acceder a la informaciC3n de los genes como un dataframe
gene_info <- as.data.frame(rowData(BRCA_data))

# Acceder a la informaciC3n clC-nica y de las muestras como un dataframe
clinical_info <- as.data.frame(colData(BRCA_data))

### Usar datos de la actividad
library(readr)

clinical_info <- read_tsv("clinical_info_TCGA-BRCA.tsv", show_col_types = FALSE,
                          na=c("NA","","[Not Applicable]","[Not Available]","[Unknown]",
                               "[Not Evaluated]","Indeterminate","[Discerpancy]"))

brca_exp_mat <- read_tsv("BRCA_exp_matrix.tsv", show_col_types = FALSE,
                         na=c("NA","","[Not Applicable]","[Not Available]","[Unknown]",
                              "[Not Evaluated]","Indeterminate","[Discerpancy]"))

# Convertirlos en data frames, que son mas manejables que el formato tibble
brca_exp_mat <- as.data.frame(brca_exp_mat)
clinical_info <- as.data.frame(clinical_info)

## Exploracion del dataset
dim(brca_exp_mat)
dim(clinical_info)

# Comprobar primeras filas y columnas
head(brca_exp_mat[, 1:5]) # Primera columna es una muestra, hay que desplazar el dataset 
#una columna a la derecha
head(brca_exp_mat[, 1150:1155])

head(clinical_info[, 1:5]) # Primera columna es una muestra, hay que desplazar el dataset 
#una columna a la derecha

# Comprobar duplicados 
table(duplicated(colnames(brca_exp_mat)))
table(duplicated(brca_exp_mat$gene_ID))
# Comprobar valores NA en los nombres de las filas
sum(is.na(brca_exp_mat[, 1]))

# Cargamos librerias
library(tidyr) # para modificar la ultima columna
library(stringr) # para eliminar el espacio que hay alrededor del \t

# Guardar las columnas que tenemos
current_colnames <- colnames(brca_exp_mat)
# Definir la ultima columna para separar los valores
last_col_name <- current_colnames[length(current_colnames)]
# Limpiar el espacio en blanco alrededor del tabulador
brca_exp_mat[[last_col_name]] <- str_trim(brca_exp_mat[[last_col_name]])
brca_exp_mat[[last_col_name]] <- str_replace_all(brca_exp_mat[[last_col_name]], 
                                                 "\\s*\t\\s*", "\t")
# Separar la ultima columna en dos con nombres temporales
brca_exp_mat <- brca_exp_mat %>%
  separate_wider_delim(
    cols = all_of(last_col_name),
    delim = "\t",
    names = c("second_to_last", "very_last"),
    too_few = "align_start"
  )

# Desplazar las columnas hacia la derecha y nombrar a la columna 1 "gene_ID"
old_names <- current_colnames[1:(length(current_colnames) - 2)]
new_colnames <- c("gene_ID", old_names, 
                  current_colnames[length(current_colnames) - 1],
                  current_colnames[length(current_colnames)])

colnames(brca_exp_mat) <- new_colnames

brca_exp_mat <- as.data.frame(brca_exp_mat)

# Establecer la columna de genes como el nombre de las filas
rownames(brca_exp_mat) <- brca_exp_mat[, 1]
brca_exp_mat <- brca_exp_mat[, -1]

# Verificar
dim(brca_exp_mat)
head(colnames(brca_exp_mat))
tail(colnames(brca_exp_mat))

### Desplazar los nombres de las columnas hacia la derecha
# Coger las n-1 columnas
col_names_correct <- colnames(clinical_info)[1:(length(colnames(clinical_info))-1)]
# Crear nombres nuevos: un marcador para la primera columna y el resto desplazados
temporal_names <- c("delete", col_names_correct)
# Asignar los nuevos nombres a las columnas
colnames(clinical_info) <- temporal_names
# Eliminar la primera columna
clinical_info <- clinical_info[, -1]

# Verificar
head(clinical_info[, 1:5])
head(clinical_info[, 110: 113])

# Comprobar si hay genes duplicados en la informacion clinica
table(duplicated(clinical_info$sample))
# Eliminamos los duplicados con el operador inverso ! y el vector logico de duplicados
clinical_info_filtered <- clinical_info[!duplicated(clinical_info$sample),]

# Comprobar que todas las muestras de clinical_info estC)n en brca_exp_mat y viceversa
table(clinical_info_filtered$sample%in%colnames(brca_exp_mat))
table(colnames(brca_exp_mat)%in%clinical_info_filtered$sample)

# Ver las 17 muestras que sobran en clinical_info respecto a brca_exp_mat
clinical_patients <- clinical_info_filtered$sample
clinical_patients[!clinical_patients %in% colnames(brca_exp_mat)]

# Eliminar las muestras sobrantes
clinical_info_filtered <- clinical_info[clinical_info$sample %in% colnames(brca_exp_mat), ]

# Comprobar que tienen la misma dimension y estan en el mismo orden
all(colnames(brca_exp_mat)==clinical_info_filtered$sample)

# Comprobar primeras filas y columnas
head(brca_exp_mat[, 1:5])

# Comprobamos las dimensiones de las matrices y sus clases
dim(clinical_info_filtered) # Filas: muestras; columnas: variables clinicas
class(clinical_info_filtered)

dim(brca_exp_mat)
class(brca_exp_mat)

# Ver un resumen de ambos datasets
summary(brca_exp_mat[, 1:5])
summary(clinical_info_filtered[, 1:5])# Comprobar si hay genes duplicados en la informacion clinica
table(duplicated(clinical_info$sample))
# Eliminamos los duplicados con el operador inverso ! y el vector logico de duplicados
clinical_info_filtered <- clinical_info[!duplicated(clinical_info$sample),]

# Comprobar que todas las muestras de clinical_info estC)n en brca_exp_mat y viceversa
table(clinical_info_filtered$sample%in%colnames(brca_exp_mat))
table(colnames(brca_exp_mat)%in%clinical_info_filtered$sample)

# Ver las 17 muestras que sobran en clinical_info respecto a brca_exp_mat
clinical_patients <- clinical_info_filtered$sample
clinical_patients[!clinical_patients %in% colnames(brca_exp_mat)]

# Eliminar las muestras sobrantes
clinical_info_filtered <- clinical_info[clinical_info$sample %in% colnames(brca_exp_mat), ]

# Comprobar que tienen la misma dimension y estan en el mismo orden
all(colnames(brca_exp_mat)==clinical_info_filtered$sample)

# Comprobar primeras filas y columnas
head(brca_exp_mat[, 1:5])

# Comprobamos las dimensiones de las matrices y sus clases
dim(clinical_info_filtered) # Filas: muestras; columnas: variables clinicas
class(clinical_info_filtered)

dim(brca_exp_mat)
class(brca_exp_mat)

# Ver un resumen de ambos datasets
summary(brca_exp_mat[, 1:5])
summary(clinical_info_filtered[, 1:5])

# Ver quetipo de muestras hay
table(clinical_info_filtered$er_status_by_ihc)
# Ver cuantos valores NA hay
sum(is.na(clinical_info_filtered$er_status_by_ihc))
sum(is.na(clinical_info_filtered$age_at_diagnosis))

# Eliminar muestras que tengan tumor_status == NA
clinical_info_clean <- clinical_info_filtered[!is.na(clinical_info_filtered$er_status_by_ihc), ]

# Ver si quedan valores NA
sum(is.na(clinical_info_clean$er_status_by_ihc))
sum(is.na(clinical_info_clean$age_at_diagnosis))

# Obtener el ID de las muestras en ambos datasets
samples_exp <- colnames(brca_exp_mat)
samples_clinic <- clinical_info_clean$sample

# Buscar muestras presentes en ambos datasets
common_samples <- intersect(samples_exp, samples_clinic)

cat("Muestras de la matriz de expresiC3n:", length(samples_exp), "\n")
cat("Muestras de los datos clC-nicos:", length(samples_clinic), "\n")
cat("Muestras comunes:", length(common_samples), "\n")

# Crear subset de las muestras comunes
brca_exp_clean <- brca_exp_mat[, common_samples]

# Verificar que las dimensiones cuadran
ncol(brca_exp_clean)
nrow(clinical_info_clean)

# Verificar que el orden es el mismo
all(colnames(brca_exp_clean) == clinical_info_clean$sample)

table(clinical_info_clean$er_status_by_ihc)
table(clinical_info_clean$age_at_diagnosis)

### Normalizacion: Limpiar conteo quitando genes poco expresados
# Convertimos el data frame en una matriz
exp_matrix <- as.matrix(brca_exp_clean)
class(exp_matrix) <- "numeric" # Forzar a que la clase sea numerica para poder crear la 
#matriz de expresion

# Creamos una matriz de expresiC3n
de.exp <- DGEList(counts=exp_matrix) # Este objeto espera una matriz. 

# Nos quedamos con genes que estC)n expresados en, al menos, 224 muestras
keep.exprs <- rowSums(cpm(de.exp)>1)>=224 # Vector booleano
# Vemos cuantos genes estan poco expresados
table(keep.exprs)

# Creamos un subset de las filas (genes) que tienen expresion en mas de 224 muestras
de.exp <- de.exp[keep.exprs,] # Objeto DGEList
# Vemos cuantos genes quedan despues de la filtracion
dim(de.exp)

# Normalizacion TMM
limma <- calcNormFactors(de.exp, method = "TMM")
# Extraemos la matriz de conteos normalizada
tmmexp <- cpm(limma, log = T, prior.count = 3) # counts*1000000/(libs.size*fact)

# Guardar tabla 
write.table(tmmexp, "data_norm.tsv", col.names = T, row.names = T, sep = "\t", quote = F)

# Hacemos un histograma de la expresion de los genes normalizados
hist(tmmexp, 1000, main = paste("Histograma de la matriz de conteos normalizada "))

# Graficar siempre las mismas 15 muestras
set.seed(8)
samples_indices = sample(1:ncol(de.exp), 15)

# Crear matriz de conteos sin normalizar
no_tmmexp<-cpm(de.exp,log = T,prior.count = 3)

# Crear subsets de 15 muestras
no_tmmexp_subset <- no_tmmexp[, samples_indices]
tmmexp_subset <- tmmexp[, samples_indices]

# Boxplot de datos no normalizados vs. normalizados
par(mfrow=c(1,2), oma = c(0, 0, 2, 0)) # oma aC1ade margen externo 
boxplot(no_tmmexp_subset, las=2,  main="")
title(main="A. Datos crudos",ylab="Log-cpm")
boxplot(tmmexp_subset, las=2,  main="")
title(main="B. Datos normalizados",ylab="Log-cpm")
mtext("Efecto de la normalizaciC3n TMM sobre 15 muestras aleatorias", 
      side = 3, line = 0, outer = TRUE, cex = 1.3, font = 2) # AC1ade titulo principal
par(mfrow=c(1,1))

# PCA despues de normalizacion
pca_results <- prcomp(t(tmmexp), scale. = TRUE)

# Calcular varianza explicada
var_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100

# Preparar datos para ggplot
pca_data <- data.frame(
  PC1 = pca_results$x[,1],
  PC2 = pca_results$x[,2],
  sample = rownames(pca_results$x)
)

# AC1adir metadatos clC-nicos
pca_data <- merge(pca_data, clinical_info_filtered[,c("sample", "tumor_status", "gender", 
                                                      "age_at_diagnosis", "bcr_patient_barcode", "er_status_by_ihc")], 
                  by.x = "sample", by.y = "sample")

# GrC!fico PCA
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = er_status_by_ihc)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA - ER+ vs ER-",
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)"),
       color = "Tipo de muestra") +
  theme_bw() +
  theme(legend.position = "right")

# Mostrar grafico
print(pca_plot)

# Crear matriz de diseC1o
design <- model.matrix(~0+clinical_info_clean$er_status_by_ihc + 
                         clinical_info_clean$age_at_diagnosis)

# Limpiar los nombres de las columnas
colnames(design) <- gsub("clinical_info_clean$er_status_by_ihc", "",
                         colnames(design), fixed = T) 
colnames(design) <- gsub("clinical_info_clean$age_at_diagnosis", "age",
                         colnames(design), fixed = T)

# Hacemos la grC!fica de voom
fit <- voom(limma, design, plot = T) 
# RegresiC3n lineal con la expresiC3n de los datos y su varianza
fit <- lmFit(fit, design)

# Contrastes
# Comparacion de tumores vs tejido normal
cont.wt <- makeContrasts(er_p_vs_n=Negative - Positive,  levels=design)
# Hacemos la regresion lineal con los contrastes
fit2 <- contrasts.fit(fit, cont.wt)

# Aproximacion bayesiana
fit2 <- eBayes(fit2)
plotSA(fit2)

# Obtener tabla de genes diferencialmente expresados
de_res<-topTable(fit2,coef = "er_p_vs_n",number = Inf) #inf para que devuelva todos los genes
table(de_res$adj.P.Val<0.05)

# Guardar tabla 
write.table(de_res, "p_values.tsv", col.names = T, row.names = T, sep = "\t", quote = F)

table(de_res$adj.P.Val<0.05 & abs(de_res$logFC)>1)

### Volcano plot
# Si se cumplen ambas condiciones, serC! significante. Si ademas se cumple que logFC > 0, 
#sera upregulated. 
de_res$Significance <- ifelse(de_res$adj.P.Val < 0.5 & abs(de_res$logFC) > 1,
                              ifelse(de_res$logFC > 0, "Upregulated", "Downregulated"), 
                              "Not Significant") 

# Guardamos el nombre de los genes en una nueva columna
de_res$ID<-rownames(de_res)

# Ordenar los datos por significancia (adj.P.Val) para seleccionar el top 10
top_genes <- de_res[order(de_res$adj.P.Val), ][1:10, ]

#Definimos los parametros del grafico
volcano_plot <- ggplot(de_res, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +  # Puntos del grC!fico
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", 
                                "Not Significant" = "grey")) +
  theme_minimal() +  # Tema limpio
  labs(title = "Volcano Plot", 
       x = "Log2 Fold Change", 
       y = "-log10 Adjusted P-Value") +
  # LC-neas de FC
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.7) +  
  # LC-nea de p-valor
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.7) +  
  theme(legend.position = "top") +
  geom_text_repel(data = top_genes, aes(label = ID), 
                  size = 3, col="black",
                  box.padding = 0.3, 
                  point.padding = 0.2, 
                  max.overlaps = 10) 

#Pintamos el plot
volcano_plot

### heatmap
# Seleccionar los 10 genes mC!s upregulados y los 10 genes mC!s downregulados
top_upregulated <- de_res[de_res$Significance == "Upregulated", ]
top_upregulated <- top_upregulated[order(-top_upregulated$logFC), ][1:10, ]

top_downregulated <- de_res[de_res$Significance == "Downregulated", ]
top_downregulated <- top_downregulated[order(top_downregulated$logFC), ][1:10, ]

# Combinar los 20 genes seleccionados 
top_genes <- rbind(top_upregulated, top_downregulated)

# Extraer la matriz de expresiC3n para estos genes
heatmap_data <- tmmexp[rownames(tmmexp) %in% top_genes$ID, ]

# Escalar los datos por filas (genes) para mejor visualizaciC3n
scaled_data <- t(scale(t(heatmap_data)))  # Escalado por fila

# Crear el heatmap
pheatmap(scaled_data,
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Escala de colores
         cluster_rows = TRUE,  # Agrupar genes
         cluster_cols = TRUE,  # Agrupar muestras
         show_rownames = TRUE,  # Mostrar nombres de genes
         show_colnames = TRUE,  # Mostrar nombres de muestras
         fontsize_row = 8,  # TamaC1o de la fuente para los nombres de los genes
         fontsize_col = 10,  # TamaC1o de la fuente para los nombres de las muestras
         main = "Heatmap del Top 20 de Genes Up- y Downregulados")  # TC-tulo del heatmap

# Llamamos al archivo de la base de datos GO, donde Hs es Homo Sapiens
# utilizando la ontologia molecular function
hsGO <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont = "MF", computeIC=FALSE)
# Creamos la lista de genes con la que trabajar
genes <- rownames(top_genes)
# Calcular la matriz de similitud con el metodo Wang combinado con BMA
mgeneSim(genes, semData = hsGO, measure = "Wang", combine = "BMA", verbose=FALSE)

# Medir la correlaciC3n de los 10 genes upregulados con los 10 genes downregulados
gs1 <- rownames(top_upregulated)
gs2 <- rownames(top_downregulated)
GOSemSim::clusterSim(
  gs1, gs2, semData = hsGO, measure = "Wang", combine = "BMA")

# --------------------------------
### Analisis de enriquecimiento
# --------------------------------

### GSEA
# Guardamos la lista de genes significativos
signif_genes_df <- de_res[de_res$Significance == "Downregulated" | de_res$Significance == "Upregulated", ]
# Lo convertimos en vector
signif_genes <- sign(signif_genes_df$logFC) * (-log10(signif_genes_df$P.Value))
# Le damos nombre a los genes
names(signif_genes) <- rownames(signif_genes_df)
# Ordenar en orden descendente
signif_genes <- sort(signif_genes, decreasing = TRUE)
# Comprobar el objeto
str(signif_genes)

# Analisis de enriquecimiento de genes
ego <- gseGO(geneList = signif_genes,
             OrgDb = org.Hs.eg.db,
             ont = "BP", # Biological Process
             keyType = "SYMBOL",
             minGSSize = 50, # minimo 50 genes asociados para el termino
             maxGSSize = 300, # maximo 300 genes asociados para el termino
             pvalueCutoff = 0.05,
             verbose = FALSE)
# min: valor minimo de genes con los que trabajo. idem max
head(ego, 3)
dim(ego)

goplot(ego)

## Usar base de datos KEGG
# AC1adir ENTREZ IDs a los datos
id_conversion <- bitr(
  geneID = signif_genes_df$ID,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
# Fusionar con los datos originales
genes_entrez <- merge(signif_genes_df, id_conversion, 
                      by.x = "ID", by.y = "SYMBOL", all.x = FALSE)

# Comprobar en cuantos genes ha funcionado la conversion
cat("Exito de la conversion:", nrow(genes_entrez)/nrow(signif_genes_df)*100, "%\n")

# Lo convertimos en un vector
genes_entrez_vec <- sign(genes_entrez$logFC) * (-log10(genes_entrez$P.Value))
# Le damos nombre a los genes
names(genes_entrez_vec) <- rownames(genes_entrez)
# Ordenar en orden descendente
genes_entrez_vec <- sort(genes_entrez_vec, decreasing = TRUE)
# Comprobar el objeto
str(genes_entrez_vec)

# ------
### Probar
entrezIds <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(signif_genes_df), 
                                   columns = "ENTREZID", keytype = "SYMBOL")

entrezIds <- entrezIds[!is.na(entrezIds$ENTREZID),]
entrezIds <- entrezIds[!duplicated(entrezIds$SYMBOL),]
rownames(entrezIds) <- entrezIds$SYMBOL
###

# GSEA - KEGG vias metabolicas
kegg1 <- gseKEGG(geneList = entrezIds, 
                 organism = "hsa",
                 minGSSize = 10,
                 maxGSSize = 300,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 verbose = FALSE)

head(kegg1, 3)

## ORA
# CreaciC3n del vector de carC!cteres
signif_genes_names <- names(signif_genes)

# Enriquecimiento en GO, utilizando el mC)todo FDR
ego2 <- enrichGO(gene = signif_genes_names,
                 universe = rownames(de_res),
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 keyType = "SYMBOL",
                 pAdjustMethod = "BH", # Benjamin-Hochberg
                 pvalueCutoff = 0.001, # Reducir el numero de genes asociados
                 qvalueCutoff = 0.05,
                 readable = TRUE)
head(ego2, 3)
dim(ego2)

goplot(ego2)

# ORA - KEGG vias metabolicas
kegg2 <- enrichKEGG(gene = genes_entrez$ENTREZID, 
                 organism = "hsa",
                 minGSSize = 50,
                 maxGSSize = 300,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")

head(kegg2)

dev.new()
dotplot(kegg2, showCategory=20) + ggtitle("KEGG Pathways ORA")

# Plot a cnet plot now
cnetplot(kegg2R, showCategory = 10)

sessionInfo()