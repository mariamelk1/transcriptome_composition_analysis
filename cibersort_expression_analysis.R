# cibersort_expression_analysis.R
# Autor: Mariam El Khattabi
# Descripción: Análisis de composición celular con CIBERSORTx y análisis de expresión diferencial ajustado por covariables y proporción celular

###############################       CARGAR LIBRERÍAS      ##########################
library(edgeR)
library(biomaRt)
library(limma)
library(MASS)
library(ggplot2)

###############################       CONFIGURACIÓN INICIAL     ##########################
# Define directorios genéricos
data_dir <- "data/"             # Donde están los datos (y_raw2, targets, etc.)
output_dir <- "results/"        # Carpeta para guardar resultados

# Crea carpetas si no existen
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

###############################       NORMALIZACIÓN      ##########################
# CIBERSORTx requiere expresión normalizada pero no transformada (sin log)

load(file.path(data_dir, "y_raw2.RData"))  # Objeto DGEList
dge <- calcNormFactors(y_raw2)
cpm_matrix <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)

##############################  CAMBIAR IDs ENSEMBL A SÍMBOLOS DE GEN #####################
# Limpieza de IDs y anotación con biomaRt
expr_CPM <- cpm_matrix
expr_CPM_out <- cbind(Gene = rownames(expr_CPM), expr_CPM)
expr_CPM_out$Gene <- sub("\\..*$", "", expr_CPM_out$Gene)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map_biomart <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"), 
  values = expr_CPM_out$Gene,
  mart = ensembl
)

expr_annotated <- merge(gene_map_biomart, expr_CPM_out,
                        by.x = "ensembl_gene_id", by.y = "Gene")
colnames(expr_annotated)[colnames(expr_annotated) == "hgnc_symbol"] <- "Gene"
expr_annotated$ensembl_gene_id <- NULL

write.table(expr_annotated, file = file.path(output_dir, "expr_cibersort_ready.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

##############################     AGRUPAR COMPOSICIÓN CELULAR     ######################
# Leer resultados de CIBERSORTx y agrupar en 6 tipos celulares

composicion <- read.table(file.path(data_dir, "CIBERSORTx_Results.csv"),
                          header = TRUE, sep = ",", fill = TRUE, quote = "\"", row.names = 1)

grupo_mapeo <- list(
  Bcell = c("B.cells.naive", "B.cells.memory", "Plasma.cells"),
  NK = c("NK.cells.resting", "NK.cells.activated"),
  CD4Tcells = c("T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated", 
                "T.cells.follicular.helper", "T.cells.regulatory..Tregs.", "T.cells.gamma.delta"),
  CD8Tcells = c("T.cells.CD8"),
  Monocytes = c("Monocytes", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", 
                "Dendritic.cells.resting", "Dendritic.cells.activated"),
  Granulocytes = c("Mast.cells.resting", "Mast.cells.activated", "Eosinophils", "Neutrophils")
)

suma_por_grupo <- data.frame(row.names = rownames(composicion))
for (grupo in names(grupo_mapeo)) {
  tipos <- grupo_mapeo[[grupo]]
  tipos_existentes <- tipos[tipos %in% colnames(composicion)]
  suma_por_grupo[[grupo]] <- rowSums(composicion[, tipos_existentes, drop = FALSE])
}
write.csv(suma_por_grupo, file = file.path(output_dir, "composicion_agrupada.csv"))

##############################   AÑADIR COMPOSICIÓN A METADATOS   ##########################
load(file.path(data_dir, "targets_filtered2.RData"))  # contiene targets_filtered2
composicion_ordenada <- suma_por_grupo[targets_filtered2$Sample_Name, ]
targets_con_composicion <- cbind(targets_filtered2, composicion_ordenada)
save(targets_con_composicion, file = file.path(output_dir, "targets_composicion_celular.RData"))

#########################       SELECCIÓN DE MODELO (AIC)       ###########################
load(file.path(output_dir, "targets_composicion_celular.RData"))
targets_con_composicion$Diagnosis <- factor(targets_con_composicion$Diagnosis)

modelo <- glm(Diagnosis ~ Pais + RNAGender + Age + BATCH + CD4Tcells + Bcell + NK + CD8Tcells + Monocytes + Granulocytes,
              data = targets_con_composicion, family = binomial)
stepwise_result <- capture.output(stepAIC(modelo))
writeLines(stepwise_result, file.path(output_dir, "AKAIKE_step.txt"))

#################     ANÁLISIS DE EXPRESIÓN DIFERENCIAL     #########################

load(file.path(data_dir, "y_raw2.RData"))  # se vuelve a usar y_raw2

modelo.def <- model.matrix(~ Diagnosis + Age + Pais + BATCH + CD4Tcells + NK + CD8Tcells,
                           data = targets_con_composicion)

dge <- calcNormFactors(y_raw2)  
v <- voom(dge, modelo.def, plot = TRUE)
fit <- lmFit(v, modelo.def)
fit <- eBayes(fit)

resultados <- topTable(fit, coef = "DiagnosisRA", number = Inf, adjust.method = "BH", p.value = 0.05)
write.csv(resultados, file = file.path(output_dir, "genes_significativos.csv"))

#################     VOLCANO PLOT      #########################
resultados$significativo <- ifelse(resultados$adj.P.Val < 0.05, "Sí", "No")

ggplot(resultados, aes(x = logFC, y = -log10(P.Value), color = significativo)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("No" = "grey", "Sí" = "red")) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-valor") +
  theme_minimal() +
  ggtitle("Volcano Plot (p-valor ajustado)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")
