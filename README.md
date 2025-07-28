# Transcriptome Composition Analysis

Este repositorio contiene un script R reproducible utilizado durante mis prácticas en el IPBLN (CSIC) para realizar un análisis de expresión diferencial en muestras de pacientes con artritis reumatoide (RA), teniendo en cuenta la composición celular estimada a partir de datos de transcriptoma.

## 🔍 Contexto del análisis

El objetivo fue identificar genes diferencialmente expresados entre pacientes y controles, corrigiendo posibles sesgos debidos a diferencias en la composición celular de las muestras.

### 📊 Pasos principales del análisis:

1. **Control de calidad previo (no incluido aquí):**
   - Eliminación de muestras outlier y filtrado de baja calidad.

2. **Normalización de la matriz de expresión:**
   - Se normalizó la matriz de expresión en CPM (`counts per million`), sin aplicar transformación log.

3. **Conversión de IDs de gen:**
   - Se convirtieron los IDs de Ensembl a símbolos de gen (HGNC) usando `biomaRt`.

4. **Estimación de composición celular con CIBERSORTx:**
   - Esta herramienta estima la proporción relativa de 22 tipos celulares inmunes a partir de datos de expresión génica.
   - Este paso es fundamental ya que la variación en la proporción celular puede introducir sesgos en los análisis de expresión.

5. **Agrupación de tipos celulares:**
   - Para simplificar el análisis e integrarlo con datos de metilación, se agruparon los 22 tipos celulares estimados por CIBERSORT en 6 grupos biológicamente relevantes (B cells, NK, CD4 T, CD8 T, monocitos/macrófagos, granulocitos).

6. **Fusión con metadatos clínicos (`targets`):**
   - Las proporciones celulares agregadas se integraron en la tabla de diseño experimental (targets), junto con covariables como sexo, país, edad y batch.

7. **Selección de modelo con AIC:**
   - Se usó `stepAIC()` para seleccionar el modelo logístico más parsimonioso que explica el diagnóstico (RA/control) ajustando por composición celular y covariables técnicas/biológicas.

8. **Análisis de expresión diferencial (DEA):**
   - Se aplicó el pipeline `voom + limma`, que transforma los datos de RNA-seq para estabilizar la varianza dependiente de la media (corrigiendo la heterocedasticidad inherente a los datos de conteo).
   - Este método es adecuado para comparar niveles de expresión en modelos lineales, especialmente al incluir covariables como composición celular.

9. **Visualización:**
   - Se generó un Volcano Plot para visualizar los genes significativamente regulados (FDR < 0.05).

---

## 📁 Archivos generados

- `expr_cibersort_ready.txt`: matriz de expresión lista para CIBERSORT
- `composicion_agrupada.csv`: proporciones celulares resumidas
- `AKAIKE_step.txt`: resultado de selección de modelo por AIC
- `genes_significativos.csv`: tabla de genes diferencialmente expresados
- Volcano Plot (generado con ggplot2)

---

## 📦 Requisitos

- R (>= 4.0.0)
- Paquetes:
  - `edgeR`, `limma`, `biomaRt`, `ggplot2`, `MASS`

---

## 🧬 Aplicación y reutilización

Este script es adaptable para otros proyectos de análisis transcriptómico donde sea necesario controlar por la composición celular en estudios de expresión diferencial.

---

## 👩‍💻 Autora

**Mariam El Khattabi**  
Máster en Bioinformática, Universidad Claude Bernard Lyon 1  
Trabajo realizado durante prácticas en el Instituto de Parasitología y Biomedicina "López-Neyra" (CSIC), 2025.

