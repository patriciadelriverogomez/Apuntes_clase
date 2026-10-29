##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################

### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

# APARTADO 1. CARGA Y EXPLORACION DE LOS DATOS

# Cargar el archivo de entrada
datos <- read.table(
  "archivo-entrada.txt",
  header = TRUE,
  sep = "\t",
  dec = ".",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Mostrar las primeras filas
head(datos)

# Mostrar las ultimas filas
tail(datos)

# Mostrar las dimensiones de la tabla
dim(datos)

# Comprobar los nombres de las columnas
names(datos)

# Comprobar la estructura de la tabla
str(datos)

### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

# APARTADO 2. MEDIA Y DESVIACION ESTANDAR

# Seleccionar las tres réplicas de cada grupo
control <- datos[, c("C2-TPM", "C1-TPM", "C3-TPM")]
condicion <- datos[, c("1W5-TPM", "1W6-TPM", "1W5-COV")] 

# Calcular la media y la desviación estándar del Control
datos$Media_Control <- rowMeans(control, na.rm = TRUE)

datos$DE_Control <- apply(control, 1, sd, na.rm = TRUE)

# Calcular la media y la desviación estándar de la Condición
datos$Media_Condicion <- rowMeans(condicion, na.rm = TRUE)

datos$DE_Condicion <- apply(condicion, 1, sd, na.rm = TRUE)

# Comprobar los resultados
head(datos)

### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

# APARTADO 3. HISTOGRAMAS

# Definir los colores
color_control <- "steelblue"
color_condicion <- "tomato"

# Mostrar dos histogramas en una misma imagen
par(mfrow = c(1, 2))

# Histograma del control
hist(
  datos$Media_Control[is.finite(datos$Media_Control)],
  main = "Media de expresion: Control",
  xlab = "Expresion media (TPM)",
  ylab = "Frecuencia",
  col = color_control,
  border = "white"
)

# Histograma de la condicion
hist(
  datos$Media_Condicion[is.finite(datos$Media_Condicion)],
  main = "Media de expresion: Condicion",
  xlab = "Expresion media (TPM)",
  ylab = "Frecuencia",
  col = color_condicion,
  border = "white"
)

# Restaurar la configuracion grafica
par(mfrow = c(1, 1))

### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

# APARTADO 4. BOXPLOT

# Definir los colores
color_control <- "steelblue"
color_condicion <- "tomato"

# Crear los boxplots
boxplot(
  datos$Media_Control,
  datos$Media_Condicion,
  names = c("Control", "Condicion"),
  col = c(color_control, color_condicion),
  main = "Comparacion de la expresion genica",
  ylab = "Expresion media (TPM)",
  xlab = "Grupo",
  outline = TRUE
)

# Añadir la leyenda
legend(
  "topleft",
  legend = c("Control", "Condicion"),
  fill = c(color_control, color_condicion),
  title = "Grupo",
  bty = "n"
)

### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

# APARTADO 5A. NORMALIDAD EXPLORATORIA

# Evaluar la normalidad del control
shapiro_control <- shapiro.test(
  datos$Media_Control
)

# Evaluar la normalidad de la condicion
shapiro_condicion <- shapiro.test(
  datos$Media_Condicion
)

# Mostrar los resultados
print(shapiro_control)
print(shapiro_condicion)

# APARTADO 5B. TEST DE COMPARACION DE MEDIAS PARA CADA GEN (segun normalidad)

# Comparación estadística por gen mediante t-test de Welch
pvalores <- sapply(seq_len(nrow(datos)), function(i) {

  # Valores de las tres réplicas del Control
  x <- as.numeric(unlist(
    datos[i, c("C2-TPM", "C1-TPM", "C3-TPM")],
    use.names = FALSE
  ))
  
  # Valores de las tres réplicas de la Condición
  y <- as.numeric(unlist(
    datos[i, c("1W5-TPM", "1W6-TPM", "1W5-COV")],
    use.names = FALSE
  ))

  # Eliminar valores no válidos
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  # Comprobar que hay suficientes observaciones
  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }

  # Realizar el contraste estadístico
  resultado <- tryCatch(
    t.test(
      x,
      y,
      alternative = "two.sided",
      var.equal = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(resultado)) {
    NA_real_
  } else {
    resultado$p.value
  }
})

# Añadir los p-valores a la tabla
datos$p_value <- pvalores

# Corregir los p-valores por comparaciones múltiples
datos$p_ajustado <- p.adjust(datos$p_value, method = "BH")

# Mostrar los resultados
head(datos[, c(
  "Gene Name",
  "Media_Control",
  "Media_Condicion",
  "p_value",
  "p_ajustado"
)])
 
### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

# APARTADO 6. GENES DIFERENCIALMENTE EXPRESADOS

# Seleccionar genes con p-value inferior a 0.05
genes_significativos <- datos[
  !is.na(datos$p_value) & datos$p_value < 0.05,
]

# Mostrar el listado de genes significativos
genes_significativos[, c(
  "Gene Name",
  "Media_Control",
  "Media_Condicion",
  "p_value",
  "p_ajustado"
)]

# Contar los genes significativos
numero_significativos <- nrow(
  genes_significativos
)

cat(
  "Número de genes significativos (p < 0.05):",
  numero_significativos,
  "\n"
)

# Contar los genes sobreexpresados
numero_sobreexpresados <- sum(
  genes_significativos$Media_Condicion >
    genes_significativos$Media_Control,
  na.rm = TRUE
)

# Contar los genes reprimidos
numero_reprimidos <- sum(
  genes_significativos$Media_Condicion <
    genes_significativos$Media_Control,
  na.rm = TRUE
)

cat(
  "Genes sobreexpresados:",
  numero_sobreexpresados,
  "\n"
)

cat(
  "Genes reprimidos:",
  numero_reprimidos,
  "\n"
)

# Contar genes significativos tras la correccion BH
numero_ajustados <- sum(
  !is.na(datos$p_ajustado) &
    datos$p_ajustado < 0.05
)

cat(
  "Genes con p ajustado inferior a 0.05:",
  numero_ajustados,
  "\n"
)

# Guardar los resultados en un archivo de texto
write.table(
  genes_significativos,
  "genes_significativos.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

# APARTADO 7. HEATMAP DE GENES SIGNIFICATIVOS

# Instalar el paquete si no esta instalado
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# Cargar el paquete
library(pheatmap)

# Seleccionar las seis columnas de expresion
expresion <- as.matrix(
  datos[, c(
    "C2-TPM",
    "C1-TPM",
    "C3-TPM",
    "1W5-TPM",
    "1W6-TPM",
    "1W5-COV"
  )]
)

# Asegurar que los valores son numericos
storage.mode(expresion) <- "numeric"

# Asignar nombres de genes a las filas
rownames(expresion) <- make.unique(
  as.character(datos[["Gene Name"]])
)

# Asignar nombres a las muestras
colnames(expresion) <- c(
  "Control_1",
  "Control_2",
  "Control_3",
  "Condicion_1",
  "Condicion_2",
  "Condicion_3"
)

# Transformar los valores de expresion
expresion_log <- log2(expresion + 1)

# Seleccionar los genes significativos
indices_significativos <- which(
  !is.na(datos$p_value) & datos$p_value < 0.05
)

genes_heatmap <- expresion_log[
  indices_significativos,
  ,
  drop = FALSE
]

# Generar el heatmap si existen genes significativos
if (nrow(genes_heatmap) > 0) {

  pheatmap(
    genes_heatmap,
    scale = "row",
    cluster_rows = nrow(genes_heatmap) > 1,
    cluster_cols = ncol(genes_heatmap) > 1,
    show_rownames = nrow(genes_heatmap) <= 60,
    show_colnames = TRUE,
    main = "Expresión de genes significativos"
  )

  # Guardar el heatmap como imagen PNG
  pheatmap(
    genes_heatmap,
    scale = "row",
    cluster_rows = nrow(genes_heatmap) > 1,
    cluster_cols = ncol(genes_heatmap) > 1,
    show_rownames = nrow(genes_heatmap) <= 60,
    show_colnames = TRUE,
    main = "Expresión de genes significativos",
    filename = "heatmap_genes_significativos.png",
    width = 10,
    height = 8
  )

} else {
  message(
    "No hay genes con p-valor < 0.05. No se puede generar el heatmap."
  )
}