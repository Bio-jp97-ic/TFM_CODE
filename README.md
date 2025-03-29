# TFM_CODE
Código sobre el árbol filogenetico sobre miRs y EN
# NOMBRES: Paul Jara-Ortega  
#Crear archivo fasta. Definir el nombre y las secuencias de los miRs
seq_id <- c("hsa-miR-16-5p", "hsa-miR-16-1-3p", "hsa-miR-29a-5p", 
            "hsa-miR-29a-3p", "hsa-miR-29b-1-5p", "hsa-miR-29b-3p", 
            "hsa-miR-29c-5p", "hsa-miR-29c-3p", "hsa-miR-33a-5p", 
            "hsa-miR-33a-3p", "hsa-miR-33b-5p", "hsa-miR-33b-3p", 
            "hsa-miR-34a-5p", "hsa-miR-34a-3p", "hsa-miR-34c-5p", 
            "hsa-miR-34c-3p", "hsa-miR-98-5p", "hsa-miR-101-5p", 
            "hsa-miR-101-3p", "hsa-miR-101-2-5p", "hsa-miR-124-5p", 
            "hsa-miR-124-3p", "hsa-miR-125b-5p", "hsa-miR-125b-1-3p", 
            "hsa-miR-125b-2-3p", "hsa-miR-128-3p", "hsa-miR-128-1-5p", 
            "hsa-miR-128-2-5p", "hsa-miR-139-5p", "hsa-miR-139-3p", 
            "hsa-miR-146a-5p", "hsa-miR-146a-3p", "hsa-miR-186-5p", 
            "hsa-miR-186-3p", "hsa-miR-188-3p", "hsa-miR-206", 
            "hsa-miR-132-5p", "hsa-miR-132-3p", "hsa-miR-212-3p", 
            "hsa-miR-212-5p", "hsa-miR-324-3p", "hsa-miR-330-5p", 
            "hsa-miR-330-3p", "hsa-miR-339-5p", "hsa-miR-922", 
            "hsa-miR-937-3p", "hsa-miR-937-5p", "hsa-let-7d-5p", 
            "hsa-let-7d-3p", "hsa-miR-7-5p", "hsa-miR-7-1-3p", 
            "hsa-miR-7-2-3p", "hsa-miR-22-5p", "hsa-miR-22-3p", 
            "hsa-miR-30e-5p", "hsa-miR-30e-3p", "hsa-miR-34b-5p", 
            "hsa-miR-34b-3p", "hsa-miR-126-5p", "hsa-miR-126-3p", 
            "hsa-miR-205-5p", "hsa-miR-205-3p", "hsa-miR-214-5p", 
            "hsa-miR-214-3p", "hsa-miR-221-5p", "hsa-miR-221-3p", 
            "hsa-miR-342-3p", "hsa-miR-15b-5p", "hsa-miR-15b-3p", 
            "hsa-miR-17-5p", "hsa-miR-17-3p", "hsa-miR-18a-5p", 
            "hsa-miR-18a-3p", "hsa-miR-19a-5p", "hsa-miR-19a-3p", 
            "hsa-miR-19b-1-5p", "hsa-miR-19b-3p", "hsa-miR-19b-2-5p", 
            "hsa-miR-20a-5p", "hsa-miR-20a-3p", "hsa-miR-92a-1-5p", 
            "hsa-miR-92a-3p", "hsa-miR-92a-2-5p", "hsa-miR-155-5p", 
            "hsa-miR-155-3p", "hsa-miR-223-5p", "hsa-miR-223-3p", 
            "hsa-miR-448", "hsa-miR-150-5p", "hsa-miR-150-3p")

seq_cont <- c("UAGCAGCACGUAAAUAUUGGCG", "CCAGUAUUAACUGUGCUGCUGA", 
              "ACUGAUUUCUUUUGGUGUUCAG", "UAGCACCAUCUGAAAUCGGUUA", 
              "GCUGGUUUCAUAUGGUGGUUUAGA", "UAGCACCAUUUGAAAUCAGUGUU", 
              "UGACCGAUUUCUCCUGGUGUUC", "UAGCACCAUUUGAAAUCGGUUA", 
              "GUGCAUUGUAGUUGCAUUGCA", "CAAUGUUUCCACAGUGCAUCAC", 
              "GUGCAUUGCUGUUGCAUUGC", "CAGUGCCUCGGCAGUGCAGCCC", 
              "UGGCAGUGUCUUAGCUGGUUGU", "CAAUCAGCAAGUAUACUGCCCU", 
              "AGGCAGUGUAGUUAGCUGAUUGC", "AAUCACUAACCACACGGCCAGG", 
              "UGAGGUAGUAAGUUGUAUUGUU", "CAGUUAUCACAGUGCUGAUGCU", 
              "UACAGUACUGUGAUAACUGAA", "UCGGUUAUCAUGGUACCGAUGC", 
              "CGUGUUCACAGCGGACCUUGAU", "UAAGGCACGCGGUGAAUGCCAA", 
              "UCCCUGAGACCCUAACUUGUGA", "ACGGGUUAGGCUCUUGGGAGCU", 
              "UCACAAGUCAGGCUCUUGGGAC", "UCACAGUGAACCGGUCUCUUU", 
              "CGGGGCCGUAGCACUGUCUGAGA", "GGGGGCCGAUACACUGUACGAGA", 
              "UCUACAGUGCACGUGUCUCCAGU", "UGGAGACGCGGCCCUGUUGGAGU", 
              "UGAGAACUGAAUUCCAUGGGUU", "CCUCUGAAAUUCAGUUCUUCAG", 
              "CAAAGAAUUCUCCUUUUGGGCU", "GCCCAAAGGUGAAUUUUUUGGG", 
              "CUCCCACAUGCAGGGUUUGCA", "UGGAAUGUAAGGAAGUGUGUGG", 
              "ACCGUGGCUUUCGAUUGUUACU", "UAACAGUCUACAGCCAUGGUCG", 
              "UAACAGUCUCCAGUCACGGCC", "ACCUUGGCUCUAGACUGCUUACU", 
              "CCCACUGCCCCAGGUGCUGCUGG", "UCUCUGGGCCUGUGUCUUAGGC", 
              "GCAAAGCACACGGCCUGCAGAGA", "UCCCUGUCCUCCAGGAGCUCACG", 
              "GCAGCAGAGAAUAGGACUACGUC", "AUCCGCGCUCUGACUCUCUGCC", 
              "GUGAGUCAGGGUGGGGCUGG", "AGAGGUAGUAGGUUGCAUAGUU", 
              "CUAUACGACCUGCUGCCUUUCU", "UGGAAGACUAGUGAUUUUGUUGUU", 
              "CAACAAAUCACAGUCUGCCAUA", "CAACAAAUCCCAGUCUACCUAA", 
              "AGUUCUUCAGUGGCAAGCUUUA", "AAGCUGCCAGUUGAAGAACUGU", 
              "UGUAAACAUCCUUGACUGGAAG", "CUUUCAGUCGGAUGUUUACAGC", 
              "UAGGCAGUGUCAUUAGCUGAUUG", "CAAUCACUAACUCCACUGCCAU", 
              "CAUUAUUACUUUUGGUACGCG", "UCGUACCGUGAGUAAUAAUGCG", 
              "UCCUUCAUUCCACCGGAGUCUG", "GAUUUCAGUGGAGUGAAGUUC", 
              "UGCCUGUCUACACUUGCUGUGC", "ACAGCAGGCACAGACAGGCAGU", 
              "ACCUGGCAUACAAUGUAGAUUU", "AGCUACAUUGUCUGCUGGGUUUC", 
              "UCUCACACAGAAAUCGCACCCGU", "UAGCAGCACAUCAUGGUUUACA", 
              "CGAAUCAUUAUUUGCUGCUCUA", "CAAAGUGCUUACAGUGCAGGUAG", 
              "ACUGCAGUGAAGGCACUUGUAG", "UAAGGUGCAUCUAGUGCAGAUAG", 
              "ACUGCCCUAAGUGCUCCUUCUGG", "AGUUUUGCAUAGUUGCACUACA", 
              "UGUGCAAAUCUAUGCAAAACUGA", "AGUUUUGCAGGUUUGCAUCCAGC", 
              "UGUGCAAAUCCAUGCAAAACUGA", "AGUUUUGCAGGUUUGCAUUUCA", 
              "UAAAGUGCUUAUAGUGCAGGUAG", "ACUGCAUUAUGAGCACUUAAAG", 
              "AGGUUGGGAUCGGUUGCAAUGCU", "UAUUGCACUUGUCCCGGCCUGU", 
              "GGGUGGGGAUUUGUUGCAUUAC", "UUAAUGCUAAUCGUGAUAGGGGUU", 
              "CUCCUACAUAUUAGCAUUAACA", "CGUGUAUUUGACAAGCUGAGUU", 
              "UGUCAGUUUGUCAAAUACCCCA", "UUGCAUAUGUAGGAUGUCCCAU", 
              "UCUCCCAACCCUUGUACCAGUG", "CUGGUACAGGCCUGGGGGACAG")

# Cambiar todas las "U" por "T"
seq_cont <- gsub("U", "T", seq_cont)

# Ver los resultados
seq_cont

# Crear un objeto RNAStringSet con las secuencias
sequences <- DNAStringSet(seq_cont)

# Asignar los identificadores de las secuencias
names(sequences) <- seq_id

# Escribir el archivo FASTA
writeXStringSet(sequences, "miRs_sequences.fasta")

# Cargar librerias
library(Biostrings)
library(msa)
library(ape)
library(seqinr)
library(ggplot2)
library(gplots)
library(ggdendro)
library(ggtree)
library(ggtext)

#Cargar las secuencias 
# Cargar las secuencias desde un archivo FASTA
miRs_seqs <- readDNAStringSet("miRs_sequences.fasta")

# Verificar las primeras 5 secuencias cargadas
head(miRs_seqs, 5)


# Realizar el alineamiento utilizando msa
alineamiento <- msa(miRs_seqs, method ="Muscle")
alineamiento

# Obtener la matriz de alineamiento
alignment_matrix <- as.matrix(alineamiento)
View(alignment_matrix)

# Eliminar los guiones de la matriz de alineamiento
alignment_matrix_cleaned <- gsub("-", "", alignment_matrix)

# Ver el resultado
View(alignment_matrix_cleaned)

# Definir la función calculate_conservation
calculate_conservation <- function(alignment_matrix_cleaned) {
  num_positions <- ncol(alignment_matrix_cleaned)
  conservation <- numeric(num_positions)
  
  for (i in 1:num_positions) {
    # Extraer la columna i de la matriz de alineamiento
    column <- alignment_matrix_cleaned[, i]
    
    # Filtrar solo los caracteres válidos (A, T, C, G)
    column <- column[column %in% c("A", "T", "C", "G")]
    
    # Si la columna tiene solo gaps o caracteres inválidos, asignar conservación 0
    if (length(column) == 0) {
      conservation[i] <- 0
    } else {
      
      # Calcular la frecuencia de cada base
      freq_A <- sum(column == "A") / length(column)
      freq_T <- sum(column == "T") / length(column)
      freq_C <- sum(column == "C") / length(column)
      freq_G <- sum(column == "G") / length(column)
      
      # Calcular la entropía
      entropy <- -sum(sapply(c(freq_A, freq_T, freq_C, freq_G), function(p) if (p > 0) p * log(p) else 0))
      conservation[i] <- 1 - entropy / log(4)  # log(4) base 2
    }
  }
  
  return(conservation)
}

# Calcular la conservación de residuos para el ejemplo de matriz de alineamiento
conservation <- calculate_conservation(alignment_matrix_cleaned)
conservation

# Visualizar la conservación por posición
# Convertir los datos a un data frame para ggplot2
df <- data.frame(Posicion = 1:length(conservation), Conservacion = conservation)

# Crear el gráfico con ggplot2
ggplot(df, aes(x = Posicion, y = Conservacion)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(aes(color = ifelse(Conservacion > 0.8, "Alta", ifelse(Conservacion < 0.2, "Baja", "Media"))), size = 3) +
  scale_color_manual(values = c("Alta" = "green", "Media" = "blue", "Baja" = "red")) +
  labs(title = "Conservación de Residuos por Posición",
       x = "Posición en el Alineamiento",
       y = "Conservación",
       color = "Conservación") +
  theme_minimal() +
  theme(legend.position = "top")

# Convertir a DNAbin
alignment_matrix_cleaned <- as.DNAbin(alignment_matrix_cleaned)

# Calcular la matriz de distancias evolutivas utilizando un modelo de sustitución (por ejemplo, el modelo de distancia de Hamming o el modelo de Kimura)
distances <- dist.dna(alignment_matrix_cleaned, model = "raw", pairwise.deletion = TRUE)

# Ver la matriz de distancias
print(distances)

# Convertir la matriz de distancias a un data frame para ggplot
dist_df <- as.data.frame(as.matrix(distances))

# Agregar nombres a las filas y columnas
rownames(dist_df) <- colnames(dist_df) <- seq_id

# Mejorar el gráfico del heatmap
library(gplots)

# Crear el heatmap
heatmap.2(as.matrix(dist_df), 
          trace = "none", 
          col = colorRampPalette(c("red", "blue"))(100),
          main = "Matriz de Distancias Evolutivas",
          xlab = "Secuencias",
          ylab = "Secuencias",
          dendrogram = "both", 
          margins = c(10, 10))

# Verificar si hay valores NA, NaN o Inf en la matriz de distancias
any(is.na(distances))    # Verifica si hay valores NA
any(is.nan(distances))   # Verifica si hay valores NaN
any(is.infinite(distances))  # Verifica si hay valores Inf 

# Reemplazar los valores NA y NaN con 0
distances[is.na(distances) | is.nan(distances)] <- 0

# Verificar si hay valores NA, NaN o Inf en la matriz de distancias
any(is.na(distances))    # Verifica si hay valores NA
any(is.nan(distances))   # Verifica si hay valores NaN

# Si distances es una matriz, conviértela en un objeto de clase 'dist'
dist_object <- as.dist(distances)

# Calcular el árbol filogenético con 'hclust'
arbol <- hclust(dist_object)

# Convertir el objeto hclust a un objeto 'phylo'
tree <- as.phylo(arbol)

dev.off()  # Esto cierra cualquier dispositivo gráfico activo

# Crear el grafico del arbol filogenetico
ggtree(tree) + 
  geom_tree(aes(color = branch)) +   # Colorea las ramas según el puntaje de "branch"
  scale_color_gradient(low = "blue", high = "red") +  # Define una escala de colores (ajusta según lo que desees)
  geom_tiplab(size = 1) +            # Etiquetas en las puntas (tamaño de fuente reducido)
  geom_nodelab(size = 1) +           # Etiquetas en los nodos internos (tamaño de fuente reducido)
  theme_tree() +                     # Aplica tema básico de árbol
  theme(legend.position = "none") + 
  theme(axis.title.x = element_text(size = 5),   # Ajusta tamaño de título del eje X
        axis.text.x = element_text(size = 5),    # Ajusta tamaño de las etiquetas del eje X
        axis.line.x = element_line(color = "black", size = 0.1), # Dibuja la línea del eje X
        axis.ticks.x = element_line(color = "black", linewidth = 0.5)) # Ajusta el tamaño de las marcas en el eje X

ggtree(tree, layout = "circular") + 
  geom_tree(aes(color = branch)) +   # Colorea las ramas
  scale_color_gradient(low = "blue", high = "red") +
  geom_tiplab(size = 1) +            # Etiquetas en las puntas (tamaño de fuente reducido)
  geom_nodelab(size = 1) +           # Etiquetas en los nodos internos (tamaño de fuente reducido)
  theme_tree()

# Realizar el análisis bootstrap utilizando la función boot.phylo de ape
set.seed(123)  # Para reproducibilidad
bootstrap_results <- boot.phylo(dist_object, function(x) hclust(x)$merge, B = 100)

# Opciones de graficos
ggtree(tree) + 
  geom_tree(aes(color = branch)) +   # Colorea las ramas según el puntaje de "branch"
  scale_color_gradient(low = "blue", high = "red") +  # Define una escala de colores (ajusta según lo que desees)
  geom_tiplab(size = 1) +            # Etiquetas en las puntas (tamaño de fuente reducido)
  geom_nodelab(size = 1) +           # Etiquetas en los nodos internos (tamaño de fuente reducido)
  theme_tree() +                     # Personalización general del tema
  theme(legend.position = "none") +   # Elimina la leyenda de la escala de colores
  scale_x_continuous(expand = c(0, 0)) +  # Reduce el espacio extra en el eje X (si necesario)
  coord_cartesian(clip = "off")      # Asegura que las etiquetas no se corten

ggtree(tree) + 
  geom_tree(aes(color = branch)) +   # Colorea las ramas según el puntaje de "branch"
  scale_color_gradient(low = "blue", high = "red") +  # Define una escala de colores (ajusta según lo que desees)
  geom_tiplab(size = 2) +            # Etiquetas en las puntas (tamaño de fuente reducido)
  geom_nodelab(size = 2) +           # Etiquetas en los nodos internos (tamaño de fuente reducido)
  theme_tree()

# Visualización del árbol
ggtree(tree) + 
  geom_tree(aes(color = branch)) +   # Colorea las ramas
  geom_tiplab(size = 1) +            # Etiquetas en las puntas (tamaño de fuente reducido)
  geom_nodelab(size = 1) +           # Etiquetas en los nodos internos (tamaño de fuente reducido)
  theme_tree()                       # Personalización general del tema

# Visualización del árbol
ggtree(tree, layout = "circular") + 
  geom_tree(aes(color = branch)) +   # Colorea las ramas
  geom_tiplab(size = 1) +            # Etiquetas en las puntas (tamaño de fuente reducido)
  geom_nodelab(size = 1) +           # Etiquetas en los nodos internos (tamaño de fuente reducido)
  theme_tree()                       # Personalización general del tema
