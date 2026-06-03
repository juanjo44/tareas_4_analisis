rm(list = ls())
options(scipen = 999)

library(FactoMineR)
library(factoextra)
library(ggplot2)

dir.create("salidas_tarea4", showWarnings = FALSE)

coffee <- read.csv(
  "coffeData.csv",
  sep = ";",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "latin1"
)

names(coffee) <- iconv(names(coffee), from = "latin1", to = "UTF-8")
names(coffee) <- tolower(gsub("[^[:alnum:]]+", "", names(coffee)))

print(names(coffee))

rownames(coffee) <- coffee[[1]]
coffee[[1]] <- NULL

vars_esperadas <- c(
  "species", "countryorigin",
  "fragrancearoma", "flavor", "aftertaste", "saltacid",
  "mouthfeel", "balance", "bittersweet",
  "uniformcup", "cleancup", "cupperpoints",
  "qualityscore"
)

faltantes <- setdiff(vars_esperadas, names(coffee))
print(faltantes)

if (length(faltantes) > 0) {
  stop(paste("Faltan variables en la base:", paste(faltantes, collapse = ", ")))
}

datos <- coffee[, vars_esperadas]
names(datos) <- c(
  "Species", "CountryOrigin",
  "FragranceAroma", "Flavor", "Aftertaste", "SaltAcid",
  "Mouthfeel", "Balance", "BitterSweet",
  "UniformCup", "CleanCup", "CupperPoints",
  "qualityscore"
)

datos$Species <- as.factor(datos$Species)
datos$CountryOrigin <- as.factor(datos$CountryOrigin)

vars_numericas <- setdiff(names(datos), c("Species", "CountryOrigin"))
datos[vars_numericas] <- lapply(datos[vars_numericas], as.numeric)

datos <- na.omit(datos)

res_mfa <- MFA(
  datos,
  group = c(1, 1, 7, 3, 1),
  type = c("n", "n", "s", "s", "s"),
  name.group = c("Species", "CountryOrigin", "Sabor", "Presentacion", "QualityScore"),
  num.group.sup = c(1),
  graph = FALSE
)

eig <- as.data.frame(res_mfa$eig)
colnames(eig) <- c("Autovalor", "PorcentajeVar", "PorcentajeAcum")
eig$Dimension <- paste0("Dim", seq_len(nrow(eig)))
eig <- eig[, c("Dimension", "Autovalor", "PorcentajeVar", "PorcentajeAcum")]

n_dim_70 <- which(eig$PorcentajeAcum >= 70)[1]
if (is.na(n_dim_70)) n_dim_70 <- 2

group_coord <- as.data.frame(res_mfa$group$coord)
group_contrib <- as.data.frame(res_mfa$group$contrib)
group_cos2 <- as.data.frame(res_mfa$group$cos2)
group_rv <- as.data.frame(res_mfa$group$RV)

ind_coord <- as.data.frame(res_mfa$ind$coord)
ind_contrib <- as.data.frame(res_mfa$ind$contrib)
ind_cos2 <- as.data.frame(res_mfa$ind$cos2)

write.csv(eig, "salidas_tarea4/01_valores_propios.csv", row.names = FALSE)
write.csv(group_coord, "salidas_tarea4/02_coord_grupos.csv")
write.csv(group_contrib, "salidas_tarea4/03_contrib_grupos.csv")
write.csv(group_cos2, "salidas_tarea4/04_cos2_grupos.csv")
write.csv(group_rv, "salidas_tarea4/05_RV_grupos.csv")
write.csv(ind_coord, "salidas_tarea4/06_coord_individuos.csv")
write.csv(ind_contrib, "salidas_tarea4/07_contrib_individuos.csv")
write.csv(ind_cos2, "salidas_tarea4/08_cos2_individuos.csv")

p1 <- fviz_screeplot(res_mfa, addlabels = TRUE) +
  ggtitle("AFM - Valores propios")
ggsave("salidas_tarea4/01_screeplot.png", p1, width = 9, height = 6, dpi = 300)

p2 <- fviz_mfa_group(res_mfa, repel = TRUE) +
  ggtitle("AFM - Nube de grupos")
ggsave("salidas_tarea4/02_grupos.png", p2, width = 9, height = 6, dpi = 300)

p3 <- fviz_mfa_var(
  res_mfa, "quanti.var",
  col.var = "contrib",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE
) + ggtitle("AFM - Variables cuantitativas")
ggsave("salidas_tarea4/03_variables_cuantitativas.png", p3, width = 10, height = 7, dpi = 300)

p4 <- fviz_contrib(res_mfa, choice = "group", axes = 1, top = 10) +
  ggtitle("Contribuciones de grupos - Dim 1")
ggsave("salidas_tarea4/04_contrib_grupos_dim1.png", p4, width = 9, height = 6, dpi = 300)

p5 <- fviz_contrib(res_mfa, choice = "group", axes = 2, top = 10) +
  ggtitle("Contribuciones de grupos - Dim 2")
ggsave("salidas_tarea4/05_contrib_grupos_dim2.png", p5, width = 9, height = 6, dpi = 300)

p6 <- fviz_mfa_ind(
  res_mfa,
  col.ind = "cos2",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = FALSE
) + ggtitle("AFM - Individuos")
ggsave("salidas_tarea4/06_individuos.png", p6, width = 9, height = 6, dpi = 300)

p7 <- fviz_mfa_ind(
  res_mfa,
  habillage = datos$Species,
  addEllipses = TRUE,
  repel = FALSE
) + ggtitle("Individuos y variable suplementaria Species")
ggsave("salidas_tarea4/07_individuos_species.png", p7, width = 9, height = 6, dpi = 300)

p8 <- fviz_cos2(res_mfa, choice = "ind", axes = 1, top = 20) +
  ggtitle("Cos2 individuos - Dim 1")
ggsave("salidas_tarea4/08_cos2_ind_dim1.png", p8, width = 9, height = 6, dpi = 300)

p9 <- fviz_cos2(res_mfa, choice = "ind", axes = 2, top = 20) +
  ggtitle("Cos2 individuos - Dim 2")
ggsave("salidas_tarea4/09_cos2_ind_dim2.png", p9, width = 9, height = 6, dpi = 300)

coord_cluster <- as.data.frame(res_mfa$ind$coord[, 1:min(5, ncol(res_mfa$ind$coord)), drop = FALSE])

set.seed(123)
k_max <- min(10, nrow(coord_cluster) - 1)

p10 <- fviz_nbclust(coord_cluster, kmeans, method = "wss", k.max = k_max) +
  ggtitle("Método del codo sobre coordenadas del AFM")
ggsave("salidas_tarea4/10_codo_kmeans.png", p10, width = 9, height = 6, dpi = 300)

res_hcpc <- HCPC(res_mfa, nb.clust = -1, graph = FALSE)

cluster_data <- res_hcpc$data.clust
write.csv(cluster_data, "salidas_tarea4/09_cluster_individuos.csv")

png("salidas_tarea4/11_dendrograma_hcpc.png", width = 1800, height = 1200, res = 200)
plot(res_hcpc, choice = "tree")
dev.off()

p11 <- fviz_cluster(res_hcpc, repel = FALSE, geom = "point") +
  ggtitle("Clusters sobre el espacio factorial")
ggsave("salidas_tarea4/12_clusters_hcpc.png", p11, width = 9, height = 6, dpi = 300)

coord1 <- res_mfa$global.pca$ind$coord[, 1]
indice_0_100 <- (coord1 - min(coord1)) / (max(coord1) - min(coord1)) * 100

ranking <- data.frame(
  Cafe = rownames(datos),
  Indice_0_100 = round(indice_0_100, 2),
  stringsAsFactors = FALSE
)

ranking <- ranking[order(ranking$Indice_0_100, decreasing = TRUE), ]

write.csv(ranking, "salidas_tarea4/13_ranking_indice.csv", row.names = FALSE)
write.csv(head(ranking, 10), "salidas_tarea4/14_top10_indice.csv", row.names = FALSE)

sink("salidas_tarea4/00_resumen_resultados.txt")
cat("TAREA 4 - AFM - GRUPO 4\n")
cat("========================================\n")
cat("Variable suplementaria: Species\n")
cat("Grupos: 1, 1, 7, 3, 1\n")
cat("Numero de observaciones analizadas:", nrow(datos), "\n")
cat("Numero de variables analizadas:", ncol(datos), "\n\n")

cat("VALORES PROPIOS\n")
print(head(eig, 10))
cat("\nNumero de dimensiones para superar 70% acumulado:", n_dim_70, "\n")
cat("Porcentaje acumulado en ese punto:", round(eig$PorcentajeAcum[n_dim_70], 2), "%\n\n")

cat("MATRIZ RV ENTRE GRUPOS\n")
print(group_rv)
cat("\n")

cat("NUMERO DE CLUSTERS EN HCPC:\n")
print(length(unique(cluster_data$clust)))
cat("\nTAMANOS DE CLUSTER:\n")
print(table(cluster_data$clust))
cat("\nTOP 10 DEL INDICE (0-100)\n")
print(head(ranking, 10))
sink()

print(head(eig, 10))
print(group_rv)
print(table(cluster_data$clust))
print(head(ranking, 10))