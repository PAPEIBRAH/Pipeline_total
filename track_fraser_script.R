############################################################
# Charger le fichier FRASER
############################################################

data <- read.delim(
  "~/Bureau/FRASER_output/MMG_GTF50/all_results_fraser_56_samples.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

############################################################
# Vérifier les données
############################################################

cat("Nombre total de lignes :", nrow(data), "\n")
cat("Nombre de colonnes :", ncol(data), "\n\n")

str(data)

############################################################
# Filtrer les événements significatifs
############################################################

data_sig <- subset(data, padjust < 0.05)

############################################################
# Résumé
############################################################

cat("\nNombre d'événements avec padjust < 0.05 :", nrow(data_sig), "\n")

############################################################
# Afficher les premières lignes
############################################################

head(data_sig)

############################################################
# Sauvegarder les résultats
############################################################

write.table(
  data_sig,
  file = "~/Bureau/FRASER_output/MMG_GTF50/all_results_fraser_56_samples_padjust_lt_0.05.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nFichier sauvegardé avec succès.\n")







############################################################
# Charger le fichier filtré
############################################################

data <- read.delim(
  "~/Bureau/FRASER_output/MMG_GTF50/all_results_fraser_56_samples_padjust_lt_0.05.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


############################################################
# Nom UCSC sans espaces
############################################################

data$name <- paste0(
  data$sampleID,
  "_",
  "_PSI=",
  round(data$psiValue,3),
  "_dPSI=",
  round(data$deltaPsi,3),
  "_padj=",
  signif(data$padjust,3),
  "_count=",
  data$counts
)


############################################################
# Couleur
############################################################

data$itemRgb <- ifelse(
  data$deltaPsi > 0,
  "255,0,0",
  "0,0,255"
)


############################################################
# BED12
############################################################

bed12 <- data.frame(
  
  chrom = data$seqnames,
  
  chromStart = data$start - 1,
  
  chromEnd = data$end,
  
  name = data$name,
  
  score = pmin(1000, round(abs(data$deltaPsi)*1000)),
  
  strand = ".", 
  
  thickStart = data$start - 1,
  
  thickEnd = data$end,
  
  itemRgb = data$itemRgb,
  
  blockCount = 1,
  
  blockSizes = data$width,
  
  blockStarts = 0
  
)


############################################################
# Sauvegarde
############################################################

write.table(
  bed12,
  file="~/Bureau/FRASER_output/MMG_GTF50/FRASER_UCSC_track_clean.bed",
  sep="\t",
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE
)



































