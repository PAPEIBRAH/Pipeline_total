#!/usr/bin/env Rscript

library(stringr)

# ============================
# Paramètres
# ============================

input_dir <- "/home/pape/NAS_1/Pape/Results_rseqc_56_echantillons_gtf50_20-07-27/jonction_rseqc"

output_bed <- "/home/pape/Bureau/RSeQC_novel_junctions_track.bed"

# Orange saumon
COLOR <- "255,121,77"

# Seuil nombre de reads
MIN_COUNT <- 10


# ============================
# Liste des fichiers RSeQC
# ============================

files <- list.files(
  input_dir,
  pattern = "\\.junction\\.xls$",
  full.names = TRUE
)


# ============================
# Création BED
# ============================

con <- file(output_bed, open = "w")


# Header IGV
writeLines(
  'track name="RSeQC_Novel_Junctions" description="Novel splice junctions (read_count > 10)" visibility=pack itemRgb="On"',
  con
)


total <- 0


# ============================
# Boucle sur les échantillons
# ============================

for (f in sort(files)) {
  
  
  sample <- basename(f) |>
    str_replace("\\.junction\\.xls$", "")
  
  
  message("Traitement : ", sample)
  
  
  # Lecture RSeQC
  junction <- read.table(
    f,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  
  
  # Nettoyage annotation
  junction$annotation <- trimws(junction$annotation)
  
  
  # Nettoyage noms colonnes
  colnames(junction) <- trimws(colnames(junction))
  
  
  # Filtre novel + profondeur
  junction_filt <- junction[
    junction$annotation %in% c("complete_novel", "partial_novel") &
      junction$read_count > MIN_COUNT,
  ]
  
  
  if (nrow(junction_filt) == 0) {
    next
  }
  
  
  total <- total + nrow(junction_filt)
  
  
  
  # ============================
  # Ecriture BED12
  # ============================
  
  for (i in seq_len(nrow(junction_filt))) {
    
    
    chrom <- junction_filt$chrom[i]
    
    
    start <- as.integer(
      junction_filt$intron_st.0.based.[i]
    )
    
    
    end <- as.integer(
      junction_filt$intron_end.1.based.[i]
    )
    
    
    count <- as.integer(
      junction_filt$read_count[i]
    )
    
    
    annotation <- junction_filt$annotation[i]
    
    
    # Nom affiché IGV
    name <- paste0(
      sample,
      "_",
      annotation,
      "_count=",
      count
    )
    
    
    # Score BED
    score <- min(count,1000)
    
    
    bed_line <- paste(
      chrom,
      start,
      end,
      name,
      score,
      ".",
      start,
      end,
      COLOR,
      1,
      end-start,
      0,
      sep="\t"
    )
    
    
    writeLines(
      bed_line,
      con
    )
    
  }
  
}


close(con)


cat(
  "Nombre total de jonctions exportées :",
  total,
  "\n"
)


cat(
  "BED créé :",
  output_bed,
  "\n"
)
