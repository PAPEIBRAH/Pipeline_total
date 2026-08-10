###############################################################
# PCA FRASER ABERRANT SPLICING PROFILE
# Significant FRASER events only
# Patient 2024A181 highlighted
###############################################################

library(tidyverse)
library(ggrepel)
library(svglite)
library(plotly)
library(htmlwidgets)


###############################################################
# PATHS
###############################################################

fraser_file <-
  "/home/pape/Bureau/FRASER_output/MMG_GTF50/all_results_fraser_56_samples_padjust_lt_0.05.tsv"


output_dir <-
  "/home/pape/Bureau/FRASER_output/MMG_GTF50/2024A181_main_figures"


patient <- "2024A181"



###############################################################
# LOAD FRASER
###############################################################

fraser <- read.delim(
  
  fraser_file,
  
  sep="\t",
  
  stringsAsFactors=FALSE
  
)



cat(
  "Significant FRASER events:",
  nrow(fraser),
  "\n"
)



###############################################################
# CLEAN deltaPsi
###############################################################

fraser$deltaPsi <- as.numeric(
  
  gsub(
    ",",
    ".",
    fraser$deltaPsi
    
  )
  
)



###############################################################
# CREATE FEATURES PER SAMPLE
###############################################################

patient_features <- fraser %>%
  
  group_by(sampleID) %>%
  
  summarise(
    
    total_events = n(),
    
    unique_genes = n_distinct(hgncSymbol),
    
    mean_abs_deltaPsi =
      mean(
        abs(deltaPsi),
        na.rm=TRUE
      ),
    
    negative_deltaPsi =
      sum(
        deltaPsi < 0,
        na.rm=TRUE
      ),
    
    positive_deltaPsi =
      sum(
        deltaPsi > 0,
        na.rm=TRUE
      ),
    
    n_chr =
      n_distinct(seqnames),
    
    n_event_types =
      n_distinct(type)
    
  )



###############################################################
# KEEP SAMPLE NAMES
###############################################################

sample_names <- patient_features$sampleID



rownames(patient_features) <- sample_names



patient_features$sampleID <- NULL



cat(
  "Samples detected:",
  length(sample_names),
  "\n"
)



###############################################################
# REMOVE CONSTANT VARIABLES
###############################################################

feature_sd <- apply(
  
  patient_features,
  
  2,
  
  sd
  
)



patient_features <- patient_features[
  
  ,
  
  feature_sd > 0
  
]



cat(
  "Variables used:",
  ncol(patient_features),
  "\n"
)



print(
  colnames(patient_features)
)



###############################################################
# SCALE
###############################################################

pca_input <- scale(
  
  as.matrix(patient_features)
  
)



# RESTORE SAMPLE NAMES

rownames(pca_input) <- sample_names



###############################################################
# PCA
###############################################################

pca <- prcomp(
  
  pca_input,
  
  scale.=FALSE
  
)



###############################################################
# PCA DATAFRAME
###############################################################

pca_df <- data.frame(
  
  sampleID = sample_names,
  
  PC1 = pca$x[,1],
  
  PC2 = pca$x[,2],
  
  PC3 = pca$x[,3],
  
  stringsAsFactors=FALSE
  
)



###############################################################
# DEFINE GROUPS
###############################################################

pca_df$group <- "Cohort"



pca_df$group[
  pca_df$sampleID == patient
] <- "2024A181"



cat(
  "\nGroup distribution:\n"
)

print(
  table(pca_df$group)
)



cat(
  "\nPatient coordinates:\n"
)

print(
  pca_df[pca_df$sampleID==patient,]
)



###############################################################
# VARIANCE
###############################################################

var_exp <-
  
  summary(pca)$importance[2,1:3]*100



###############################################################
# PCA 2D
###############################################################

p <- ggplot(
  
  pca_df,
  
  aes(
    
    PC1,
    
    PC2
    
  )
  
)+
  
  
  geom_point(
    
    aes(
      color=group
    ),
    
    size=4
    
  )+
  
  
  geom_point(
    
    data=subset(
      pca_df,
      sampleID==patient
    ),
    
    aes(
      color=group
    ),
    
    size=7
    
  )+
  
  
  geom_text_repel(
    
    aes(
      label=sampleID
    ),
    
    size=3,
    
    max.overlaps=Inf
    
  )+
  
  
  scale_color_manual(
    
    values=c(
      
      "Cohort"="grey70",
      
      "2024A181"="red"
      
    )
    
  )+
  
  
  theme_classic(
    
    base_size=15
    
  )+
  
  
  labs(
    
    title=
      "PCA of FRASER significant splicing abnormalities",
    
    x=paste0(
      "PC1 (",
      round(var_exp[1],1),
      "%)"
    ),
    
    y=paste0(
      "PC2 (",
      round(var_exp[2],1),
      "%)"
    ),
    
    color=""
    
  )



svglite(
  
  file.path(
    
    output_dir,
    
    "PCA_FRASER_abnormal_profile_2024A181.svg"
    
  ),
  
  width=10,
  
  height=8
  
)


print(p)


dev.off()



###############################################################
# PCA 3D
###############################################################

fig <- plot_ly(
  
  pca_df,
  
  x=~PC1,
  
  y=~PC2,
  
  z=~PC3,
  
  type="scatter3d",
  
  mode="markers+text",
  
  text=~sampleID,
  
  textposition="top center",
  
  color=~group,
  
  colors=c(
    
    "Cohort"="grey70",
    
    "2024A181"="red"
    
  ),
  
  marker=list(
    
    size=7
    
  )
  
)



fig <- fig %>%
  
  layout(
    
    title=
      "3D PCA of FRASER significant splicing abnormalities",
    
    scene=list(
      
      xaxis=list(
        title=paste0(
          "PC1 (",
          round(var_exp[1],1),
          "%)"
        )
      ),
      
      yaxis=list(
        title=paste0(
          "PC2 (",
          round(var_exp[2],1),
          "%)"
        )
      ),
      
      zaxis=list(
        title=paste0(
          "PC3 (",
          round(var_exp[3],1),
          "%)"
        )
        
      )
      
    )
    
  )



saveWidget(
  
  fig,
  
  file.path(
    
    output_dir,
    
    "PCA_3D_FRASER_abnormal_profile_2024A181.html"
    
  ),
  
  selfcontained=TRUE
  
)



###############################################################
# SAVE RESULTS
###############################################################

write.table(
  
  pca_df,
  
  file.path(
    
    output_dir,
    
    "PCA_FRASER_coordinates_with_names.tsv"
    
  ),
  
  sep="\t",
  
  quote=FALSE,
  
  row.names=FALSE
  
)



write.table(
  
  patient_features,
  
  file.path(
    
    output_dir,
    
    "FRASER_patient_features.tsv"
    
  ),
  
  sep="\t",
  
  quote=FALSE
  
)



cat(
  "\nPCA FRASER completed successfully\n"
)






















features <- read.delim(
  
  "/home/pape/Bureau/FRASER_output/MMG_GTF50/2024A181_main_figures/FRASER_patient_features.tsv",
  
  sep="\t",
  
  stringsAsFactors=FALSE
  
)


head(features)

###############################################################
# EXPORT PCA FEATURES TABLE
###############################################################

patient_features_export <- patient_features %>%
  
  rownames_to_column(
    "sampleID"
  )


write.table(
  
  patient_features_export,
  
  file.path(
    
    output_dir,
    
    "PCA_features_used_for_FRASER_profile.tsv"
    
  ),
  
  sep="\t",
  
  quote=FALSE,
  
  row.names=FALSE
  
)








###############################################################
# PCA FRASER ABERRANT SPLICING PROFILE
# Significant FRASER events only
# Patient 2024A181 highlighted
###############################################################

library(tidyverse)
library(ggrepel)
library(svglite)
library(plotly)
library(htmlwidgets)


###############################################################
# PATHS
###############################################################

fraser_file <-
  "/home/pape/Bureau/FRASER_output/MMG_GTF50/all_results_fraser_56_samples_padjust_lt_0.05.tsv"


output_dir <-
  "/home/pape/Bureau/FRASER_output/MMG_GTF50/2024A181_main_figures"


patient <- "2024A181"



###############################################################
# LOAD FRASER SIGNIFICANT EVENTS
###############################################################

fraser <- read.delim(
  
  fraser_file,
  
  sep="\t",
  
  stringsAsFactors=FALSE
  
)


cat(
  "Significant FRASER events:",
  nrow(fraser),
  "\n"
)



###############################################################
# CLEAN deltaPsi
###############################################################

fraser$deltaPsi <- as.numeric(
  
  gsub(
    ",",
    ".",
    fraser$deltaPsi
    
  )
  
)



###############################################################
# CREATE FEATURES PER SAMPLE
###############################################################

patient_features <- fraser %>%
  
  group_by(sampleID) %>%
  
  summarise(
    
    total_events = n(),
    
    unique_genes = n_distinct(hgncSymbol),
    
    mean_abs_deltaPsi =
      mean(
        abs(deltaPsi),
        na.rm=TRUE
      ),
    
    negative_deltaPsi =
      sum(
        deltaPsi < 0,
        na.rm=TRUE
      ),
    
    positive_deltaPsi =
      sum(
        deltaPsi > 0,
        na.rm=TRUE
      ),
    
    n_chr =
      n_distinct(seqnames),
    
    n_event_types =
      n_distinct(type)
    
  )



###############################################################
# CHECK SAMPLE NAMES
###############################################################

cat("\nSample names detected:\n")

print(patient_features$sampleID)



###############################################################
# SAVE TABLE WITH SAMPLE IDS
###############################################################

write.table(
  
  patient_features,
  
  file.path(
    
    output_dir,
    
    "PCA_features_used_for_FRASER_profile.tsv"
    
  ),
  
  sep="\t",
  
  quote=FALSE,
  
  row.names=FALSE
  
)



###############################################################
# PREPARE PCA MATRIX
###############################################################

# garder les noms des échantillons

sample_names <- patient_features$sampleID



# enlever sampleID pour PCA

pca_matrix <- patient_features %>%
  
  select(
    -sampleID
  )



# enlever variables constantes

feature_sd <- apply(
  
  pca_matrix,
  
  2,
  
  sd
  
)


pca_matrix <- pca_matrix[
  
  ,
  
  feature_sd > 0
  
]



cat(
  "\nVariables used for PCA:",
  ncol(pca_matrix),
  "\n"
)



###############################################################
# SCALE
###############################################################

pca_input <- scale(
  
  as.matrix(pca_matrix)
  
)


rownames(pca_input) <- sample_names



###############################################################
# PCA
###############################################################

pca <- prcomp(
  
  pca_input,
  
  scale.=FALSE
  
)



###############################################################
# CREATE PCA DATAFRAME
###############################################################

pca_df <- data.frame(
  
  sampleID = sample_names,
  
  PC1 = pca$x[,1],
  
  PC2 = pca$x[,2],
  
  PC3 = pca$x[,3],
  
  stringsAsFactors=FALSE
  
)



###############################################################
# GROUPS
###############################################################

pca_df$group <- "Cohort"


pca_df$group[
  pca_df$sampleID == patient
] <- "2024A181"



cat("\nGroup distribution:\n")

print(table(pca_df$group))



###############################################################
# VARIANCE
###############################################################

var_exp <-
  
  summary(pca)$importance[2,1:3]*100



###############################################################
# PCA 2D
###############################################################

p <- ggplot(
  
  pca_df,
  
  aes(
    
    PC1,
    
    PC2
    
  )
  
)+
  
  
  geom_point(
    
    aes(
      color=group
    ),
    
    size=4
    
  )+
  
  
  geom_point(
    
    data=subset(
      pca_df,
      sampleID==patient
    ),
    
    aes(
      color=group
    ),
    
    size=7
    
  )+
  
  
  geom_text_repel(
    
    aes(
      label=sampleID
    ),
    
    size=3.5,
    
    max.overlaps=Inf
    
  )+
  
  
  scale_color_manual(
    
    values=c(
      
      "Cohort"="grey70",
      
      "2024A181"="red"
      
    )
    
  )+
  
  
  theme_classic(
    
    base_size=15
    
  )+
  
  
  labs(
    
    title="PCA of FRASER significant splicing abnormalities",
    
    x=paste0(
      "PC1 (",
      round(var_exp[1],1),
      "%)"
    ),
    
    y=paste0(
      "PC2 (",
      round(var_exp[2],1),
      "%)"
    ),
    
    color=""
    
  )



svglite(
  
  file.path(
    
    output_dir,
    
    "PCA_FRASER_profile_2024A181.svg"
    
  ),
  
  width=10,
  
  height=8
  
)


print(p)


dev.off()



###############################################################
# PCA 3D
###############################################################

fig <- plot_ly(
  
  pca_df,
  
  x=~PC1,
  
  y=~PC2,
  
  z=~PC3,
  
  type="scatter3d",
  
  mode="markers+text",
  
  text=~sampleID,
  
  textposition="top center",
  
  color=~group,
  
  colors=c(
    
    "Cohort"="grey70",
    
    "2024A181"="red"
    
  ),
  
  marker=list(
    
    size=7
    
  )
  
)



fig <- fig %>%
  
  layout(
    
    title="3D PCA of FRASER significant splicing abnormalities",
    
    scene=list(
      
      xaxis=list(
        title=paste0(
          "PC1 (",
          round(var_exp[1],1),
          "%)"
        )
      ),
      
      yaxis=list(
        title=paste0(
          "PC2 (",
          round(var_exp[2],1),
          "%)"
        )
      ),
      
      zaxis=list(
        title=paste0(
          "PC3 (",
          round(var_exp[3],1),
          "%)"
        )
        
      )
      
    )
    
  )



saveWidget(
  
  fig,
  
  file.path(
    
    output_dir,
    
    "PCA_3D_FRASER_profile_2024A181.html"
    
  ),
  
  selfcontained=TRUE
  
)



###############################################################
# SAVE PCA COORDINATES
###############################################################

write.table(
  
  pca_df,
  
  file.path(
    
    output_dir,
    
    "PCA_coordinates_FRASER_profile.tsv"
    
  ),
  
  sep="\t",
  
  quote=FALSE,
  
  row.names=FALSE
  
)



cat(
  "\nPCA FRASER completed successfully\n"
)





