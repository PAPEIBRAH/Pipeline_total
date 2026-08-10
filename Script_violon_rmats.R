###############################################################
# rMATS aberrant event burden - combined violin plot
# Patient inside cohort distribution
###############################################################

library(tidyverse)
library(svglite)

###############################################################
# PATH
###############################################################

rmats_dir <- "/home/pape/Bureau/rmats_2024A181_vs_controls"

output_dir <- file.path(
  rmats_dir,
  "Figures_combined"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################
# Samples
###############################################################

patient <- "2024A181"

controls <- c(
  "2022A45",
  "2022A46",
  "2022A47",
  "2022A48",
  "2022A50",
  "2022A51",
  "2022A53",
  "2022A54",
  "2022A55",
  "2023A548",
  "2023A549",
  "2023A550",
  "2023A551",
  "2023A553",
  "2023A554",
  "2023A555",
  "2025A01",
  "2025A02",
  "2025A49",
  "2025A50"
)


samples <- c(
  patient,
  controls
)



###############################################################
# rMATS files
###############################################################

rmats_files <- c(
  SE="SE.MATS.JCEC.txt",
  A3SS="A3SS.MATS.JCEC.txt",
  A5SS="A5SS.MATS.JCEC.txt",
  MXE="MXE.MATS.JCEC.txt",
  RI="RI.MATS.JCEC.txt"
)



###############################################################
# Function
###############################################################


get_event_burden <- function(file,event_type){
  
  
  cat("\nProcessing:",event_type,"\n")
  
  
  x <- read.delim(
    file,
    check.names=FALSE,
    stringsAsFactors=FALSE
  )
  
  
  ###############################################################
  # remove duplicate columns
  ###############################################################
  
  x <- x[,!duplicated(colnames(x))]
  
  
  ###############################################################
  # significant events
  ###############################################################
  
  x <- x %>%
    filter(
      FDR < 0.05
    )
  
  
  cat(
    "Significant events:",
    nrow(x),
    "\n"
  )
  
  
  if(nrow(x)==0){
    
    return(NULL)
    
  }
  
  
  
  ###############################################################
  # Count aberrant events per sample
  ###############################################################
  
  
  burden <- rep(
    0,
    length(samples)
  )
  
  names(burden) <- samples
  
  
  
  for(i in 1:nrow(x)){
    
    
    ###############################################################
    # Patient PSI
    ###############################################################
    
    psi_patient <- suppressWarnings(
      as.numeric(
        x$IncLevel1[i]
      )
    )
    
    
    
    ###############################################################
    # Controls PSI
    ###############################################################
    
    psi_controls <- suppressWarnings(
      as.numeric(
        strsplit(
          x$IncLevel2[i],
          ","
        )[[1]]
      )
    )
    
    
    
    all_psi <- c(
      psi_patient,
      psi_controls
    )
    
    
    
    ###############################################################
    # median control reference
    ###############################################################
    
    ctrl_median <- median(
      psi_controls,
      na.rm=TRUE
    )
    
    
    
    ###############################################################
    # aberrant if delta PSI >20%
    ###############################################################
    
    aberrant <- abs(
      all_psi - ctrl_median
    ) > 0.20
    
    
    
    burden <- burden + aberrant
    
  }
  
  
  
  data.frame(
    
    sampleID=names(burden),
    
    events=as.numeric(burden),
    
    type=event_type,
    
    stringsAsFactors=FALSE
    
  )
  
  
  
}



###############################################################
# Run all event types
###############################################################


all_results <- map2_dfr(
  
  rmats_files,
  
  names(rmats_files),
  
  ~get_event_burden(
    file.path(
      rmats_dir,
      .x
    ),
    .y
  )
  
)



###############################################################
# metadata
###############################################################


all_results <- all_results %>%
  
  mutate(
    
    group="Cohort",
    
    patient_flag=
      ifelse(
        sampleID==patient,
        "2024A181",
        "Control"
      )
    
  )



###############################################################
# save table
###############################################################


write.table(
  
  all_results,
  
  file.path(
    output_dir,
    "rMATS_event_burden_all_types.tsv"
  ),
  
  sep="\t",
  
  quote=FALSE,
  
  row.names=FALSE
  
)



###############################################################
# Plot
###############################################################


p <- ggplot(
  
  all_results,
  
  aes(
    
    x=type,
    
    y=log10(events+1)
    
  )
  
)+
  
  
  geom_violin(
    
    fill="grey80",
    
    color="black",
    
    alpha=0.5,
    
    trim=FALSE
    
  )+
  
  
  geom_jitter(
    
    aes(
      color=patient_flag
    ),
    
    width=0.08,
    
    size=3,
    
    alpha=0.8
    
  )+
  
  
  scale_color_manual(
    
    values=c(
      
      Control="black",
      
      `2024A181`="red"
      
    )
    
  )+
  
  
  theme_classic(
    
    base_size=18
    
  )+
  
  
  labs(
    
    x="rMATS event type",
    
    y="log10(number of aberrant events + 1)",
    
    color=""
    
  )+
  
  
  theme(
    
    legend.position="right"
    
  )



###############################################################
# Save SVG
###############################################################


svglite(
  
  file.path(
    
    output_dir,
    
    "Figure_rMATS_event_burden_violin_combined.svg"
    
  ),
  
  width=9,
  
  height=6
  
)


print(p)


dev.off()



cat(
  "Finished successfully\n"
)
