###############################################################
# OUTRIDER analysis - MMG cohort
# Robust version (memory/thread controlled)
###############################################################

library(OUTRIDER)
library(rtracklayer)
library(dplyr)
library(BiocParallel)


###############################################################
# COMPUTING SETTINGS
###############################################################

workers <- 24

param <- SnowParam(
    workers = workers,
    type = "SOCK",
    progressbar = TRUE
)

register(param)

cat("Running OUTRIDER with", workers, "workers\n")


###############################################################
# PATHS
###############################################################

ctsFile <- "/home/pape/Bureau/OUTRIDER/counts_for_outrider_clean.tsv"

gtfFile <- "/home/pape/Bureau/gencode.v50.annotation.gtf"

outDir <- "/home/pape/Bureau/OUTRIDER/Results"


dir.create(
    outDir,
    recursive = TRUE,
    showWarnings = FALSE
)



###############################################################
# READ COUNTS
###############################################################

cat("Reading count matrix...\n")


ctsTable <- read.table(
    ctsFile,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    row.names = 1
)


# Remove ENSG version
rownames(ctsTable) <- sub(
    "\\..*",
    "",
    rownames(ctsTable)
)


cat(
    "Genes:",
    nrow(ctsTable),
    "\n"
)

cat(
    "Samples:",
    ncol(ctsTable),
    "\n"
)


# Check duplicated samples

if(anyDuplicated(colnames(ctsTable)) > 0){

    stop(
        "ERROR: duplicated sample names detected"
    )

}


###############################################################
# IMPORT GTF
###############################################################

cat("Reading GTF...\n")


gtf <- rtracklayer::import(gtfFile)

gtf_df <- as.data.frame(gtf)



gene_annotation <- gtf_df %>%
    filter(type == "gene") %>%
    select(
        gene_id,
        gene_name,
        width
    ) %>%
    distinct()


gene_annotation$gene_id <- sub(
    "\\..*",
    "",
    gene_annotation$gene_id
)



###############################################################
# CREATE OUTRIDER OBJECT
###############################################################

cat("Creating OUTRIDER object...\n")


ods <- OutriderDataSet(
    countData = as.matrix(ctsTable)
)


saveRDS(
    ods,
    file.path(
        outDir,
        "01_initial_ods.rds"
    )
)



###############################################################
# ADD ANNOTATION
###############################################################

mcols(ods)$basepairs <-
    gene_annotation$width[
        match(
            rownames(ods),
            gene_annotation$gene_id
        )
    ]


mcols(ods)$gene_name <-
    gene_annotation$gene_name[
        match(
            rownames(ods),
            gene_annotation$gene_id
        )
    ]



###############################################################
# FILTER EXPRESSION
###############################################################

cat("Filtering genes...\n")


ods <- filterExpression(
    ods,
    gtfFile,
    filterGenes = FALSE,
    savefpkm = TRUE
)


ods <- ods[
    mcols(ods)$passedFilter,
]


cat(
    "Genes after filtering:",
    nrow(ods),
    "\n"
)


saveRDS(
    ods,
    file.path(
        outDir,
        "02_filtered_ods.rds"
    )
)



###############################################################
# NORMALIZATION
###############################################################

cat("Estimating size factors...\n")


ods <- estimateSizeFactors(
    ods
)


saveRDS(
    ods,
    file.path(
        outDir,
        "03_normalized_ods.rds"
    )
)



###############################################################
# CONTROL CONFOUNDERS
###############################################################

cat("Controlling confounders...\n")


ods <- controlForConfounders(
    ods,
    BPPARAM = param
)


saveRDS(
    ods,
    file.path(
        outDir,
        "04_confounded_corrected_ods.rds"
    )
)



###############################################################
# P VALUES
###############################################################

cat("Computing p-values...\n")


ods <- computePvalues(
    ods,
    alternative = "two.sided",
    method = "BY",
    BPPARAM = param
)


saveRDS(
    ods,
    file.path(
        outDir,
        "05_pvalues_ods.rds"
    )
)



###############################################################
# Z SCORES
###############################################################

cat("Computing Z scores...\n")


ods <- computeZscores(
    ods
)



###############################################################
# RESULTS
###############################################################

cat("Exporting results...\n")


res <- results(
    ods,
    all = TRUE
)


write.csv(
    res,
    file.path(
        outDir,
        "OUTRIDER_results_all.csv"
    ),
    row.names = FALSE
)



saveRDS(
    ods,
    file.path(
        outDir,
        "OUTRIDER_final_object.rds"
    )
)



###############################################################
# END
###############################################################

cat("\n#################################\n")
cat("OUTRIDER FINISHED SUCCESSFULLY\n")
cat("#################################\n")
