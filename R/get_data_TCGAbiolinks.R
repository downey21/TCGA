
# -*- coding: utf-8 -*-

rm(list = ls())

# BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)

# "CPTAC-3": single cell

TCGAbiolinks:::getGDCprojects()$project_id

TCGAbiolinks:::getProjectSummary("TCGA-LUAD")

query <- TCGAbiolinks::GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling"
)

head(TCGAbiolinks::getResults(query))
colnames(TCGAbiolinks::getResults(query))

# 1. "id"                        # Unique file ID (UUID), e.g., "53f2835e-e13a-4fb6-90a5-448a1a726249"
# 2. "data_format"               # File format, e.g., "TSV", "TXT"
# 3. "cases"                     # Associated case ID, e.g., "TCGA-64-1678-01A-01R-0946-07"
# 4. "access"                    # Data access level, e.g., "open" (public) or "controlled" (restricted)
# 5. "file_name"                 # Name of the file, e.g., "b1d9364c-d703-4884-b96d-20d8084040a8.rna_seq.star_splice_junctions.tsv.gz"
# 6. "submitter_id"              # Submitter-provided file ID, e.g., "e9aded74-2b22-4dab-a583-96d50c30be68"
# 7. "data_category"             # Data category, e.g., "Transcriptome Profiling"
# 8. "type"                      # Simplified file type, e.g., "gene_expression" or "mirna_expression"
# 9. "platform"                  # Platform used for data generation, e.g., "Illumina"
# 10. "file_size"                # File size in bytes, e.g., 2032798
# 11. "created_datetime"         # Date and time when the file was created, e.g., "2021-12-13T19:43:50.041252-06:00"
# 12. "md5sum"                   # MD5 checksum to verify file integrity, e.g., "b05d932361f3b6309c457c5076992a04"
# 13. "updated_datetime"         # Date and time when the file was last updated, e.g., "2024-07-30T12:17:53.416802-05:00"
# 14. "file_id"                  # File ID (same as "id"), e.g., "53f2835e-e13a-4fb6-90a5-448a1a726249"
# 15. "data_type"                # Data type (important for GDCquery), e.g., "Gene Expression Quantification"
# 16. "state"                    # File release state, e.g., "released"
# 17. "experimental_strategy"    # Experimental strategy, e.g., "RNA-Seq", "miRNA-Seq"
# 18. "version"                  # Data version number, e.g., "1"
# 19. "data_release"             # GDC data release version, e.g., "32.0 - 42.0"
# 20. "project"                  # Project ID, e.g., "TCGA-LUAD"
# 21. "analysis_id"              # Unique ID for the analysis pipeline run, e.g., "24e654d9-146a-4531-99ed-6aed2eb79845"
# 22. "analysis_state"           # Analysis state, e.g., "released"
# 23. "analysis_submitter_id"    # Submitter-provided ID for the analysis, e.g., "b1d9364c-d703-4884-b96d-20d8084040a8_star__counts"
# 24. "analysis_workflow_link"   # Link to the analysis workflow (e.g., Docker URL), e.g., "quay.io/ncigdc"
# 25. "analysis_workflow_type"   # Type of analysis workflow, e.g., "STAR - Counts"
# 26. "analysis_workflow_version" # Version of the workflow used, e.g., "eecca6e2e475aab335bd7c365025ca0d0daa144b"
# 27. "sample_type"              # Sample type, e.g., "Primary Tumor", "Solid Tissue Normal"
# 28. "is_ffpe"                  # Whether the sample is FFPE (Formalin-Fixed Paraffin-Embedded), TRUE/FALSE, or NA
# 29. "cases.submitter_id"       # Submitter-provided case (patient) ID, e.g., "TCGA-64-1678"
# 30. "sample.submitter_id"      # Submitter-provided sample ID, e.g., "TCGA-64-1678-01A"

# TCGA-LUAD (Lung Adenocarcinoma) RNA-seq expression
query <- TCGAbiolinks::GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    access = "open",
    barcode = c(
    "TCGA-44-2665-01A",
    "TCGA-64-1678-01A",
    "TCGA-97-7937-01A",
    "TCGA-05-4244-01A",
    "TCGA-05-4249-01A"
    )
)

head(TCGAbiolinks::getResults(query))
colnames(TCGAbiolinks::getResults(query))

# platform         # Filter by platform used for data generation, e.g., "IlluminaHiSeq_RNASeq"
# barcode          # Filter by specific sample barcodes (submitter IDs)
# sample.type      # Filter by sample type, e.g., "Primary Tumor", "Solid Tissue Normal"

# Download
TCGAbiolinks::GDCdownload(
    query,
    directory = "./data",
    files.per.chunk = 1,
    method = "api"
)

# Read
data <-
    TCGAbiolinks::GDCprepare(
        query,
        save = TRUE,
        save.filename = "./data/TCGA_LUAD_example.RData",
        directory = "./data"
    )

# mutant_variant_classification
# Which types of mutations should be included if using mutation data.
# Default important variants: "Frame_Shift_Del", "Missense_Mutation", "Nonsense_Mutation", etc.
# Only relevant if you are working with somatic mutation data.

# 데이터 정보 확인
data
head(colnames(data)) # Sample ID
head(rownames(data)) # Gene ID
SummarizedExperiment::assayNames(data) # expression matrix name

# Gene expression matrix (assay)
SummarizedExperiment::assay(data, "unstranded")[1:5, 1:6]

# Clinical information (colData)
SummarizedExperiment::colData(data)

# Feature (gene) information (rowData)
SummarizedExperiment::rowData(data)

# Data from linkedOmics
TCGA_LUAD_protein <- TCGAbiolinks::getLinkedOmicsData(
    project = "TCGA-LUAD",
    dataset = "RPPA (Gene Level)"
)

head(TCGA_LUAD_protein)

help(getLinkedOmicsData)
