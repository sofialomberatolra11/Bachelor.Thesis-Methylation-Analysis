
### epiCMIT ANALYSIS ###

# Required packages
library(remotes)
library(GenomicRanges)
library(data.table)
library(DT)

# Install epiCMIT package (run once if needed)
# remotes::install_github("fridamw/epiCMIT")

### PREPARE DNA METHYLATION MATRIX ###

# Assign CpG probe names as rownames
rownames(data.annotation) <- probe_names

# Remove suffixes from duplicated CpG probe names
clean_cpgs <- sub("(_.*)$", "", rownames(data.annotation))

# Identify duplicated CpGs
duplicated_cpgs <- duplicated(clean_cpgs)

# Remove duplicated CpGs
data.annotation <- data.annotation[
  !duplicated_cpgs,
]

# Reassign cleaned CpG names
rownames(data.annotation) <- clean_cpgs[
  !duplicated_cpgs
]

# Convert to dataframe
data.annotation <- as.data.frame(data.annotation)

# Assign CpG names
rownames(data.annotation) <- data.annotation$V1

# Create DNA methylation matrix
DNAm_matrix <- data.annotation[, -1]

# Clean rownames
rownames(DNAm_matrix) <- sub(
  "_.*$",
  "",
  rownames(DNAm_matrix)
)

### CONVERT DATA TO epiCMIT FORMAT ###

DNAm.epiCMIT <- DNAm.to.epiCMIT(
  DNAm = DNAm_matrix,
  
  # Genome assembly
  DNAm.genome.assembly = "hg19",
  
  # Map probes to Illumina 450K epiCMIT CpGs
  map.DNAm.to = "Illumina.450K.epiCMIT",
  
  # Minimum number of CpGs required
  min.epiCMIT.CpGs = 800
)

DNAm.epiCMIT

### CALCULATE epiCMIT SCORES ###

epiCMIT.Illumina.par <- epiCMIT(
  DNAm.epiCMIT = DNAm.epiCMIT,
  
  # Do not export CpG annotation metadata
  return.epiCMIT.annot = FALSE,
  
  # Export results
  export.results = TRUE,
  
  # Output directory
  export.results.dir = "results/",
  
  # Output filename
  export.results.name = "epiCMIT_results"
)

head(epiCMIT.Illumina.par)

### CREATE RESULTS TABLE ###

epiCMIT.Illumina.results.par <- cbind(
  epiCMIT.Illumina.par$epiCMIT.scores,
  
  epiCMIT.CpGs =
    epiCMIT.Illumina.par$epiCMIT.run.info$epiCMIT.CpGs,
  
  epiCMIT.hyper.CpGs =
    epiCMIT.Illumina.par$epiCMIT.run.info$epiCMIT.hyper.CpGs,
  
  epiCMIT.hypo.CpGs =
    epiCMIT.Illumina.par$epiCMIT.run.info$epiCMIT.hypo.CpGs
)

head(epiCMIT.Illumina.results.par)

### DISPLAY RESULTS TABLE ###

DT::datatable(
  epiCMIT.Illumina.results.par,
  options = list(
    scrollX = TRUE,
    scrollY = TRUE
  ),
  rownames = FALSE
)

### EXPORT RESULTS ###

write.csv(
  epiCMIT.Illumina.results.par,
  file = "results/epiCMIT_results.csv",
  row.names = TRUE
)

