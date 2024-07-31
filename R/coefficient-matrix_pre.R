#' Process Viruses Annotation
#'
#' This function processes annotation information by replacing empty annotations with 'Viruses' and filtering based on length.
#' @param input_file Path to the input file (without the `_megahit_assemble_nr.edited.tsv` suffix).
#' @param output_dir Directory to save the output file.
#' @param length_threshold Length threshold for filtering.
#' @export
#' @importFrom stringr str_split
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
process_viruses_annotation <- function(input_file, output_dir, length_threshold) {
  library(stringr)
  library(dplyr)
  library(magrittr)

  # Set working directory
  setwd(output_dir)

  # Read data from the input file
  data1 <- read.table(paste0(input_file, "_megahit_assemble_nr.edited.tsv"), sep = "\t", quote = "")
  data1 <- data1[, -1]
  colnames(data1) <- c("Contigs", "Length", "Accession", "Species", "Similarly", "aa_length", "E-value", "Species_annotation")

  # Split and replace empty annotation information with 'Viruses'
  split_b <- str_split(data1$Species_annotation, "\\|", simplify = TRUE)
  data1$Species_annotation <- split_b[, 2]
  data1$Species_annotation[is.na(data1$Species_annotation)] <- "Viruses"

  # Filter data
  data2 <- data1 %>% filter(Species_annotation == "Viruses")
  data_res <- data2[which(data2$Length > as.numeric(length_threshold)), ]

  # Write the result to a file
  write.table(data_res, "RSEM_pre.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = TRUE)
}
