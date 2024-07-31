#' Process TPM files
#'
#' This function merges the TPM of all contigs from all libraries and generates a CSV output.
#' @param input_dir The directory containing the RdRp files.
#' @return A CSV file containing the merged TPM data.
#' @export
#' @importFrom dplyr full_join
#'
process_rdrp_files <- function(input_dir) {
  # Load necessary library
  library(dplyr)

  # List all files in the directory
  filelist <- list.files(paste0(input_dir, "/total.rdrp.megahit.fa_contigs/"))
  files <- paste(paste0(input_dir, "/total.rdrp.megahit.fa_contigs/"), filelist, sep = "")
  data <- list()

  # Read the first file to initialize finalData
  finalData <- read.table(files[1], header = T, sep = "\t")
  finalData <- finalData[c(1, which(colnames(finalData) %in% "TPM"))]
  TPMtoName <- strsplit(filelist[1], split = "_", fixed = T)[[1]][1]
  colnames(finalData)[2] <- TPMtoName
  by_name <- colnames(finalData)[1]
  finalData <- finalData[!duplicated(finalData[, 1]), ]

  # Loop through the remaining files and merge them into finalData
  for (i in 2:(length(files))) {
    data <- read.table(files[i], header = T, sep = "\t")
    data <- data[c(1, which(colnames(data) %in% "TPM"))]
    TPMtoName <- strsplit(filelist[i], split = "_", fixed = T)[[1]][1]
    colnames(data)[2] <- TPMtoName
    data <- data[!duplicated(data[, 1]), ]
    finalData <- full_join(finalData, data, by = by_name)
  }

  # Define the output file path
  output_file <- paste0(input_dir, "/RSEM.csv")

  # Write the final combined data to a CSV file
  write.csv(finalData, file = output_file, row.names = F)
}
