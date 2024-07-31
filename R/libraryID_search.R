#' Search the library with the highest abundance for the given contig
#'
#' @param input_contigs Name of the input contig
#' @param rsem_file Path to the RSEM file, default is "RSEM.csv"
#'
#' @return library_ID
#' @export
#'
search_libraryID <- function(input_contigs, rsem_file = "RSEM.csv") {
  # Load data
  data <- read.table(rsem_file, sep = ",", header = TRUE, row.names = 1)
  colnames(data) <- data[1, ]
  data <- data[-1, ]
  rowname <- rownames(data)
  colname <- colnames(data)
  data <- as.data.frame(lapply(data, as.numeric))
  rownames(data) <- rowname
  colnames(data) <- colname

  # Filter rows and columns with non-zero values
  data <- data[which(rowSums(data) > 0), ]
  data <- data[, which(colSums(data) > 0)]

  # Search library_ID
  library_ID <- colname[which.max(c(t(data[input_contigs, ])))]

  # Save result
  save(library_ID, file = paste0("library_ID/", library_ID))
}
