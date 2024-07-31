#' Abundance Correlation Clustering Result Filtering
#'
#' This function filters the clustering results based on abundance correlation.
#' @param cor_threshold Numeric value for correlation threshold.
#' @param sample_id Character string representing the sample ID.
#' @param tpm_threshold Numeric value for TPM threshold.
#' @param rdrp_multi_threshold Numeric value for RdRp multi threshold.
#' @param non_rdrp_multi_threshold Numeric value for non-RdRp multi threshold.
#' @export
#' @importFrom data.table fread
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer
#' @importFrom openxlsx write.xlsx
#'
abundance_cor_cluster_filter <- function(cor_threshold, sample_id, tpm_threshold, rdrp_multi_threshold = 100, non_rdrp_multi_threshold = 20) {
  # Load necessary libraries
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(data.table)

  cor_threshold <- as.numeric(cor_threshold)
  tpm_threshold <- as.numeric(tpm_threshold)
  rdrp_multi_threshold <- as.numeric(rdrp_multi_threshold)
  non_rdrp_multi_threshold <- as.numeric(non_rdrp_multi_threshold)

  # Read correlation p-value data
  p_data <- fread("cor.p.1.csv", header = FALSE, data.table = FALSE)
  colnames(p_data) <- p_data[1, ]
  rownames(p_data) <- p_data[, 1]
  p_data <- p_data[-1, -1]

  # Read correlation r-value data
  r_data <- fread("cor.r.1.csv", header = FALSE, data.table = FALSE)
  colnames(r_data) <- r_data[1, ]
  rownames(r_data) <- r_data[, 1]
  r_data <- r_data[-1, -1]

  # Convert p-values and r-values to numeric
  contigs <- colnames(p_data)
  p_data <- as.data.frame(lapply(p_data, as.numeric))
  r_data <- as.data.frame(lapply(r_data, as.numeric))
  p_data <- ifelse(p_data < 0.05, 1, 0)
  r_data <- ifelse(r_data > cor_threshold, 1, 0)

  data <- p_data * r_data
  data <- apply(data, 2, sum)
  index <- which(data >= 1)

  if (length(index) > 1) {
    # Read and process RSEM data
    rsem <- fread("RSEM_rdrp.csv", header = FALSE, data.table = FALSE)
    colnames(rsem) <- rsem[1, ]
    rsem <- rsem[-1, ]
    colnames(rsem)[1] <- "Name"
    cor_contigs_TPM <- as.data.frame(t(rsem[which(rsem$Name %in% contigs[index]), ]))
    colname <- as.data.frame(t(cor_contigs_TPM))$Name
    colnames(cor_contigs_TPM) <- colname
    cor_contigs_TPM <- cor_contigs_TPM[-1, ]
    cor_contigs_TPM <- as.data.frame(lapply(cor_contigs_TPM, as.numeric))
    cor_contigs_TPM <- ifelse(cor_contigs_TPM > 0, 1, 0)
    confidence <- as.data.frame(colSums(cor_contigs_TPM))
    colnames(confidence) <- "Frequency"
    rownames(confidence) <- colname

    # Add sequences length and cut status
    seq_length <- fread("re.fasta_length.txt", sep = "\t")
    seq_length <- seq_length[which(seq_length$V1 %in% rownames(confidence)), ]
    seq_length <- seq_length[!duplicated(seq_length$V1), ]
    confidence[["real_length"]] <- seq_length$V2
    confidence[["cutted"]] <- seq_length$V3

    # Add RdRp information
    RdRp_blastx <- fread(paste0(sample_id, ".megahit.fa.rdrp.tsv"), sep = "\t", header = FALSE, data.table = FALSE)
    RdRp_blastx <- RdRp_blastx[which(RdRp_blastx$V1 %in% rownames(confidence)), ]
    RdRp_blastx <- RdRp_blastx[!duplicated(RdRp_blastx$V1), ]
    confidence[["RdRp(yes/no)"]] <- "no"
    confidence[RdRp_blastx$V1, ]$`RdRp(yes/no)` <- "yes"

    # Order and add RdRp details
    RdRp_contig_Ordered <- rownames(confidence)[which(confidence$`RdRp(yes/no)` == "yes")]
    RdRp_blastx$V1 <- factor(RdRp_blastx$V1, levels = RdRp_contig_Ordered)
    RdRp_blastx <- RdRp_blastx[order(RdRp_blastx$V1), ]
    confidence$RdRp_blastx <- ""
    confidence[RdRp_contig_Ordered, ]$RdRp_blastx <- RdRp_blastx$V3
    confidence$RdRp_identity <- ""
    confidence[RdRp_contig_Ordered, ]$RdRp_identity <- RdRp_blastx$V5

    # Add NR blast information
    NR_blastx <- fread(paste0(sample_id, "_megahit_assemble_nr.tsv"), sep = "\t", header = FALSE, data.table = FALSE)
    NR_blastx$V4 <- apply(NR_blastx, 1, function(x, c1, c2) {
      gsub(paste0(x[c1], " "), "", x[c2], fixed = TRUE)
    }, "V3", "V4")
    NR_blastx <- NR_blastx[which(NR_blastx$V1 %in% rownames(confidence)), ]
    NR_blastx <- NR_blastx[!duplicated(NR_blastx$V1), ]

    NR_contig_Ordered <- rownames(confidence)[which(rownames(confidence) %in% NR_blastx$V1)]
    NR_blastx$V1 <- factor(NR_blastx$V1, levels = NR_contig_Ordered)
    NR_blastx <- NR_blastx[order(NR_blastx$V1), ]
    confidence$NR_blastx <- ""
    confidence[NR_contig_Ordered, ]$NR_blastx <- NR_blastx$V4
    confidence$NR_identity <- ""
    confidence[NR_contig_Ordered, ]$NR_identity <- NR_blastx$V5
    confidence$NR_sstart <- ""
    confidence[NR_contig_Ordered, ]$NR_sstart <- NR_blastx[, ncol(NR_blastx) - 1]
    confidence$NR_send <- ""
    confidence[NR_contig_Ordered, ]$NR_send <- NR_blastx[, ncol(NR_blastx)]

    # Add NT blast information
    NT_blast <- fread(paste0(sample_id, "_megahit_assemble_re_nt.tsv"), sep = "\t", header = FALSE, data.table = FALSE)
    NT_blast <- NT_blast[which(NT_blast$V1 %in% rownames(confidence)), ]
    NT_blast <- NT_blast[!duplicated(NT_blast$V1), ]

    NT_contig_Ordered <- rownames(confidence)[which(rownames(confidence) %in% NT_blast$V1)]
    NT_blast$V1 <- factor(NT_blast$V1, levels = NT_contig_Ordered)
    NT_blast <- NT_blast[order(NT_blast$V1), ]
    confidence$NT_blast <- ""
    confidence[NT_contig_Ordered, ]$NT_blast <- NT_blast$V4
    confidence$NT_identity <- ""
    confidence[NT_contig_Ordered, ]$NT_identity <- NT_blast$V5
    confidence$NT_sstart <- ""
    confidence[NT_contig_Ordered, ]$NT_sstart <- NT_blast[, ncol(NT_blast) - 1]
    confidence$NT_send <- ""
    confidence[NT_contig_Ordered, ]$NT_send <- NT_blast[, ncol(NT_blast)]

    # Add TPM data
    TPM <- fread("RSEM_rdrp.csv", header = FALSE, data.table = FALSE)
    colnames(TPM) <- TPM[1, ]
    TPM <- TPM[-1, ]
    TPM <- TPM[which(TPM$Name %in% rownames(confidence)), ]

    TPM_contigs_ordered <- rownames(confidence)
    TPM[1] <- factor(TPM[, 1], levels = TPM_contigs_ordered)
    TPM <- TPM[order(TPM[, 1]), ]
    confidence$TPM <- TPM[[sample_id]]

    # Process cluster data
    cluster <- fread("cor.r.1.csv", header = FALSE, data.table = FALSE)
    cluster[1, 1] <- "contig"
    colnames(cluster) <- cluster[1, ]
    cluster <- cluster[-1, ]
    rownames(cluster) <- cluster[, 1]
    cluster <- cluster[-1]
    cluster <- cluster[which(rownames(cluster) %in% rownames(confidence)), ]
    cluster <- cluster[which(colnames(cluster) %in% rownames(confidence))]
    cluster1 <- cluster

    cluster$node1 <- rownames(cluster)
    cluster <- pivot_longer(cluster, cols = -node1, names_to = "node2", values_to = "value")
    cluster$value <- abs(as.numeric(cluster$value))
    cluster <- cluster %>%
      filter(value >= cor_threshold, node1 != node2)

    # Identify clusters
    rdrp_contigs <- rownames(cluster1)
    data <- cluster1[rdrp_contigs, ][rdrp_contigs]
    row_names <- rownames(data)
    col_names <- colnames(data)
    data <- data.frame(apply(data, 2, as.numeric))
    rownames(data) <- row_names
    colnames(data) <- col_names
    rdrp_cluster <- list()
    for (i in 1:length(rdrp_contigs)) {
      rdrp_cluster[[i]] <- unique(c(rdrp_contigs[i], rdrp_contigs[which(data[rdrp_contigs[i], ] > cor_threshold)]))
    }

    result <- list()

    # Iterate through each vector in rdrp_cluster
    while (length(rdrp_cluster) > 0) {
      # Extract the first vector from rdrp_cluster
      currentVector <- unlist(rdrp_cluster[1])
      # Remove the first vector from rdrp_cluster
      rdrp_cluster <- rdrp_cluster[-1]
      # Check if currentVector overlaps with any vector in result
      overlap <- sapply(result, function(x) any(x %in% currentVector))
      # If there is overlap, merge the vectors
      if (any(overlap)) {
        currentVector <- unique(c(currentVector, unlist(result[overlap])))
        result <- result[!overlap]
      }
      # Add the merged vector to result
      result <- c(result, list(currentVector))
    }

    confidence <- confidence[unique(c(cluster$node1, cluster$node2)), ]
    cluster.label <- rownames(confidence)

    for (i in seq_along(result)) {
      # Mark all nodes in the current cluster as belonging to the current cluster
      cluster_nodes <- unique(c(result[[i]], cluster[which(cluster$node1 %in% result[[i]]), ]$node2))
      cluster.label[cluster.label %in% cluster_nodes] <- paste0("cluster", i)
    }

    label <- unique(cluster.label)
    for (i in 1:length(label)) {
      cluster.label <- gsub(paste0("^", label[i], "$"), paste0("cluster_", i), cluster.label)
    }
    confidence$cluster <- cluster.label
    confidence <- confidence[order(confidence$cluster), ]

    confidence1 <- confidence[confidence$`RdRp(yes/no)` == "yes", ]
    confidence1$`RdRp(yes/no)` <- paste0(confidence1$`RdRp(yes/no)`, "_", confidence1$cluster)
    confidence1 <- confidence1[!duplicated(confidence1$`RdRp(yes/no)`), ]
    cor <- fread("cor.r.1.csv", header = FALSE, data.table = FALSE)
    cor[1, 1] <- "contig"
    colnames(cor) <- cor[1, ]
    cor <- cor[-1, ]
    rownames(cor) <- cor[, 1]
    cor <- cor[-1]

    confidence$cor <- ""
    confidence[rownames(confidence1), ]$cor <- "*"

    for (i in seq_len(nrow(confidence1))) {
      idx <- which(confidence$cluster == confidence1$cluster[i] & !rownames(confidence) %in% rownames(confidence1)[i])
      confidence[idx, "cor"] <- unlist(cor[rownames(confidence1)[i], rownames(confidence[idx, ])])
    }

    ### search clusters without RdRp
    confidence1 <- confidence[which(confidence$`RdRp(yes/no)` == "yes"), ]
    rdrp_cluster <- unique(confidence1$cluster)
    to_remove_clusters <- setdiff(unique(confidence$cluster), rdrp_cluster)
    if (length(to_remove_clusters) > 0) {
      confidence <- confidence[-which(confidence$cluster %in% to_remove_clusters), ]
    }

    write.xlsx(confidence, file = paste0(sample_id, ".pre.confidence_table.xlsx"), rowNames = TRUE)

    if (nrow(confidence) > 0) {
      ### search clusters with more than 1 non_virus(NR_identity > 30%)
      confidence1 <- confidence
      virus_to_remove <- "Virus|virus|Viruses|viruses|Phage|phage|Riboviria"
      confidence_to_removeCluster <- confidence1[confidence1$NR_blastx != "", ]
      confidence_to_removeCluster$NR_identity <- as.numeric(confidence_to_removeCluster$NR_identity)
      confidence_to_removeCluster <- confidence_to_removeCluster[confidence_to_removeCluster$NR_identity > 30, ]
      if (nrow(confidence_to_removeCluster) > 0) {
        confidence_to_removeCluster <- confidence_to_removeCluster[!grepl(virus_to_remove, confidence_to_removeCluster$NR_blastx), ]
        if (nrow(confidence_to_removeCluster) > 0) {
          confidence_to_removeCluster <- aggregate(confidence_to_removeCluster$cluster, by = list(cluster = confidence_to_removeCluster$cluster), length)
          confidence_to_removeCluster <- confidence_to_removeCluster[which(confidence_to_removeCluster$x >= 1), ]
          to_remove_clusters <- confidence_to_removeCluster$cluster
          if (length(to_remove_clusters) > 0) {
            confidence <- confidence[-which(confidence$cluster %in% to_remove_clusters), ]
          }
        }
      }
    }

    confidence1 <- confidence
    if (nrow(confidence1) > 0) {
      ### search clusters with more than 2 different RdRps
      confidence1 <- confidence1[which(confidence1$`RdRp(yes/no)` == "yes"), ]
      if (nrow(confidence1) > 0) {
        confidence1$NR_blastx <- paste0(confidence1$NR_blastx, "_", confidence1$NR_identity)
        cluster_RdRp_num <- aggregate(confidence1$cluster, by = list(cluster = confidence1$cluster), length)
        cluster_RdRp_num <- cluster_RdRp_num[which(cluster_RdRp_num$x > 1), ] # choose cluster with more than 1 rdrp_contigs
        if (nrow(cluster_RdRp_num) > 0) {
          confidence1 <- confidence1[which(confidence1$cluster %in% cluster_RdRp_num$cluster), ]
          confidence1$NR_blastx <- paste0(confidence1$NR_blastx, "_", confidence1$cluster)
          confidence1 <- confidence1[!duplicated(confidence1$NR_blastx), ]
          cluster <- aggregate(confidence1$cluster, by = list(cluster = confidence1$cluster), length)
          cluster <- cluster[which(cluster$x > 1), ]
          if (nrow(cluster) > 0) {
            to_remove_clusters <- cluster$cluster
            confidence <- confidence[-which(confidence$cluster %in% to_remove_clusters), ]
          }
        }
      }
    }

    if (nrow(confidence) > 0) {
      #### remove the cluster whose contigs are all rdrps
      to_remove_clusters <- vector()
      confidence1 <- confidence
      confidence1$`RdRp(yes/no)` <- paste0(confidence1$`RdRp(yes/no)`, "_", confidence1$cluster)
      confidence1 <- confidence1[!duplicated(confidence1$`RdRp(yes/no)`), ]
      cluster <- aggregate(confidence1$cluster, by = list(cluster = confidence1$cluster), length)
      cluster <- cluster[which(cluster$x == 1), ] # choose the cluster whose contigs are all rdrps(or no_rdrps)
      if (nrow(cluster) > 0) {
        confidence1 <- confidence[which(confidence$cluster %in% cluster$cluster), ]
        confidence1 <- confidence1[which(confidence1$`RdRp(yes/no)` == "yes"), ] # choose the cluster whose contigs are all rdrps
        if (nrow(confidence1) > 0) {
          to_remove_clusters <- unique(confidence1$cluster)
          confidence <- confidence[-which(confidence$cluster %in% to_remove_clusters), ]
        }
      }

      if (nrow(confidence) > 0) {
        ### remove the clusters which exist the contigs that they and their re_assemble contigs all meet the condition (multi < threshold value)
        contigs_name <- rownames(confidence)
        rdrp <- which(confidence$`RdRp(yes/no)` == "yes")
        non_rdrp <- setdiff(1:length(contigs_name), rdrp)
        contigs_name <- as.numeric(gsub(".*(cov|multi)_([0-9.]+).*", "\\2", contigs_name))

        contigs_name[rdrp] <- ifelse(contigs_name[rdrp] > rdrp_multi_threshold, 1, 0)
        contigs_name[non_rdrp] <- ifelse(contigs_name[non_rdrp] > non_rdrp_multi_threshold, 1, 0)
        multi_row <- which(contigs_name == 0)
        to_remove_clusters <- unique(confidence[multi_row, ]$cluster)
        if (length(to_remove_clusters) > 0) {
          confidence <- confidence[-which(confidence$cluster %in% to_remove_clusters), ]
        }
      }
    }

    if (nrow(confidence) > 0) {
      confidence1 <- confidence[confidence$`RdRp(yes/no)` == "yes", ]
      confidence1$`RdRp(yes/no)` <- paste0(confidence1$`RdRp(yes/no)`, "_", confidence1$cluster)
      confidence1 <- confidence1[!duplicated(confidence1$`RdRp(yes/no)`), ]
      cor <- fread("cor.r.1.csv", header = FALSE, data.table = FALSE)
      cor[1, 1] <- "contig"
      colnames(cor) <- cor[1, ]
      cor <- cor[-1, ]
      rownames(cor) <- cor[, 1]
      cor <- cor[-1]

      confidence$cor <- ""
      confidence[rownames(confidence1), ]$cor <- "*"

      for (i in seq_len(nrow(confidence1))) {
        idx <- which(confidence$cluster == confidence1$cluster[i] & !rownames(confidence) %in% rownames(confidence1)[i])
        confidence[idx, "cor"] <- unlist(cor[rownames(confidence1)[i], rownames(confidence[idx, ])])
      }
    }

    ### remove the cluster which has contigs whose TPM < threshold
    if (nrow(confidence) > 0) {
      confidence1 <- confidence
      confidence1$TPM <- as.numeric(confidence1$TPM)
      confidence1 <- confidence1[confidence1$TPM < tpm_threshold, ]
      if (nrow(confidence1) > 0) {
        remove_clusters <- unique(confidence1$cluster)
        confidence <- confidence[!confidence$cluster %in% remove_clusters, ]
      }
    }

    ### remove the cluster which has contigs whose Frequency < 3
    if (nrow(confidence) > 0) {
      confidence1 <- confidence
      confidence1$Frequency <- as.numeric(confidence1$Frequency)
      confidence1 <- confidence1[confidence1$Frequency < 3, ]
      if (nrow(confidence1) > 0) {
        remove_clusters <- unique(confidence1$cluster)
        confidence <- confidence[!confidence$cluster %in% remove_clusters, ]
      }
    }

    write.xlsx(confidence, file = paste0(sample_id, ".final.confidence_table.xlsx"), rowNames = TRUE)
    write.table(rownames(confidence), file = "Cor_contigs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
  }
}
