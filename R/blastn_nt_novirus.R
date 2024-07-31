#' Trim Contamination from Nucleotide Sequences
#'
#' This function trims contamination from nucleotide sequences using BLAST.
#'
#' @param input_file Path to the input nucleotide sequences file.
#' @param db Path to the BLAST database.
#' @param out_tsv Path to the output TSV file.
#' @param threads Number of threads to use.
#' @param evalue E-value threshold for BLAST.
#' @param output_nt_fasta Path to the output nucleotide FASTA file.
#' @return Path to the output nucleotide FASTA file.
#' @importFrom data.table fread
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr str_detect

trim_contamination <- function(input_file, db, out_tsv, threads, evalue, output_nt_fasta) {
  # Call BLAST
  system(paste(
    "blastn -query", input_file,
    "-db", db,
    "-out", out_tsv, "-evalue", evalue,
    '-outfmt "6 qacc qlen sseqid slen pident length qstart qend sstart send evalue bitscore qcovs"',
    "-num_threads", threads,
    "-max_target_seqs 5"
  ))

  # Read BLAST output
  blastn_table <- fread(out_tsv,
    header = FALSE, sep = "\t",
    col.names = c(
      "qacc", "qlen", "NT_sseqid", "NT_slen", "NT_pident", "NT_length",
      "qstart", "qend", "NT_sstart", "NT_send", "NT_evalue", "NT_bitscore", "NT_qcovs"
    )
  )

  # Extract unique sequence IDs
  ID_list <- unique(blastn_table$qacc)
  dict_fasta <- readDNAStringSet(input_file)

  # Helper function to write sequences to a FASTA file
  write_fasta <- function(id, seq, file) {
    cat(paste0(">", id, "\n", seq, "\n"), file = file, append = TRUE)
  }

  # Process sequences and write to output
  for (seq_id in names(dict_fasta)) {
    seq <- as.character(dict_fasta[[seq_id]])
    if (seq_id %in% ID_list) {
      line <- blastn_table[blastn_table$qacc == seq_id, ]
      qstart <- min(line$qstart)
      qend <- max(line$qend)
      origin_qlen <- nchar(seq)
      now_start <- min(qstart, qend)
      now_end <- max(qstart, qend)
      left <- now_start - 1
      right <- origin_qlen - now_end
      if (left > right && left > 600) {
        left_seq <- substr(seq, 1, left)
        write_fasta(seq_id, left_seq, output_nt_fasta)
      } else if (right > left && right > 600) {
        right_seq <- substr(seq, now_end + 1, origin_qlen)
        write_fasta(seq_id, right_seq, output_nt_fasta)
      }
    } else {
      write_fasta(seq_id, seq, output_nt_fasta)
    }
  }
  return(output_nt_fasta)
}

#' Run Trim Contamination
#'
#' This function parses command line arguments and runs the trim contamination process.
#'
#' @export
#' @importFrom argparse ArgumentParser

run_trim_contamination <- function() {
  parser <- ArgumentParser(description = "Centroid based sequence clustering")
  parser$add_argument("--input", type = "character", required = TRUE, metavar = "PATH", help = "Path to nucleotide sequences")
  parser$add_argument("--out_fasta", type = "character", required = TRUE, metavar = "BASENAME", help = "Path to output file")
  parser$add_argument("--out_tsv", type = "character", metavar = "PATH", help = "")
  parser$add_argument("--threads", type = "integer", metavar = "INT", help = "")
  parser$add_argument("--evalue", type = "character", metavar = "EVALUE", help = "")
  parser$add_argument("--db", type = "character", required = TRUE, metavar = "PATH", help = "Path to the blastn database")
  args <- parser$parse_args()

  # Run the trimming process
  trim_contamination(args$input, args$db, args$out_tsv, args$threads, args$evalue, args$out_fasta)
}

# Entry point for the command line script
if (sys.nframe() == 0) {
  run_trim_contamination()
}
