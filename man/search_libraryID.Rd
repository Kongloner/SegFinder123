% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/libraryID_search.R
\name{search_libraryID}
\alias{search_libraryID}
\title{Search the library with the highest abundance for the given contig}
\usage{
search_libraryID(input_contigs, rsem_file = "RSEM.csv")
}
\arguments{
\item{input_contigs}{Name of the input contig}

\item{rsem_file}{Path to the RSEM file, default is "RSEM.csv"}
}
\value{
library_ID
}
\description{
Search the library with the highest abundance for the given contig
}
