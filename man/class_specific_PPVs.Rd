% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.r
\name{class_specific_PPVs}
\alias{class_specific_PPVs}
\title{Class specific positive predictive values (PPV)}
\usage{
class_specific_PPVs(
  fasta.file_path,
  taxonomy.file_path,
  temp.directory,
  count.reads_training,
  read_length,
  simulation_runs
)
}
\arguments{
\item{fasta.file_path}{FASTA file path (e.g. FASTA file viral.genomic.fna downloaded from NCBI: ftp://ftp.ncbi.nih.gov/refseq/release/viral).}

\item{taxonomy.file_path}{Taxonomy file path (e.g. taxonomy file taxonomy_viruses_available.csv delivered within this package).}

\item{temp.directory}{Results directory containing accuracy results, generalised confusion matrix results and accuracy as well as loss graphics of the simulation runs.}

\item{count.reads_training}{Number of sampled viruses and artificially generated reads per virus taxonomy, e.g. order (training data).}

\item{read_length}{Read length of all artificially generated reads.}

\item{simulation_runs}{Simulation runs}
}
\description{
Step 1 of the taxa classification and estimation procedure.
This function creates artificial reads based on NCBI viral reference genomes FASTA file,
computes input features and artificial neural network (ANN) model
and finally stores the trained model and classification results in R-intern variables.
}
\examples{
# Please specify your file paths and directories!
fasta.file_path = "~/ag_bioinf/genomes/viruses_na/refseq/viral.genomic.fna" # Download from NCBI!
taxonomy.file_path = "inst/extdata/taxonomy_viruses_available.csv" # Relative file path
temp.directory = "~/ag_bioinf/research/metagenomics/temp" # Results directory
count.reads_training = 100
read_length = 150
simulation_runs = 10
\dontrun{
class_specific_PPVs ( fasta.file_path, taxonomy.file_path, temp.directory, count.reads_training, read_length )
}
}
\seealso{
\code{\link{ANN_new_sample}}, \code{\link{estimation_taxa_distribution}}
}
