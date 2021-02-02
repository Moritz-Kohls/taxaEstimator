\name{taxaEstimator-package}
\alias{taxaEstimator-package}
\title{This is a short tutorial with an example of how to use this package.}
\description{
The execution order of the functions of the estimation procedure is as follows:

\itemize{
\item \code{\link{class_specific_PPVs}}
\item \code{\link{ANN_new_sample}}
\item \code{\link{estimation_taxa_distribution}}

Step 1 of the taxa classification and estimation procedure.
This function creates artificial reads based on NCBI viral reference genomes FASTA file,
computes input features and artificial neural network (ANN) model
and finally stores the trained model and classification results in R-intern variables.

Step 2 of the taxa classification and estimation procedure.
This function loads an ANN model which was trained on artificial data,
computes input features of the new, adjusted sample file (FASTQ or SAM)
and stores the predicted classes of its read sequences.

Step 3 of the taxa classification and estimation procedure.
This function loads the predicted classes of the ANN model which was trained on artificial data,
computes prior as well as posterior taxa distribution estimations
and saves a graphics result file containing the estimation of the predicted classes.
}
