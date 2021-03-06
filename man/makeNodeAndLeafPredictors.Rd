% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeda-functions.R
\name{makeNodeAndLeafPredictors}
\alias{makeNodeAndLeafPredictors}
\title{Make a matrix with predictors for each leaf and node}
\usage{
makeNodeAndLeafPredictors(leafPredictors, tree)
}
\arguments{
\item{leafPredictors}{A predictor matrix for the leaves: rows are
samples, columns are leaves.}

\item{tree}{A phylogenetic tree describing the relationships
between the species/leaves.}
}
\value{
A predictor matrix for leaves and nodes together: rows are
    samples, columns are leaf/node predictors.
}
\description{
Make a matrix with one predictor for each leaf and node in the
tree, where the node predictors are the sum of the leaf predictors
descending from them.
}
