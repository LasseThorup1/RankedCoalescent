#' Kingman covariance matrix for n = 25
#'
#' A precomputed covariance matrix under the Kingman coalescent
#' for ranked unlabeled rooted bifurcating trees with 25 leaves.
#'
#' This matrix is used for statistical inference and hypothesis testing
#' in the RankedCoalescent package.
#'
#' @format A numeric matrix of dimension 253x253
#' @details
#' The matrix represents the covariances between all non-fixed F-matrix entries
#'
#' It was precomputed to avoid expensive computations during runtime.
#'
#' @source Generated using covariance_matrix_F
"KingmanCovMat25"

#' Kingman mean of non fixed entries for n = 25
#' #' A precomputed mean vector under the Kingman coalescent
#' for ranked unlabeled rooted bifurcating trees with 25 leaves.
#'
#' This vector is used for statistical inference and hypothesis testing
#' in the RankedCoalescent package.
#'
#' @format A vector of length 253
#' @details
#' The vector represents the expectation all non-fixed F-matrix entries
#'
#' It was precomputed to avoid expensive computations during runtime.
#'
#' @source Generated using mean_tree_DPH_nofixed
"KingmanMeanNoFixed25"

#' PMF of external branch length under the Kingman for n = 25
#' A precomputed vector of probabilities under the Kingman coalescent.
#'
#' This vector is used for statstical inference and hypothesis testing
#' in the RankedCoalescent package.
#'
#' @format A vector of length 133
#' @details The vector represents the PMF with entry i = P(E = minE - 1 + i)
#'
#' @source Generated using dist_E
"ExternalBranchLengthDist25"

#' KingmanFrechetMean25_1
#'
#' The first of two Kingman Frechet mean F-matrices with 25 lineages
#'
#' @format A matrix of dimension 24 x 24
#' @details the matrix is the first of two Frechet means in the Kingman model
#'
#' @source Generated using ViTreebi
"KingmanFrechetMean25_1"

#' KingmanFrechetMean25_2
#'
#' The second of two Kingman Frechet mean F-matrices with 25 lineages
#'
#' @format A matrix of dimension 24 x 24
#' @details the matrix is the first of two Frechet means in the Kingman model
#'
#' @source Generated using ViTreebi
"KingmanFrechetMean25_2"

#' MeanVarCov_SE
#'
#' The mean, variance and covariance of external branch length and sum of
#' non-fixed entries for n between 4-25.
#'
#' @format A data frame with columns n, MeanS, MeanE, VarS, VarE, CovSE
#' @details the data frame has 22 columns and 6 rows.
#'
#' @source Generated using covariance_matrix_SE, mean_S, and mean_E
"MeanVarCov_SE"

#' StSpM25
#'
#' The state space of the Markov embedding with n = 25 lineages.
#'
#' @format A matrix with 121,392 rows and 24 columns
#' @details Each row is a state. The state represents a column in an Fmatrix.
#'
#' @source Generated using the ranked_coalescent function
"StSpM25"

#' BCP.ssp
#'
#' The size of state spaces for the block counting process for n ranging
#' between 1 and 25.
#'
#' @format a vector of length 25.
#' @details the n'th entry is the number of states in BCP with n lineages.
#'
#' @source Hard coded
"BCPssp"

