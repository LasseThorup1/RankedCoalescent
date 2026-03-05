#' binary_matrix
#'
#' Returns a matrix of binary numbers up to a specified number of digits
#'
#' @param d Maximum number of digits of binary number
#'
#' @returns Matrix of binary numbers
#' @export

binary_matrix <- function(d) {
  nums <- 0:(2^d - 1)
  mat <- matrix(as.integer(intToBits(nums)), nrow = length(nums), byrow = TRUE)
  mat[, d:1]
}

#' generate_binary_files
#'
#' Function creates three data files "binM25, diffM25, ndiff25". Needed to gene-
#' rate the transition probability matrix and state space using ranked_coal.
#' Preferably, read the binM25.rds file and apply the diffM25_from_binM25 and
#' ndiff25_from_binM25 functions to it.
#'
#' @export

generate_binary_files <- function()
{
  if(file.exists("binaryMatrixPreloads.RData")){
    load("binaryMatrixPreloads.RData")
  } else {
    binM25 <- binary_matrix(25)
    ndiff25 <- rowSums(binM25)
    # diffM25 <- t(apply(apply(binM25,1,cumsum),2,rev))
    diffM25 <- t(apply(apply(apply(binM25,1,rev),2,cumsum),2,rev))
    save(binM25,diffM25,ndiff25,file="binaryMatrixPreloads.RData")
    gc()
  }
}

#' diffM25_from_binM25
#'
#' function computes the diffM25 matrix from the binM25 matrix.
#'
#' @param binM25 a matrix of all binary numbers that can be represented
#' with 25 digits. Computed using binary_matrix(25). Alternatively can be loaded.
#'
#' @returns the diffM25 matrix.
#' @export

diffM25_from_binM25 <- function(binM25)
{
  t(apply(apply(apply(binM25,1,rev),2,cumsum),2,rev))
}

#' ndiff25_from_binM25
#'
#' Function computes the ndiffM25 vector from the binM25 matrix.
#'
#' @param binM25 a matrix of all binary numbers that can be represented
#' with 25 digits. Computed using binary_matrix(25). Alternatively can be loaded.
#'
#' @returns the ndiff25 vector.
#' @export

ndiff25_from_binM25 <- function(binM25)
{
  ndiff25 <- rowSums(binM25)
}

# R/load_data.R
#' Load precomputed binM25
#'
#' @return Numeric matrix
#' @export
load_binM25 <- function() {
  readRDS(system.file("extdata", "binM25.rds", package = "RankedCoalescent"))
}

#' Load precomputed diffM25
#'
#' @return Numeric matrix
#' @export
load_diffM25 <- function() {
  readRDS(system.file("extdata", "diffM25.rds", package = "RankedCoalescent"))
}

#' Load precomputed ndiff25
#'
#' @return Numeric matrix
#' @export
load_ndiff25 <- function() {
  readRDS(system.file("extdata", "ndiff25.rds", package = "RankedCoalescent"))
}

