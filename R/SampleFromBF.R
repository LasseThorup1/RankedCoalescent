#' gen_Fmats_from_BF
#'
#' Samples F-matrices from the Blum-Francois model. The code is taken from
#' the f-matrix package, but changed so that it outputs in the desired format.
#'
#' @param b the parameter beta in the Blum-Francois model
#' @param m the size of sample. 1000 by default
#' @param n the size of the trees (number of leaves). 25 by default
#' @returns a sample of m F-matrices sampled from the Blum-Francois model
#' @export

gen_Fmats_from_BF <- function(b, m = 1000, n = 25){
  x <- fmatrix::rEncod(m = m, n = n, b = b)
  Fmats <- lapply(x, fmatrix::Fmat_from_myencod)
}
