#' wald_test_F
#'
#' Performs the wald test based on all non-fixed entries as described in
#' Theorem 6.1 in NAME OF PAPER.
#'
#' @param Fmats a list of (n-1) x (n-1) F-matrices
#' @param M0 the mean F-matrix under the null hypothesis. Can be
#' computed using mean_tree_DPH_nofixed
#' @param Sigma0invSqrt The inverse square root of the covariance matrix of
#' all non fixed entries under the null hypothesis. Can be computed using
#' covariance_matrix_F
#' @param IndexList A list of all non fixed indices.
#' Made using `non_fixed_indices`
#' @returns the p-value of the test
#' @export

wald_test_F <- function(Fmats, M0, Sigma0invSqrt,IndexList){
  # FmatsNonFixed <- lapply(Fmats, function(Fmat){F_no_fixed(Fmat, IndexList)})
  # FmatsMat <- rlist::list.rbind(FmatsNonFixed)
  m <- length(Fmats)
  Kn <- length(IndexList)
  # meanObs <- colMeans(FmatsMat)
  meanObs <- F_no_fixed(meanF_from_sample(Fmats), IndexList)
  W <- sqrt(m/Kn)*sum(Sigma0invSqrt%*%(meanObs-M0))
  pval <- 2*(1-pnorm(abs(W)))
  return(pval)
}

#' wald_test_SE
#'
#' Performs the wald test based on the balance indices S, E as described in
#' Theorem 6.1 in NAME OF PAPER.
#'
#' @param Fmats a list of (n-1) x (n-1) F-matrices
#' @param Mu0SE a vector of mean(S), mean(E) under the null hypothesis. Can be
#' computed using mean_E and mean_S
#' @param Sigma0SEinvSqrt The inverse square root of the covariance matrix
#' of E and S under the null hypothesis. Compute using covariance_matrix_SE
#' @param IndexList A list of all non fixed indices.
#' Made using `non_fixed_indices`
#' @returns the p-value of the test
#' @export

wald_test_SE <- function(Fmats, Mu0SE, Sigma0SEinvSqrt,IndexList)
{
  S <- unlist(lapply(Fmats, function(Fmat){sum(F_no_fixed(Fmat, IndexList))}))
  E <- unlist(lapply(Fmats, function(Fmat){sum(Fmat[nrow(Fmat),])}))
  m <- length(Fmats)
  W <- sqrt(m/2) * sum(Sigma0SEinvSqrt %*% (c(mean(S),mean(E))-Mu0SE))
  pval <- 2*(1-pnorm(abs(W)))
  return(pval)
}

#' GoF_E
#'
#' Performs the goodnees of fit test based on the balance index E as described
#' in Theorem 6.1 in NAME OF PAPER.
#'
#' @param Fmats a list of (n-1) x (n-1) F-matrices
#' @param Bins The bins needed to compute the GoF test statistics. Made using
#' make_bins function.
#' @param ExpectedBinned the expected number of observation within each bin
#' under the null hypothesis. Can be computed using the make_bins function
#' @param MinE the minimum external branch length
#' @param MaxE the maximum external branch length
#' @returns the p-value of the GoF test
#' @export

GoF_E <- function(Fmats, Bins, ExpectedBinned, MinE, MaxE){
  range <- MinE:MaxE
  ExternalBranchLengths <- unlist(lapply(Fmats, function(matrix){sum(matrix[nrow(matrix),])}))
  ObservedBinned <- unlist(lapply(Bins, function(bin){sum(sapply(ExternalBranchLengths, function(x){x %in% range[bin]}))}))
  nonZeroIdx <- which(ObservedBinned > 0)
  Gtestor <- 2*sum(ObservedBinned[nonZeroIdx]*log(ObservedBinned[nonZeroIdx]/ExpectedBinned[nonZeroIdx]))
  pval <- (1 - pchisq(Gtestor, df = length(Bins)-1))
  return(pval)
}
