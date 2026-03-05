#' reward_matrix_F
#'
#' Computes the reward matrix of the non fixed F-matrix entries from
#' equation (23) in NAME OF PAPER
#'
#' @param StSpM matrix of all states
#' @param Tiers A list of tiers made using make_tiers.
#' @param IndexList A list of all non fixed indices. Made using non_fixed_indices
#' @returns A nrow(StSpM) x (ncol(StSpM)-1)*(ncol(StSpM)-2)/2 reward matrix
#' @export

reward_matrix_F <- function(StSpM, Tiers, IndexList)
{
  Kn <- length(IndexList)
  RewardMatrixF <- matrix(0, nrow = nrow(StSpM), ncol = Kn)
  n <- ncol(StSpM)+1
  for(k in 1:Kn)
  {
    i <- IndexList[[k]][1]
    j <- IndexList[[k]][2]
    for(stateidx in Tiers[[j]])
    {
      RewardMatrixF[stateidx,k] <- StSpM[stateidx,(n-i)]
    }
  }
  RewardMatrixF
}

#' inititial_distribution
#'
#' Computes the initial distribution in the DPH framework
#'
#' @param StSpM matrix of all state
#' @returns an nrow(StSpM) probability with its first entry equal to 1
#' @export

initial_distribution <- function(StSpM)
{
  c(1,rep(0,nrow(StSpM)-1))
}

#' left_products_seq
#'
#' computes the left factors init_dist * U * diag(r_{i,j}) where (i,j) is a non
#' fixed index.
#'
#' @param InitDist initial distribution of the markov embedding
#' @param SeqProbs List of transition probability matrix between consecutive
#' tiers. The k'th entry is equal to T_{k-1,k}, k = 1,...,n-2. Can be computed
#' using seq_ranked_coalescent.
#' @param RewardMatrixF reward matrix of the non-fixed entries. Can be computed
#' using reward_matrix_F.
#' @param Tiers A list of tiers made using make_tiers.
#' @param IndexList A list of all non fixed indices. Made using non_fixed_indices
#' @returns A list of all the products pi*U*diag(r_{i,j}) where i,j is a non
#' fixed index.
#' @export

left_products_seq <- function(InitDist, SeqProbs, RewardMatrixF, Tiers, IndexList)
{
  n <- length(Tiers) + 1
  PiGreen <- matrix(0, nrow = n-1, ncol = length(InitDist))
  PiGreen[1,] <- InitDist
  for(k in 1:(n-2))
  {
    PiGreen[(k+1), Tiers[[(n-(k+1))]]] <- PiGreen[k,Tiers[[(n-k)]]] %*% SeqProbs[[k]]
  }
  LeftVector <- colSums(PiGreen)
  LeftProdList <- list()
  Kn <- ncol(RewardMatrixF)
  for(j in 1:Kn)
  {
    CurrentTier <- IndexList[[j]][2]
    LeftProd <- numeric(length(InitDist))
    LeftProd[Tiers[[CurrentTier]]] <- LeftVector[Tiers[[CurrentTier]]] * RewardMatrixF[Tiers[[CurrentTier]], j]
    LeftProdList[[j]] <- LeftProd
  }
  return(LeftProdList)
}

#' right_products_seq
#'
#' Computes the right factors U * diag(r_{i,j}) * e for all non fixed entries
#' i,j.
#'
#' @param SeqProbs List of transition probability matrix between consecutive
#' tiers. The k'th entry is equal to T_{k-1,k}, k = 1,...,n-2. Can be computed
#' using seq_ranked_coalescent.
#' @param RewardMatrixF reward matrix of the non-fixed entries. Can be computed
#' using reward_matrix_F.
#' @param Tiers A list of tiers made using make_tiers.
#' @param IndexList A list of all non fixed indices. Made using non_fixed_indices
#' @returns A list of all the products U*diag(r_{i,j})*e where i,j is a non
#' fixed index.
#' @export

right_products_seq <- function(SeqProbs, RewardMatrixF, Tiers, IndexList)
{
  n <- length(Tiers) + 1
  nstates <- nrow(RewardMatrixF)
  Kn <- ncol(RewardMatrixF)
  right_products <- list()
  for(j in 1:Kn)
  {
    currentTier <- IndexList[[j]][2]
    rightProdMatJ <- matrix(0, nrow = n-1, ncol = nstates)
    rightProdMatJ[1,] <- RewardMatrixF[,j]
    for(k in 1:(n-1-currentTier))
    {
      rightProdMatJ[(k+1),Tiers[[(k+currentTier)]]] <- SeqProbs[[(n-currentTier-k)]] %*%  rightProdMatJ[k,Tiers[[(k-1+currentTier)]]]
    }
    right_products[[j]] <- colSums(rightProdMatJ)
  }
  return(right_products)
}

#' mean_tree_DPH
#'
#' Computes the mean F-matrix with `M_{ij} = E[F_{ij}]` where the columns of
#' `F` are sampled from the ranked coalescent with transition probabilities
#' used in `LeftProducts`.
#'
#' @param n the number of leaves/lineages.
#' @param LeftProducts A list of the products init_dist * U * diag(r_{i,j})
#' where (i,j) is a non-fixed index. Computed via left_products_seq
#' @param IndexList A list of all non fixed indices. Made using `non_fixed_indices`
#' @returns The mean F matrix of the distribution with TPM used to compute
#' `LeftProducts`.
#' @export

mean_tree_DPH <- function(n, LeftProducts, IndexList)
{
  MeanTree <- matrix(0, n-1, n-1)
  diag(MeanTree) <- 2:n
  for(i in 1:(n-2))
  {
    MeanTree[(i+1), i] <- i
  }
  Kn <- length(IndexList)
  for(k in 1:Kn)
  {
    MeanTree[IndexList[[k]][1], IndexList[[k]][2]] <- sum(LeftProducts[[k]])
  }
  return(MeanTree)
}

#' mean_tree_DPH_nofixed
#'
#' Computes the non fixed entries of the mean F-matrix with `M_{ij} = E[F_{ij}]`
#' where the columns of `F` are sampled from the ranked coalescent with
#' transition probabilities used in `LeftProducts`.
#'
#' @param LeftProducts A list of the products init_dist * U * diag(r_{i,j})
#' where (i,j) is a non-fixed index. Computed via left_products_seq
#' @param IndexList A list of all non fixed indices. Made using non_fixed_indices
#' @returns The mean F matrix of the distribution with TPM used to compute
#' LeftProducts.
#' @export

mean_tree_DPH_nofixed <- function(LeftProducts, IndexList)
{
  Kn <- length(IndexList)
  MeanTreeNoFixed <- numeric(length(Kn))
  for(i in 1:Kn)
  {
    MeanTreeNoFixed[i] <- sum(LeftProducts[[i]])
  }
  return(MeanTreeNoFixed)
}

#' covariance_matrix_F
#'
#' Computes the covariance matrix of the non fixed F-matrix indices.
#'
#' @param LeftProducts A list of the products init_dist * U * diag(r_{i,j})
#' where (i,j) is a non-fixed index. Computed via left_products_seq
#' @param RightProducts A list of the products U * diag(r_{i,j}) * e
#' where (i,j) is a non-fixed index. Computed via right_products_seq
#' @param IndexList A list of all non fixed indices. Made using non_fixed_indices
#' @param RewardMatrixF reward matrix of the non-fixed entries. Can be computed
#' using reward_matrix_F.
#' @returns the covariance matrix of all non fixed entries
#' @export

covariance_matrix_F <- function(LeftProducts, RightProducts, IndexList, RewardMatrixF)
{
  Kn <- length(IndexList)
  CovMatNoFixed <- matrix(0, Kn, Kn)
  for(i in 1:Kn)
  {
    for(j in 1:i)
    {
      term1 <- crossprod(LeftProducts[[i]], RightProducts[[j]])
      term2 <- crossprod(LeftProducts[[j]], RightProducts[[i]])
      term3 <- -crossprod(LeftProducts[[i]],RewardMatrixF[,j])
      term4 <- -sum(LeftProducts[[i]])*sum(LeftProducts[[j]])
      CovMatNoFixed[i,j] <- term1 + term2 + term3 + term4
    }
  }
  for(j in 2:Kn)
  {
    for(i in 1:(j-1))
    {
      CovMatNoFixed[i,j] <- CovMatNoFixed[j,i]
    }
  }
  return(CovMatNoFixed)
}

#' mean_S
#'
#' computes the mean of the sum of non fixed entries
#'
#' @param MeanTreeNoFixed vector of means of all non fixed entries
#' @returns the mean of the sum of non fixed entries
#' @export

mean_S <- function(MeanTreeNoFixed){sum(MeanTreeNoFixed)}

#' var_S
#'
#' computes the variance of the sum of non fixed entries
#'
#' @param CovMatNoFixed the covariance matrix of all non fixed entries
#' @returns the variance of the sum of the non fixed entries
#' @export

var_S <- function(CovMatNoFixed){sum(CovMatNoFixed)}

#' mean_E
#'
#' computes the mean external branch length
#'
#' @param MeanTree mean F-matrix
#' @returns the mean external branch length
#' @export

mean_E <- function(MeanTree){sum(MeanTree[nrow(MeanTree),])}

#' var_E
#'
#' computes the variance of the external branch length
#'
#' @param CovMatNoFixed the covariance matrix of all non fixed entries
#' @param IndexList A list of all non fixed indices. Made using non_fixed_indices
#' @returns the variance of the external branch length
#' @export

var_E <- function(CovMatNoFixed, IndexList)
{
  Kn <- length(IndexList)
  n <- (5 + sqrt(1 + 8*Kn))/2
  res <- 0
  for(i in 1:Kn)
  {
    for(j in 1:Kn)
    {
      if(IndexList[[i]][1] == (n-1) & IndexList[[j]][1] == (n-1))
      {
        res <- res + CovMatNoFixed[i,j]
      }
    }
  }
  res
}

#' covariance_matrix_SE
#'
#' Computes the 2 x 2 covariance matrix of the external branch length and sum
#' of non-fixed entries.
#'
#' @param CovMatNoFixed the covariance matrix of all non fixed entries
#' @param IndexList A list of all non fixed indices. Made using non_fixed_indices
#' @returns the covaraince matrix of the external branch length and sum of
#' non-fixed entries
#' @export

covariance_matrix_SE <- function(CovMatNoFixed, IndexList)
{
  Kn <- length(IndexList)
  n <- (5 + sqrt(1 + 8*Kn))/2
  res <- 0
  for(i in 1:Kn)
  {
    for(j in 1:Kn)
    {
      if(IndexList[[j]][1] == (n-1))
      {
        res <- res + CovMatNoFixed[i,j]
      }
    }
  }
  VarS <- var_S(CovMatNoFixed)
  VarE <- var_E(CovMatNoFixed, IndexList)
  return(matrix(c(VarS, res, res, VarE), 2, 2))
}

#' min_S
#'
#' Computes the minimum sum of non-fixed entries
#'
#' @param n the number of leaves
#' @returns minimum sum of non fixed entries
#' @export

min_S <- function(n)
{
  res <- 0
  for(j in 1:(n-1))
  {
    if(j < floor(n/2))
    {
      res <- res + (j+1)*(j+2)/2 - (2*j+1)
    }
    else
    {
      res <- res + (j+1)*(j+2)/2 - (2*j+1) - (2*j+1-n)*(2*j+2-n)/2
    }
  }
  res
}

#' max_S
#'
#' Computes the minimum sum of non-fixed entries
#'
#' @param n the number of leaves
#' @returns maximum sum of non fixed entries
#' @export

max_S <- function(n){
  res <- 0
  for(j in 1:(n-3))
  {
    res <- res + j * (n-2-j)
  }
  res
}

#' min_E
#'
#' Computes the minimum sum of non-fixed entries
#'
#' @param n the number of leaves
#' @returns minimum external branch length
#' @export

min_E <- function(n){
  1/4 * n^2 + 1/2 * n + 1/4 * (n %% 2)
}

#' max_E
#'
#' Computes the minimum sum of non-fixed entries
#'
#' @param n the number of leaves
#' @returns maximum external branch length
#' @export

max_E <- function(n){
  1/2 * n^2 - 1/2 * n + 1
}

#' reward_E_BCP
#'
#' Computes the reward transform for E based on the BCP
#'
#' @param StSpMBCP State space matrix of the ranked BCP
#' @returns The reward vector of the external branch length
#' @export

reward_E_BCP <- function(StSpMBCP)
{
  StSpMBCP[,1]
}

#' DPH_E
#'
#' Computes the sub-TPM and initial distribution of E
#'
#' @param BCP An object of the discrete_phase_type class. Computed using
#' the DPH function from PhaseTypeR
#' @param RewardE The reward transform of the external branch length. Computed
#' using reward_E_BCP
#' @returns The sub-TPM and initial distribution of E
#' @export

DPH_E <- function(BCP, RewardE)
{
  Y <- PhaseTypeR::reward_phase_type(BCP, RewardE)
  return(Y)
}

#' dist_E
#'
#' Computes the probabilities P(E = k) for all k between min_E and max_E
#'
#' @param TPME Sub-TPM of E. Can be computed using PhaseTypeR
#' @param InitDist initial distribution of E. Can be computed using PhaseTypeR
#' @param MinE the minimum external branch length
#' @param MaxE the maximum external branch length
#' @returns A vector of length max_E - min_E, with entry i = P(E = min_E - 1 + i)
#' @export

dist_E <- function(TPME, InitDist, MinE, MaxE)
{
  t <- rep(1, length(InitDist)) - rowSums(TPME)
  res <- InitDist %*% TPME
  ProbMat <- res
  for(i in 1:(MaxE-1))
  {
    res <- res %*% TPME
    ProbMat <- rbind(ProbMat,res)
  }
  Probs <- ProbMat[(MinE-1):(MaxE-1),] %*% t
}

#' binned_probs_E
#'
#' Bins the probabilities computed above. Used in the GoF test
#'
#' @param m number of observations in sample
#' @param Probs distribution of E
#' @returns A list of the bins as well as the expected number in each bin
#' @export

binned_probs_E <- function(m, Probs)
{
  Expected <- m * Probs
  counter <- 1
  BinList <- list()
  BinCounter <- 0
  while(T)
  {
    BinCounter <- BinCounter + 1
    currentSum <- 0
    currentBin <- c()
    if(counter > length(Expected))
    {
      break
    }
    for(i in counter:length(Expected))
    {
      currentSum <- currentSum + Expected[counter]
      currentBin <- c(currentBin, counter)
      counter <- counter + 1 #Number of variables used
      if(currentSum > 5)
      {
        BinList[[BinCounter]] <- currentBin
        break
      }
      if(counter == (length(Expected)+1))
      {
        BinList[[BinCounter]] <- currentBin
      }
    }
  }
  ExpectedBinned <- unlist(lapply(BinList, function(bin){sum(Expected[bin])}))
  if(ExpectedBinned[length(ExpectedBinned)] < 1)
  {
    BinList <- BinList[-length(BinList)]
    minItem <- min(BinList[[length(BinList)]])
    BinList[[length(BinList)]] <- minItem:length(Expected)
  }
  return(list("Bins" = BinList, "ExpectedBinned" = ExpectedBinned))
}


