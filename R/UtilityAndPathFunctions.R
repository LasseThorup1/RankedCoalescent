#' check_path
#'
#' Checks if a path has probability > 0
#'
#' @param path a vector of length n-1
#' @param TPM a transition probability matrix
#' @returns TRUE/FALSE depending on whether the path is possible
#' @export

check_path <- function(path, TPM)
{
  for(t in 2:length(path))
  {
    if(TPM[path[(t-1)],path[t]] == 0)
    {
      return(FALSE)
    }
  }
  return(TRUE)
}

#' make_tiers
#'
#' Puts states into tiers
#'
#' @param StSpM matrix of all states
#' @returns A list of all n-1 tiers
#' @export

make_Tiers <- function(StSpM)
{
  n <- max(StSpM)
  Tiers <- list()
  for(k in 2:n)
  {
    tier_idx <- c()
    for(i in 1:nrow(StSpM))
    {
      if(max(StSpM[i,]) == k)
      {
        tier_idx <- c(tier_idx, i)
      }
    }
    Tiers[[(k-1)]] <- tier_idx
  }

  return(Tiers)
}

#' make_cost_tiers
#'
#' Computes c(s) for all states s
#'
#' @param Tiers A list of tiers made using make_tiers
#' @param StSpM matrix of all states
#' @param MeanTree mean F-matrix from a given distribution
#' @returns a list of all costs in the same format is tiers
#' @export

make_cost_tiers <- function(Tiers, StSpM, MeanTree)
{
  CostTiers <- list()
  n <- max(StSpM)
  for(k in 2:n)
  {
    costs <- c()
    for(row in Tiers[[(k-1)]])
    {
      costs <- c(costs, sum((rev(StSpM[row,])-MeanTree[,(k-1)])^2))
    }

    CostTiers[[(k-1)]] <- costs
  }

  return(CostTiers)
}

#' path_to_Fmat
#'
#' Translates a path to its corresponding F-matrix
#'
#' @param path a vector of length n-1
#' @param StSpM matrix of all states
#' @returns the F-matrix corresponding to the path
#' @export

path_to_Fmat <- function(path, StSpM)
{
  n <- ncol(StSpM)
  Fmat <- matrix(numeric(n^2), n, n)
  for(i in 1:length(path))
  {
    Fmat[,i] <- rev(StSpM[path[i],])
  }

  return(Fmat)
}

#' Fmat_to_path
#'
#' Translates an F-matrix to its corresponding ordered path
#'
#' @param Fmat an (n-1) x (n-1) F-matrix
#' @param StSpM matrix of all states
#' @param Tiers A list of tiers made using make_tiers. Added for efficiency
#' @returns ordered path
#' @export

Fmat_to_path <- function(Fmat, StSpM, Tiers)
{
  path <- c(1,2)
  for(t in 3:(ncol(Fmat)))
  {
    path[t] <- which(StSpM[,Tiers[[t]]] == rev(Fmat[,t]))
  }
  path
}

#' non_fixed_indices
#'
#' computes all the non-fixed entries and returns it as a list
#'
#' @param n the number of leaves / lineages in the tree
#' @returns A list of all the non fixed indices
#' @export

non_fixed_indices <- function(n)
{
  IndexList <- list()
  counter <- 0
  for(i in 3:(n-1))
  {
    for(j in 1:(i-2))
    {
      counter <- counter + 1
      IndexList[[counter]] <- c(i,j)
    }
  }
  IndexList
}

#' F_no_fixed
#'
#' computes a vector of all non fixed entries in a given F-matrix
#'
#' @param Fmat An (n-1) x (n-1) F-matrix
#' @param IndexList A list of all non fixed indices.
#' Made using `non_fixed_indices`
#' @returns A vector of the non fixed entries
#' @export

F_no_fixed <- function(Fmat, IndexList)
{
  noFixed <- numeric(length(IndexList))
  for(i in 1:length(IndexList))
  {
    idx <- IndexList[[i]]
    noFixed[i] <- Fmat[idx[1], idx[2]]
  }
  noFixed
}

#' external_branch_length
#'
#' computes the external branch length of tree based on its F-matrix
#'
#' @param Fmat An (n-1) x (n-1) F-matrix
#' @returns the external branch length of the ranked and unlabeled tree.
#' @export

external_branch_length <- function(Fmat)
{
  sum(Fmat[nrow(Fmat),])
}

#' sum_nonfixed_entries
#'
#' computes the sum of non fixed entries in an F-matrix
#'
#' @param Fmat An (n-1) x (n-1) F-matrix
#' @returns the sum of non fixed entries
#' @export

sum_nonfixed_entries <- function(Fmat)
{
  n <- nrow(Fmat) + 1
  sum(Fmat) + 1 - n*(n+1)/2 - (n-1)*(n-2)/2
}

#' DfromF
#'
#' Transforms an F-Matrix into a D-matrix
#'
#' @param Fmat Valid F-Matrix
#'
#' @returns D-Matrix corresponding to the same ranked unlabelled tree
#'
#' @export

DfromF <- function(Fmat){
  if(!is.matrix(Fmat))stop("Please provide a matrix")
  n <- nrow(Fmat)+1
  Dmat <- diag(rep(2,n-1))
  Fmat <- cbind(0,Fmat)
  for(i in 1:(n-1)){for(j in 1:(i-1)){
    Dmat[i,j] <- Fmat[i,j+1]-Fmat[i,j]
  }}
  return(Dmat)
}

#' FfromD
#'
#' Transforms a D-Matrix into an F-matrix
#'
#' @param Dmat Valid D-Matrix
#'
#' @returns F-Matrix corresponding to the same ranked unlabelled tree
#' @export

FfromD <- function(Dmat){
  if(!is.matrix(Dmat))stop("Please provide a valid matrix")
  n <- nrow(Dmat)+1
  Fmat <- diag(seq(2,n))
  for(i in 1:(n-1)){for(j in 1:(i-1)){
    Fmat[i,j] <- sum(Dmat[i,1:j])
  }}
}



#' FtoBlock
#'
#' Transforms an F-matrix into a representation of the block-counting process (giving SFS of each coalescence event)
#'
#' @param Fmat Valid F-matrix
#'
#' @returns Block Counting Matrix
#' @export


FtoBlock <- function(Fmat){
  if(!is.matrix(Fmat))stop("Please provide a matrix")
  coals <- nrow(Fmat)
  # backup <- copy(Fmat)
  # Bmat <- matrix(0,coals,coals)
  # Bmat[coals,] <- Fmat[coals,]

  Dmat <- DfromF(Fmat)
  upper <- matrix(0,coals-1,coals)
  upper[coals-1,coals] <- 1
  splitType <- c(rep(0,coals-1),2)
  # curCol <- coals-1
  # while()
  for(curCol in (coals-1):2){
    upper[,curCol] <- upper[,curCol+1]
    Dcol <- Dmat[,curCol]
    if(Dmat[coals,curCol]==2){
      upper[coals-1,curCol] <- upper[coals-1,curCol]+1
      splitType[curCol] <- 2
    } else if(Dmat[coals,curCol]==1){
      splitId <- which(Dcol==1)[1]
      splitType[curCol] <- 1+splitType[splitId]
      upper[coals-splitType[splitId],curCol] <- upper[coals-splitType[splitId],curCol]+1
      upper[coals-splitType[splitId]+1,curCol] <- upper[coals-splitType[splitId]+1,curCol]-1
    } else {
      splitId <- which(Dcol==1)[1]
      # split0 <- Dcol[Dcol==0]

      id0 <- which(Dcol==0)
      split0 <- id0[id0>splitId][1]

      splitType[curCol] <- splitType[split0]+splitType[splitId]

      upper[coals-splitType[curCol]+1,curCol] <- upper[coals-splitType[curCol]+1,curCol]+1
      upper[coals-splitType[splitId]+1,curCol] <- upper[coals-splitType[splitId]+1,curCol]-1
      upper[coals-splitType[split0]+1,curCol] <- upper[coals-splitType[split0]+1,curCol]-1
    }
  }
  return(rbind(cbind(upper[,-1],0),Fmat[coals,]))
}

#' toBool
#'
#' Ensures that some input value is boolean
#'
#' @param val some 1 dimensional object
#'
#' @returns TRUE or FALSE, FALSE only if val is some form of FALSE/0
#' @export

toBool <- function(val){
  if(!(is.null(dim(val))&length(val)==1))stop("Provide boolean values for the boolean parameters")
  return(as.integer(val)!=0)
}

#' checkN
#'
#' Checks input for being an integer between upper and lower bounds
#'
#' @param n the input
#' @param lower the lower bound of the input
#' @param upper the upper bound of the input
#'
#' @returns integer n or raises Error
#' @export

checkN <- function(n, lower = 3, upper = 27) {
  if (!(is.numeric(n) &&
        length(n) == 1 &&
        is.null(dim(n)) &&
        n == as.integer(n) &&
        n >= lower &&
        n <= upper)) {
    stop(paste0("Please provide an integer n between ", lower, " and ", upper))
  }
  else{n}
}


#' getAdditiveDiff
#'
#' For a certain tier, gives the required differenceVectors
#'
#' @param sdiag the index of the subdiagonal
#' @param sdVal the value of the subdiagonal
#'
#' @returns Additive vectors that can be used to create new states
#' @export

getAdditiveDiff <- function(sdiag,sdVal){
  d <- sdiag-1
  maxDiff <- max(min(d,sdVal),0)
  indices <- c()
  for(i in 1:(2^(d))){
    if(ndiff25[i]<=maxDiff) indices <- c(indices,i)
  }
  return(diffM25[indices,1:d])
}

#' fibP1
#'
#' Gives the (n+1)'th Fibonacci number minus 1
#'
#' @param n The index of the fibonacci number -1
#'
#' @returns Fib(n+1)-1
#' @export

fibP1 <- function(n){#
  a <- 0
  b <- 0
  c <- 0
  for(i in 1:(n-1)){
    c <- b
    b <- a+b+1
    a <- c
  }
  return(b)
}

#' mergeRep
#'
#' Difference operator returning a low dimensional summary of two vectors coalescence potential
#'
#' @param a column in which the process is
#' @param b column into which the process might transition
#'
#' @returns Vector of differences, all but the last element reversed
#' @export

mergRep <- function(a,b){
  subVec <- a-b
  c(subVec[length(subVec)],diff(rev((subVec)[-1])))
}

#' meanF_from_sample
#'
#' computes the mean F-matrix of a sample
#'
#' @param Flist A list of F-matrices
#'
#' @returns The mean F matrix of the sample
#' @export

meanF_from_sample <- function(Flist)
{
  Reduce("+", Flist) / length(Flist)
}

