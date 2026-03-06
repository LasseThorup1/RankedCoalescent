#' diff_matrix_big_enough
#' 
#' Tests whether the difference Matrix is loaded into the package environment and if it's big enough for its current task
#' 
#' @param nInput The number of leaves of the state space to create
#' 
#' @returns A boolean answering whether there is a big enough difference matrix stored
#' 

diff_matrix_big_enough <- function(nInput=25){
  if(exists("binM25")&exists("diffM25")&exists("ndiff25")){
    if(max(which(binM25[nrow(binM25),]==1))<(nInput-2)) return(FALSE)
    return(TRUE)
  }
  return(FALSE)
}



#' preload_binary_matrices
#' 
#' Loads the necessary binary numbers with up to nInput-2 digits
#' 
#' @param nInput The number of leaves of the state space to create
#' 


preload_binary_matrices <- function(nInput=25){
  if(file.exists("binM25.rds")){
    binM25 <- readRDS(system.file("extdata", "binM25.rds", package = "RankedCoalescent"))
    diffM25 <- readRDS(system.file("extdata", "diffM25.rds", package = "RankedCoalescent"))
    ndiff25 <- readRDS(system.file("extdata", "ndiff25.rds", package = "RankedCoalescent"))
    if(nInput>2){
      maxIdx <- which(binM25[,nInput-2]==1)[1]
      binM25 <- binM25[1:maxIdx,1:(nInput-2)]
      diffM25 <- diffM25[1:maxIdx,1:(nInput-2)]
      ndiff25 <- ndiff25[1:maxIdx]
      gc()
    }
  } else {
    binM25 <- binary_matrix(nInput-2)
    ndiff25 <- rowSums(binM25)
    diffM25 <- t(apply(apply(apply(binM25,1,rev),2,cumsum),2,rev))
    gc()
  }
  return(list(binM25,diffM25,ndiff25))
}




#' ranked_coal
#'
#' Generates state space and transition probability matrix for given n
#'
#' @param n Number from 3 to 23. The number of lineages at the start of the
#' process
#' @param prob If false, returns a matrix where possible transitions are 1 and
#' impossible ones are 0 instead of a transition probability matrix
#' (this is used for the ViTreebi)
#'
#' @returns A list where the first element is the state space matrix, the
#' second is the transition probability/possibility matrix
#' @export

ranked_coal <- function(n,prob=TRUE){
  #-----------------------------------------------------------------------
  #-----------             Input Checks     -------------------------------------
  #-----------------------------------------------------------------------
  prob <- toBool(prob)
  n <- checkN(n)
  if(!diff_matrix_big_enough(nInput)){
    a <- preload_binary_matrices(nInput)
    binM25 <- a[[1]]
    diffM25 <- a[[2]]
    ndiff25 <- a[[3]]
    rm(a);gc()
  }
  #-----------------------------------------------------------------------
  #-----------             Setup     -------------------------------------
  #-----------------------------------------------------------------------
  ssm <- matrix(0,fibP1(n),n-1)
  probm <- matrix(0,nrow(ssm),nrow(ssm))
  probm[1,2] <- 1
  ssm[1,1] <- n;ssm[2,1] <- n-2;ssm[2,2] <- n-1
  diag <- 3
  nextrow <- 3
  prevrow <- 2
  binaryM <- matrix(1,1)
  if(n==3)return(list(StSpM=ssm,ProbM=probm))
  for(subVal in (n-3):1){

    d <- diag-2
    maxDiff <- max(min(d,subVal),0)
    indices <- c()
    for(i in 1:(2^(d))){
      if(ndiff25[i]<=maxDiff) indices <- c(indices,i)
    }
    prevBin <- binaryM
    binaryM <- binM25[indices,1:d]
    diffM <- diffM25[indices,1:d]
    if(is.null(nrow(binaryM)))binaryM <- as.matrix(binaryM)
    newMat <- cbind(subVal-diffM,subVal,subVal+1,deparse.level=0)
    ssm[nextrow:(nextrow+nrow(newMat)-1),1:diag] <-   cbind(subVal-diffM,subVal,subVal+1,deparse.level=0)
    prevStates <- 1:nrow(prevBin)
    nextStates <- 1:length(indices)

    for(state in 1:nrow(prevBin)){
      possibleNextStates <- nextStates
      for(cell in 1:ncol(prevBin)){
        if(!prevBin[state,cell]){
          possibleNextStates <- possibleNextStates[binaryM[possibleNextStates,cell]==0]
        }
      }
      curst <- ssm[state+prevrow-1,1:(diag-1)]
      possibleNextStates <- possibleNextStates[apply(ssm[possibleNextStates+nextrow-1,1:(diag-1)],1,function(i){
        all((curst-i)>=0)
      })]

      if(prob){
        prevVals <- c(ssm[state+prevrow-1,1],prevBin[state,])
        if(subVal<n-3)prevVals <- c(prevVals,1)
        denom <- sum(prevVals)*(sum(prevVals)-1)
        nextVals <- cbind(ssm[possibleNextStates+nextrow-1,1],binaryM[possibleNextStates,])
        singletonOrderMerge <- t(apply(nextVals,1,function(i)prevVals-i))
        probs <- apply(singletonOrderMerge,1,function(nxtRw){
          if(nxtRw[1]>1){
            return(prevVals[1]*(prevVals[1]-1)/denom)
          }
          mergeIdx <- which(nxtRw>0)
          return(2*prevVals[mergeIdx[1]]*prevVals[mergeIdx[2]]/(denom))

        })
      } else probs <- 1
      probm[state+prevrow-1,possibleNextStates+nextrow-1] <- probs
    }


    prevrow <- nextrow
    nextrow <- nextrow+nrow(newMat)
    diag <- diag+1
  }
  return(list(StSpM=ssm,ProbM=probm))
}



#' ranked_coal_list
#'
#' Same process as RankedCoal, but works for larger n by returning linked
#' lists instead of a matrix
#'
#' @param n Number from 3 to 25. The number of lineages at the start of the process
#' @param ChildAndProb If true output contains two lists. childrenList is
#' a list of length = #states - n-2, where each entry is the possible next
#' states for that state. probList is a list containing the corresponding
#' transition probabilities.
#' @param TierInd if TRUE returns a tier list (import for ranked_coal_seq)
#' @returns A list where the first element is a state space list, followed by
#' lists of antecedents,children,probabilities (in the same order as children),
#' and the tier of each state
#' @export

ranked_coal_list <- function(n,ChildAndProb=FALSE,TierInd=FALSE){
  ChildAndProb <- toBool(ChildAndProb); TierInd <- toBool(TierInd)
  n <- checkN(n,upper=Inf)
  if(!diff_matrix_big_enough(nInput)){
    a <- preload_binary_matrices(nInput)
    binM25 <- a[[1]]
    diffM25 <- a[[2]]
    ndiff25 <- a[[3]]
    rm(a);gc()
  }
  #-----------------------------------------------------------------------
  #-----------             Setup     -------------------------------------
  #-----------------------------------------------------------------------
  ssm <- matrix(0,fibP1(n),n-1)
  AntecedentList <- lapply(1:fibP1(n),function(i)NULL)
  childrenList <- NULL
  probList <- NULL
  tierList <- NULL
  if(TierInd) tierList <- c(1)
  if(ChildAndProb){
    childrenList <- list(2)
    probList <- lapply(1:fibP1(n),function(i)NULL)
  }
  AntecedentList[[2]] <- 1
  ssm[1,1] <- n;ssm[2,1] <- n-2;ssm[2,2] <- n-1
  diag <- 3
  nextrow <- 3
  prevrow <- 2
  binaryM <- matrix(1,1)

  #-----------------------------------------------------------------------
  #-----------             Loop      -------------------------------------
  #-----------------------------------------------------------------------
  for(subVal in (n-3):1){
    d <- diag-2

    ##Figures out all the binary combinations that create a valid state
    maxDiff <- max(min(d,subVal),0) #you cannot substract more than there are (as q_ij >=0)
    indices <- c()
    for(i in 1:(2^(d))){
      if(ndiff25[i]<=maxDiff) indices <- c(indices,i)
    }
    prevBin <- binaryM
    binaryM <- binM25[indices,1:d] #binary substraction matrix
    diffM <- diffM25[indices,1:d] #cumulative substracted values
    if(is.null(nrow(binaryM)))binaryM <- as.matrix(binaryM)

    newMat <- cbind(subVal-diffM,subVal,subVal+1,deparse.level=0)
    #assign the new states to the matrix
    ssm[nextrow:(nextrow+nrow(newMat)-1),1:diag] <-   cbind(subVal-diffM,subVal,subVal+1,deparse.level=0)

    #-----------------------------------------------------------------------
    #-----------  Transition Matrix    -------------------------------------
    #-----------------------------------------------------------------------
    nextStates <- 1:length(indices)
    for(state in 1:nrow(prevBin)){
      possibleStates <- nextStates
      for(cell in 1:ncol(prevBin)){
        if(!prevBin[state,cell]){ # from the rules you can extract one single rule
          #that governs possibility of transition, which is diff(q_{i,j})==0 -> diff(q_{(i-1),j}==0)
          possibleStates <- possibleStates[binaryM[possibleStates,cell]==0]
        }
      }
      curst <- ssm[state+prevrow-1,1:(diag-1)]
      possibleStates <- possibleStates[apply(ssm[possibleStates+nextrow-1,1:(diag-1)],1,function(i){
        all((curst-i)>=0)
      })]
      if(ChildAndProb){
        childrenList[[state+prevrow-1]] <- possibleStates+nextrow-1
        prevVals <- c(ssm[state+prevrow-1,1],prevBin[state,])
        if(subVal<n-3)prevVals <- c(prevVals,1)
        denom <- sum(prevVals)*(sum(prevVals)-1)
        nextVals <- cbind(ssm[possibleStates+nextrow-1,1],binaryM[possibleStates,])
        singletonOrderMerge <- t(apply(nextVals,1,function(i)prevVals-i))
        probs <- apply(singletonOrderMerge,1,function(nxtRw){
          if(nxtRw[1]>1){
            return(prevVals[1]*(prevVals[1]-1)/denom)
          }
          mergeIdx <- which(nxtRw>0)
          return(2*prevVals[mergeIdx[1]]*prevVals[mergeIdx[2]]/(denom))

        })
        probList[[state+prevrow-1]] <- probs
      }
      if(TierInd) tierList <- c(tierList,length(possibleStates))
      for(posState in possibleStates) AntecedentList[[posState+nextrow-1]] <- c(AntecedentList[[posState+nextrow-1]],state+prevrow-1)

    }

    probList[[1]] <- 1



    prevrow <- nextrow
    nextrow <- nextrow+nrow(newMat)
    diag <- diag+1
  }
  return(list(StSpM=ssm,AntecedentList=AntecedentList,childrenList=childrenList,probList=probList,tierList=tierList))
}

#' ranked_coal_seq
#'
#' Generates state space matrix and transition probability matrix for given n
#' the transition probability matrix is generated as a list, where the k'th
#' entry contains the transition probability matrix from tier (k-1) to tier k.
#'
#' @param n The number of lineages at the start of the process
#'
#' @returns A list containing the state space matrix, the TierInd, and the
#' transition probability matrix in list format.
#' @export

ranked_coal_seq <- function(n){
  brc <- ranked_coal_list(n,T)
  tierList <- brc$tierList
  tierList <- rev(unname(table(apply(brc$StSpM,1,function(i)which(rev(i)>0)[1]))))
  TierInd <- c(1,cumsum(c(tierList))+1)
  probs <- brc$probList
  chList <- brc$childrenList
  seqProb <- lapply(1:(n-1),function(i)NULL)
  seqInit <- lapply(1:(n-1),function(i)NULL)

  seqProb[[1]] <- matrix(c(0,0,1,0),2,2)
  seqProb[[2]] <- matrix(c(0,0,1,0),2,2)

  dimPrev <- 1
  dimCur <- 1

  for(tier in 2:(length(tierList)-1)){
    qMat <- matrix(0,tierList[tier],tierList[tier+1])
    toIndices <- lapply(chList[TierInd[tier]:((TierInd[tier+1])-1)],function(i){
      i-TierInd[tier+1]+1
    })
    for(i in 1:length(toIndices)){
      qMat[i,toIndices[[i]]] <- probs[[TierInd[tier]+i-1]]
    }
    seqProb[[tier]] <- qMat
  }
  seqProb <- seqProb[-length(seqProb)]
  seqProb[[1]] <- 1
  return(list(StSpM=brc$StSpM,TierInd=TierInd,seqProbs=seqProb))
}

#' seq_probM_to_probM
#'
#' Functions returns a transition probability matrix given a list of transition
#' probability matrices between consecutive tiers.
#'
#' @param seqProb a list of transition probability matrices between consecutive
#' tiers.
#'
#' @returns An (Fib(n+1)-1) x (Fib(n+1) - 1) transition probability matrix
#' @export

seq_probM_to_probM <- function(seqProb){
  resultMatrix <- matrix(0,ncol(seqProb)+nrow(seqProb),ncol(seqProb)+nrow(seqProb))
  resultMatrix[1:nrow(seqProb),(nrow(seqProb)+1):ncol(resultMatrix)] <- seqProb
  return(resultMatrix)
}



#' ranked_coal_StSpM_only
#'
#' Generates state space matrix for given n
#'
#' @param n The number of lineages at the start of the process
#'
#' @returns State space matrix of ranked unlabeled trees for n lineages
#' @export

ranked_coal_StSpM_only <- function(n){
  n <- checkN(n)
  if(!diff_matrix_big_enough(nInput)){
    a <- preload_binary_matrices(nInput)
    binM25 <- a[[1]]
    diffM25 <- a[[2]]
    ndiff25 <- a[[3]]
    rm(a);gc()
  }
  #-----------------------------------------------------------------------
  #-----------             Setup     -------------------------------------
  #-----------------------------------------------------------------------
  ssm <- matrix(0,fibP1(n),n-1)
  probm <- matrix(0,nrow(ssm),nrow(ssm))
  ssm[1,1] <- n;ssm[2,1] <- n-2;ssm[2,2] <- n-1
  diag <- 3
  nextrow <- 3
  binaryM <- matrix(1,1)
  for(subVal in (n-3):1){

    d <- diag-2
    maxDiff <- max(min(d,subVal),0)
    indices <- c()
    for(i in 1:(2^(d))){
      if(ndiff25[i]<=maxDiff) indices <- c(indices,i)
    }
    binaryM <- binM25[indices,1:d]
    diffM <- diffM25[indices,1:d]
    if(is.null(nrow(binaryM)))binaryM <- as.matrix(binaryM) 
    newMat <- cbind(subVal-diffM,subVal,subVal+1,deparse.level=0)
    ssm[nextrow:(nextrow+nrow(newMat)-1),1:diag] <-   cbind(subVal-diffM,subVal,subVal+1,deparse.level=0)

    nextrow <- nextrow+nrow(newMat)
    diag <- diag+1
  }
  return(ssm)
}







#' bcp_coal
#'
#' Computes state space and transition probability matrix for the
#' Block Counting process
#'
#' @param n the number of lineages at the start of the process
#'
#' @returns a list, the first element is the state space matrix,
#' the second the transition probability matrix
#' @export


bcp_coal <- function(n){
  n <- checkN(n,upper=Inf)
  stsp <- as.data.frame(matrix(c(n,rep(0,n-1)),1))
  rownames(stsp) <- "50000"
  ProbM <- matrix(0,BCPssp[n],BCPssp[n])
  oldIdxFrom <- 1
  oldIdxTo <- 1
  stspPerTier <- 1
  while(sum(stsp[nrow(stsp),])>1){
    for(oldRow in oldIdxFrom:oldIdxTo){
      curRow <- stsp[oldRow,]
      for(index in 1:(n-1)){
        curVal <- curRow[index]
        if(curVal>0){
          if(curVal>1){
            newRow <- curRow
            newRow[index] <- newRow[index]-2
            newRow[index*2] <- newRow[index*2]+1
            checksum <- paste0(newRow,collapse="")
            locRw <- which(checksum==rownames(stsp)[oldIdxTo:nrow(stsp)])+oldIdxTo-1
            if(length(locRw)==0){
              stsp <- rbind(stsp,newRow)
              rownames(stsp)[nrow(stsp)] <- checksum
              ProbM[oldRow,nrow(stsp)] <- exp(lchoose(as.numeric(curRow[index]),2)-lchoose(sum(curRow),2))

            } else {
              ProbM[oldRow,locRw] <- ProbM[oldRow,locRw]+exp(lchoose(as.numeric(curRow[index]),2)-lchoose(sum(curRow),2))
            }

          }
          for(greaterIndex in (index+1):ncol(stsp)){
            if(curRow[greaterIndex]>=1){
              newRow <- curRow
              newRow[index] <- newRow[index]-1
              newRow[greaterIndex] <- newRow[greaterIndex]-1
              newRow[index+greaterIndex] <- newRow[index+greaterIndex]+1
              checksum <- paste0(newRow,collapse="")
              locRw <- which(checksum==rownames(stsp)[oldIdxTo:nrow(stsp)])+oldIdxTo-1
              if(length(locRw)==0) {
                stsp <- rbind(stsp,newRow)
                rownames(stsp)[nrow(stsp)] <- checksum
                ProbM[oldRow,nrow(stsp)] <- exp((sum(log(curRow[c(index,greaterIndex)])))-lchoose(sum(curRow),2))
              } else{
                ProbM[oldRow,locRw] <- ProbM[oldRow,locRw]+exp((sum(log(curRow[c(index,greaterIndex)])))-lchoose(sum(curRow),2))
              }
            }
          }
        }
      }
    }
    stsp <- dplyr::distinct(stsp)
    oldIdxFrom <- oldIdxTo+1
    oldIdxTo <- nrow(stsp)
    stspPerTier <- c(stspPerTier,oldIdxTo-oldIdxFrom+1)
  }
  rownames(stsp)=NULL
  ProbM <- ProbM[-nrow(ProbM),-ncol(ProbM)]
  stsp <- stsp[-nrow(stsp),]
  return(list(StSpM=stsp,ProbM=ProbM))
}
