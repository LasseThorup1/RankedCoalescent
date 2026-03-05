#' rRankedIslandCoal
#'
#' Simulates from the Ranked Structured Coalescent
#'
#' @param m Number of trees to simulate
#' @param n Number of lineages at start
#' @param K number of different islands
#' @param migrationRate Rate of migration
#' @param coalescenceRate Rate of coalescence
#' @param startDist If there is a fixed starting distribution, it can be given
#' using a vector of length K that sums to n. Otherwise, a random starting
#' distribution will be drawn using the multinomial distribution
#'
#' @returns A list of n ranked unlabeled trees that were drawn according to the specifications of the parameters
#' @export


rRankedIslandCoal <- function(m,n,K,migrationRate=1,coalescenceRate=1,startDist=NA){
  m <- checkN(m,1,Inf)
  n <- checkN(n,2,50)
  K <- checkN(K,1,Inf)
  if(!(is.numeric(migrationRate)&migrationRate>0))stop("Provide a valid migrationRate in R^+")
  if(!(is.numeric(coalescenceRate)&coalescenceRate>0))stop("Provide a valid coalescenceRate in R^+")
  #-----------------------------------------------------------------------
  #-----------             Setup     -------------------------------------
  #-----------------------------------------------------------------------
  resultTrees <- list()
  drawStart <- any(is.na(startDist))
  for(iter in 1:m){
    ntonsPerIsland <- matrix(0,K,n-1)
    if(drawStart){
      ntonsPerIsland[,1] <-  rmultinom(1,n,rep(1/K,K))
    } else {
      ntonsPerIsland[,1] <- startDist
    }
    pplPerIsland <- rowSums(ntonsPerIsland)
    diag <- 2
    tree <- matrix(0,n-1,n-1)
    tree[n-1,n-1] <- n
    col <- n-2
    #############
    # Sooo, I take the migration risk to be simply multiplicative with
    # the surviving branches, so exp(lambda*diagVal).
    # The coalescence risk however is multiplicative with all the possible
    # combinations, so sum over all islands: n_{island} choose 2
    # ???
    # As such, I will first evaluate whether migration or coalescence happens next?
    while(col>0){
      #####
      # Calculate the rates of coalescing and migrating, then choose what happens next
      #####
      coalRisk <- pmax(0,choose(pplPerIsland,2))/100
      if(sum(coalRisk)==0) coalescenceNext <- FALSE
      else coalescenceNext <- rexp(1,coalescenceRate*sum(coalRisk))<rexp(1,migrationRate*(n-diag+2)/100)

      if(coalescenceNext){
        #####
        # Draw an island with probability proportional to number of pairs
        # Then, merge a random pair from that island
        # The exact calculation stem from this binary matrix view we also use for statespace creation
        #####
        coalIsland <- sample.int(K,1,prob=coalRisk)

        coalNton1 <- sample.int(n-1,1,prob=ntonsPerIsland[coalIsland,])
        ntonsPerIsland[coalIsland,coalNton1] <- ntonsPerIsland[coalIsland,coalNton1]-1
        coalNton2 <- sample.int(n-1,1,prob=ntonsPerIsland[coalIsland,])
        ntonsPerIsland[coalIsland,coalNton2] <- ntonsPerIsland[coalIsland,coalNton2]-1
        pplPerIsland[coalIsland] <- pplPerIsland[coalIsland]-1

        ntonsPerIsland[coalIsland,diag] <- 1
        tree[(n-diag):(n-1),col] <- rev(cumsum(colSums(ntonsPerIsland))[1:diag])
        diag <- diag+1
        col <- col-1
      } else {
        #####
        # Draw a random person to migrate, and send them on their way
        #####
        emigratioK <- sample.int(K, 1, prob = pplPerIsland)
        emigratingNton <- sample.int(n-1, 1, prob = as.numeric(ntonsPerIsland[emigratioK, ]))
        ntonsPerIsland[emigratioK,emigratingNton] <- ntonsPerIsland[emigratioK,emigratingNton]-1

        immigrationTarget <- sample(c(1:K)[-emigratioK],1)
        ntonsPerIsland[immigrationTarget,emigratingNton] <- ntonsPerIsland[immigrationTarget,emigratingNton]+1

        pplPerIsland[emigratioK] <- pplPerIsland[emigratioK]-1
        pplPerIsland[immigrationTarget] <- pplPerIsland[immigrationTarget]+1
      }
    }
    resultTrees <- rlist::list.append(resultTrees,tree)
  }
  return(resultTrees)
}
