#' ViTreebi
#'
#' Computes the minimum attainable distance to the min tree and the paths
#' attaining this minimum.
#'
#' @param n the number of leaves / lineages in the trees
#' @param AntecedentList a list of length Fib(n+1)-1. Each entry contains the
#' possible previous states. Can be computed using BigRankedCoal
#' @param StSpM matrix of all states
#' @param CostTiers A list of costs for each state. Made using make_cost_tiers
#' @param Tiers A list of tiers made using make_tiers.
#' @returns all Frechet mean paths as well as the minimal attainable cost
#' @export

ViTreebi <- function(n, AntecedentList, StSpM, CostTiers, Tiers)
{

  m <- nrow(StSpM)

  CostMat <- matrix(rep(Inf,(n-1)*m), nrow = m, ncol = n-1)
  CostMat[1,1] <- 0
  CostMat[2,2] <- 0

  for(k in 3:(n-1))
  {
    for(j in 1:length(Tiers[[(n-k)]]))
    {
      current_tier_start <- max(Tiers[[(n-k+2)]])+1
      current_tier_end <- min(Tiers[[(n-k)]])-1

      state <- Tiers[[(n-k)]][j]

      possible_states <- AntecedentList[[state]] #This should be changed.
      CostMat[state, k] <- CostTiers[[(n-k)]][j] + min(CostMat[possible_states,(k-1)])
    }
  }

  PathList <- list()

  path_starts <- which(CostMat[,(n-1)] == min(CostMat[,(n-1)]))
  counter <- 0
  for(i in 1:length(path_starts))
  {
    current <- path_starts[i]
    possible_PathList <- list(c(current))

    for(k in 2:(n-1))
    {
      for(pos in current)
      {
        Antecedents <- AntecedentList[[pos]]
        back_arrows <- Antecedents[which(CostMat[Antecedents,(n-k)] == min(CostMat[Antecedents,(n-k)]))]
      }
      possible_PathList[[k]] <- back_arrows
      current <- back_arrows
    }

    M <- t(cbind(expand.grid(possible_PathList)))
    minimal_paths <- unname(split(M, rep(1:ncol(M), each = nrow(M))))

    for(path in minimal_paths)
    {
      counter <- counter + 1
      PathList[[counter]] <- path
    }
  }
  return(list("FrechetPaths" = PathList, "MinCost" = min(CostMat[,(n-1)])))
}
