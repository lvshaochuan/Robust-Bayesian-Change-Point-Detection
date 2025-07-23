  CP <- res$CP
  Q <- res$Q
  theta <- res$param
  Tau <- res$Tau
  Alpha <- res$Alpha
  Ui <- res$Ui
  changepoints <- Changepoints <- intensity <- jumprates <- CPnum <- Taus <- Alphas <- NULL
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

  for (j in 1:iter){
  if (res$Number != 1){  
  rates <- -diag(Q)
  rates <- rates[-length(rates)]
  rates <- c(rates, 1/(ending-starting))
  flag <- FALSE
  while (flag==FALSE){
  A <- Uniform(CP, rates)
  Uniformtimes <- A$U
  index <- which(diff(c(starting,Uniformtimes,ending)) <= diff(eventtimes[1:2]))
  Uniformtimes <- Uniformtimes[-index]
  JumpProb <- A$JumpProb[-index, ]
  if ((length(Uniformtimes) >=3) && (length(Uniformtimes) <= 300)) flag <- TRUE
  }
  }
  if (res$Number == 1){
  
  flag <- FALSE
  while (flag==FALSE){
  A <- Uniform(c(starting, ending), 1/(ending-starting))
  Uniformtimes <- A$U
  index <- which(diff(c(starting,Uniformtimes,ending)) <= diff(eventtimes[1:2]))
  Uniformtimes <- Uniformtimes[-index]
  JumpProb <- A$JumpProb[-index, ]
  if ((length(Uniformtimes) >=3) && (length(Uniformtimes) <= 300)) flag <- TRUE
  }
  }
  print(j)
 # FFBS.t <- function(Uniformtimes, JumpProb, eventtimes, obs, param, DF, Tau, Ui, Alpha)
  res <- FFBS.t(Uniformtimes, JumpProb, eventtimes, obs, param=theta, DF=1, Tau, Ui, Alpha)
  Ind <- which(diff(res$cp)<=3)
  if (length(Ind)==0){
  cp <- res$cp
  CP <- res$CP
  num <- res$Number    
  Q <- res$Q    
  theta <- res$param    
      } else{
  cp <- res$cp[-Ind]
  CP <- res$CP[-Ind]
  num <- res$Number - length(Ind)
  P <- res$Q 
  Q <- matrix(rep(0, (nrow(P) - length(Ind))*(nrow(P) - length(Ind))), nrow=nrow(P)-length(Ind))
  i <- 1
  for (j in 1:nrow(P)){
       if(any(j == Ind)) next else{
          k <- sum(j < Ind)
          Q[i,] <- P[j, (length(Ind) - k+1):(ncol(P) - k)]
          i <- i + 1
          }   
       }
  theta <- res$param[-Ind]
  }
  
  Tau <- res$Tau
  Alpha <- res$Alpha
  Ui <- res$Ui
  if (j > burnin){
  changepoints <- c(changepoints, list(cp))
  Changepoints <- c(Changepoints, list(CP))
  jumprates <- c(jumprates, list(Q))
  CPnum <- c(CPnum, num)
  intensity <- c(intensity, list(theta))
  Taus <- c(Taus, Tau)
  Alphas <- c(Alphas, Alpha)
        }
  }      
    