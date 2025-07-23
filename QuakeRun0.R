  require(MASS)
  CP <- res$CP
  Q <- res$Q
  theta <- res$theta
  changepoints <- Changepoints <- intensity <- jumprates <- CPnum <-  threshold <- NULL
  ################################################################################
  ##############################################################################
  ##############################################################################
  for (j in 1:iter){

  cat("\n iteration", j, "\n")

################################################################################
##############################################
  if (res$Number != 1){  
  rates <- -diag(Q)
  rates <- rates[-length(rates)]
  rates <- c(rates, 1/(ending-starting))
  flag <- FALSE
  while (flag==FALSE){
  A <- Uniform(CP, rates)
  Uniformtimes <- A$U
  index <- which(diff(c(starting,Uniformtimes,ending)) <= 0.001)
  Uniformtimes <- Uniformtimes[-index]
  JumpProb <- A$JumpProb[-index, ]
  if ((max(diff(Uniformtimes)) < 5) && (length(Uniformtimes) <= 300) && (length(Uniformtimes) >=3)) flag <- TRUE
  }
  }
  
  if (res$Number == 1){
  flag <- FALSE
  while (flag==FALSE){
  A <- Uniform(c(starting, ending), 1/(ending-starting))
  Uniformtimes <- A$U
  index <- which(diff(c(starting,Uniformtimes,ending)) <= 0.001)
  Uniformtimes <- Uniformtimes[-index]
  JumpProb <- A$JumpProb[-index, ]
  if ((max(diff(Uniformtimes)) < 5) && (length(Uniformtimes) <= 300) && (length(Uniformtimes) >=3) ) flag <- TRUE
  }
  }
################################################################################ 

  z <- res <- FFBS(Uniformtimes, JumpProb, eventtimes, obs)
  num <- res$Number
  cp <- res$cp
  CP <- res$CP
  Q <- res$Q 
  KK <- length(CP)
  if ((ending == CP[KK-1])&&(res$Number != 1)){
  CP <- CP[-KK]
  cp <- cp[-length(cp)]
  Q <- Q[,-(KK-1)]
  Q <- Q[-(KK-2),]
  num <- num - 1
  }  
  Lambda <- res$Lambda
  if (j > burnin){
  changepoints <- c(changepoints, list(cp))
  Changepoints <- c(Changepoints, list(CP))
  jumprates <- c(jumprates, list(Q))
  CPnum <- c(CPnum, num)
  intensity <- c(intensity, list(Lambda))
  threshold <- c(threshold, list(cuts))
        }
  
################################################################################
  } 