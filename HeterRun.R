   CP <- res$CP
   Q <- res$Q
   
   changepoints <- Changepoints <- Mus <- Sigmas <- jumprates <- CPnum <- NULL
 ################################################################################
 ################################################################################
 ################################################################################
################################################################################
 ################################################################################
   n <- length(obs)
   outliers <- rep(0, n)
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
   if ((max(diff(Uniformtimes)) < 5) && (length(Uniformtimes) <= 300) && (length(Uniformtimes) >=3)) flag <- TRUE
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
   if ((max(diff(Uniformtimes)) < 5) && (length(Uniformtimes) <= 300) && (length(Uniformtimes) >=3) ) flag <- TRUE
   }
   }
   print(j)
 
 ######################   ####################### 
 ###############   Trimming    ####################  
 cp <- res$cp
 m <- res$Number
 N <- c(cp[1]-1, diff(cp), (n+1)-cp[m-1])
 mu <- sderror <-  rep(NA, m)
 cp <- c(1, cp, ( n+1))
 
 for (i in 1:(m)){
      mus <- sderrors <- temp <- NULL
      for(r in 1:6){
 s <- sample(cp[i]:(cp[i+1]-1), size=ceiling(N[i]/2), prob=NULL)
 mu[i] <- mean(obs[s])
 sderror[i] <- sd(obs[s])
 mus <- c(mus, mu[i])
 sderrors <- c(sderrors, sderror[i])
 #tag <- which(abs(obs[cp[i]:(cp[i+1]-1)] - param[i])>2)
 index.temp <- cp[i]:(cp[i+1]-1)
#index.temp <- index.temp[-tag]
 FUNC <- function(xx){
         return(dnorm(xx, mean = mu[i], sd = sderror[i], log = TRUE))
                  }
 logL <- sum(sapply(obs[index.temp], FUN=FUNC))
 temp <- c(temp, logL)
 }
 k <- which.max(temp)
 mu[i] <- mus[k]
 sderror[i] <- sderrors[k]
}
 
 cp <- c(res$cp,n)
 trim <- rep(NA, n)
 trim[n] <- (obs[n]-mu[m])/sderror[m]
 k <- 1
 if (m>1){
 for (i in 1:n){
 if (i < cp[k]) trim[i] <- (obs[i]-mu[k])/sderror[k] else {
     k <- k+1
     if(k <= m) trim[i] <- (obs[i] -mu[k])/sderror[k]  else stop
     }
 } 
 } else {
         s <- sample(1:n, size=n/2, prob=NULL)
         trim <- (obs - mean(obs[s]))/sd(obs[s])
        }
 
 #theta <- param 
 index <- which(abs(trim) >=2)
 outliers[index] <- outliers[index] + 1
 ################################################################################
 
   res <- FFBS(Uniformtimes, JumpProb, eventtimes[-index], obs[-index], Mu0, Kappa0, Nv0, Sigma0)
   num <- res$Number
   cp <- res$cp
   CP <- res$CP
   Q <- res$Q 
   Mu <- res$Mu
   Sigma <- res$Sigma
   KK <- length(CP)
   if ((ending == CP[KK-1])&&(res$Number != 1)){
   CP <- CP[-KK]
   cp <- cp[-length(cp)]
   Q <- Q[,-(KK-1)]
   Q <- Q[-(KK-2),]
   num <- num - 1
   }  
   theta <- res$theta
   if (j > burnin){
   changepoints <- c(changepoints, list(cp))
   Changepoints <- c(Changepoints, list(CP))
   jumprates <- c(jumprates, list(Q))
   CPnum <- c(CPnum, num)
   Mus <- c(Mus, list(Mu))
   Sigmas <- c(Sigmas, list(Sigma))
         }
   } 
