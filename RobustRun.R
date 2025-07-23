  CP <- res$CP
  Q <- res$Q
  theta <- res$theta
  changepoints <- Changepoints <- intensity <- jumprates <- CPnum <- Inds <- NULL
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
  n <- length(obs)
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
  if ((max(diff(Uniformtimes) < 5)) && (length(Uniformtimes) <= 300)  && (length(Uniformtimes) >=3)) flag <- TRUE
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
  if ((max(diff(Uniformtimes) < 5)) && (length(Uniformtimes) <= 300)  && (length(Uniformtimes) >=3)) flag <- TRUE
  }
  }
  print(j)
###############   Trimming    ####################  
cp <- res$cp
m <- res$Number
N <- c(cp[1]-1, diff(cp), (n+1)-cp[m-1])
param <-  rep(NA, m)
cp <- c(1,cp)
for (i in 1:(m-1)){
param[i] <- mean(sample(obs[cp[i]:(cp[i+1]-1)], size=N[i]/2, prob=NULL))
}
param[m] <- mean(sample(obs[cp[m]:n], size=N[m]/2, prob=NULL))

cp <- c(res$cp,n)
trim <- rep(NA, n)
trim[n] <- ob



if (m>1){
for (i in 1:n){
if (i < cp[k]) trim[i] <- obs[i]-param[k] else {
    k <- k+1
    if(k <= m) trim[i] <- obs[i] -param[k]  else stop
    }
} 
} else {trim <- trim-mean(sample(obs, size=n/2, prob=NULL))}

theta <- param 
Index <- which(abs(trim) >=2)


  
######################   ####################### 
  res <- FFBS(Uniformtimes, JumpProb, eventtimes[-Index], obs[-Index])
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
  theta <- res$theta
  if (j > burnin){
  changepoints <- c(changepoints, list(cp))
  Changepoints <- c(Changepoints, list(CP))
  jumprates <- c(jumprates, list(Q))
  CPnum <- c(CPnum, num)
  intensity <- c(intensity, list(theta))
  Inds <- c(Inds, list(Index))
        }
  }  
  
  
################################################################################
################################################################################
################################################################################
 position <- NULL
 for (i in 2001:length(CPnum)){
 temp <- length(Changepoints[[i]])
 position <- c(position, Changepoints[[i]][-c(1,temp)])
 }
 

 savpar <- par(mfrow=c(3,1))
 plot(obs, xlab="", ylab="", xlim=c(0, 4050), main="(a)")
 points(Index, obs[Index], col="red", pch=10)
 hist(position/eventtimes[1], 1000, xlim=c(0, 4050), xlab="Changepoint Location", main="(b) Changepoint Locations")
 box()
 hist(CPnum[2001:3000]-1, 100, xlab="Number of Changepoints", main="(c) Number of Changepoints", col="grey")
 box()
 par(savpar) 
    
    