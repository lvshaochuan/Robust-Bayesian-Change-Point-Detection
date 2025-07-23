  CP <- res$CP
  Q <- res$Q
  theta <- res$theta
  changepoints <- Changepoints <- intensity <- jumprates <- CPnum <- NULL
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
###############   Trimming    ####################  
cp <- res$cp
m <- res$Number
N <- c(cp[1]-1, diff(cp), (n+1)-cp[m-1])

param <-  rep(NA, m)
cp <- c(1, cp, ( n+1))

for (i in 1:(m)){
     retain <- temp <- NULL
     for(r in 1:3){
param[i] <- mean(sample(obs[cp[i]:(cp[i+1]-1)], size=N[i]/2, prob=NULL))
retain <- c(retain, param[i])
tag <- which(abs(obs[cp[i]:(cp[i+1]-1)] - param[i])>2)
#index.temp <- cp[i]:(cp[i+1]-1)
#index.temp <- index.temp[-tag]
#logL <- sapply(obs[index.temp], dnorm(x, mean = param[i], sd = 1, log = TRUE))
temp <- c(temp, length(tag))
}
k <- which.min(temp)
param[i] <- retain[k]
}


cp <- c(res$cp,n)
trim <- rep(NA, n)
trim[n] <- obs[n]-param[m]
k <- 1
if (m>1){
for (i in 1:n){
if (i < cp[k]) trim[i] <- obs[i]-param[k] else {
    k <- k+1
    if(k <= m) trim[i] <- obs[i] -param[k]  else stop
    }
} 
} else {trim <- trim-mean(sample(obs, size=n/2, prob=NULL))}

theta <- param 
index <- which(abs(trim) >=2)


  
######################   ####################### 
  res <- FFBS(Uniformtimes, JumpProb, eventtimes[-index], obs[-index])
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
        }
  }      
    