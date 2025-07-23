###################   mean and sigma change points detection 

require(LaplacesDemon)

FFBS <- function(Uniformtimes, JumpProb, eventtimes, obs, Mu0, Kappa0, Nv0, Sigma0){
n <- length(eventtimes)
k <- length(Uniformtimes)
filtering <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
filtering[1, ] <- c(1, rep(0,k))

changeprobs <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
survivals <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
N <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
MuN <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
SigmaN <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
JumpProb <- rbind(JumpProb, JumpProb[k,])
Uniformtimes <- c(Uniformtimes, ending)
for (i in 2:(k+1)){
     for (j in 1:i){
      if (j < i){
    cusumprob <- 0
## caculate the cumulative distribution function of a changepoint
   for (t in j:(i-1)){
   changeprob <- 1
   for (r in j:t)
   changeprob <- changeprob*JumpProb[r,1]
   changeprob <- changeprob*JumpProb[r+1,2]
   cusumprob <- cusumprob + changeprob
   }

  survivala <- 1-cusumprob
  survivalb <- 1-cusumprob+changeprob 
  temp1 <- MarginallogLik.Norm2(eventtimes, obs, c(starting, 
  Uniformtimes)[j], Uniformtimes[i], Mu0, Kappa0, Nv0, Sigma0)
  marginalLika <- temp1$z     
  temp2 <- MarginallogLik.Norm2(eventtimes, obs, c(starting, 
  Uniformtimes)[j], Uniformtimes[i-1], Mu0, Kappa0, Nv0, Sigma0)
  marginalLikb <- temp2$z      
  changeprobs[i,j] <- changeprob
  survivals[i,j] <- survivalb
  filtering[i,j] <- filtering[i-1,j]*(survivala/survivalb)*exp(marginalLika-marginalLikb) 
  N[i,j] <- temp1$N
  MuN[i,j] <- temp1$MuN
  SigmaN[i,j] <- temp1$SigmaN     
       }
       
  if (j ==i){
  temp <- MarginallogLik.Norm2(eventtimes, obs, Uniformtimes[i-1], Uniformtimes[i], Mu0, Kappa0, Nv0, Sigma0)
  marginalLik <- temp$z
  filtering[i,i] <- exp(marginalLik) * sum(filtering[i-1,1:(j-1)]*
  changeprobs[i,1:(j-1)]/survivals[i,1:(j-1)])
  N[i,i] <- temp$N
  MuN[i,i] <- temp$MuN
  SigmaN[i,i] <- temp$SigmaN      
       }     
  }
  filtering[i,] <- filtering[i,]/sum(filtering[i,])
 }

cpn <- cps <- NULL
cpn <- min(sample(1:(k+1), size=2, prob=filtering[k+1,], replace=T))
if (cpn != (k+1)) cps <- cpn
while (cpn > 1){
P <- rep(0, k+1)
index <- cpn-1
if(index==1) break
P[1:(index-1)] <- filtering[index,1:(index-1)]*changeprobs[index,1:(index-1)]/survivals[index,1:(index-1)]
P <- P/sum(P)
cpn <- min(sample(1:(index-1), size=2, prob=P[1:(index-1)], replace=T))
cps <- c(cpn, cps)
}

############################## Parameter Update ########################

CP <- c(Uniformtimes[cps][-1])
m <- length(CP) + 1
cp <- NULL

j <- 1
logL <- 0
if (m>1){
for (i in 1:n){
if (eventtimes[i] > CP[j]){
    cp <- c(cp, i)
    j <- j+1
    if(j > length(CP)) break
    }
} 
N <- c(cp[1]-1, diff(cp), (n+1)-cp[m-1])
sojourn <- c(CP[1], diff(CP), eventtimes[n]-CP[m-1])
Q <- matrix(rep(0, m*m), nrow=m)
sigma2 <- rep(NA, m)
Mu <- rep(NA, m)
#param <- rep(NA, m)
cp <- c(1, cp)
for (i in 1:(m-1)){
Q[i , i] <- -rgamma(1, alpha + 1, zeta+sojourn[i])
Q[i , (i+1)] <- -Q[i , i]
Xbar <- mean(obs[cp[i]:(cp[i+1]-1)])
KappaN <- Kappa0 + N[i]
MuN <- (Kappa0*Mu0 + sum(obs[cp[i]:(cp[i+1]-1)]))/KappaN
NvN <- Nv0 + N[i]
SigmaN <- (Nv0*Sigma0^2 + sum((obs[cp[i]:(cp[i+1]-1)] - Xbar)^2) + N[i]*Kappa0*(Mu0 - Xbar)^2/KappaN)/NvN  
sigma2[i] <- rinvchisq(1, df=NvN, scale=sqrt(SigmaN))
Mu[i] <- rnorm(1, mean=MuN, sd=sqrt(SigmaN/KappaN))
logL <- logL + sum(dnorm(obs[cp[i]:(cp[i+1]-1)], mean=Mu[i], sd=sqrt(sigma2[i]), log=T))      
#param[i] <- mean(obs[cp[i]:(cp[i+1]-1)])
#param[i] <- rnorm(1, mean(obs[cp[i]:(cp[i+1]-1)]), sd=sigma/N[i])
##Lambda[i] <- rgamma(1, a + N[i], b + sojourn[i])
}
Xbar <- mean(obs[cp[m]:n])
KappaN <- Kappa0 + N[m]
MuN <- (Kappa0*Mu0 + sum(obs[cp[m]:n]))/KappaN
NvN <- Nv0 + N[m]
SigmaN <- (Nv0*Sigma0^2 + sum((obs[cp[m]:n] - Xbar)^2) + N[m]*Kappa0*(Mu0 - Xbar)^2/KappaN)/NvN  
sigma2[m] <- rinvchisq(1, df=NvN, scale=sqrt(SigmaN))
Mu[m] <- rnorm(1, mean=MuN, sd=sqrt(SigmaN/KappaN)) 
logL <- logL + sum(dnorm(obs[cp[m]:n], mean=Mu[m], sd=sqrt(sigma2[m]), log=T))
#param[m] <- mean(obs[cp[m]:n])
#param[m] <- rnorm(1, mean(obs[cp[m]:n]), sd=sigma/N[m])
##Lambda[m] <- rgamma(1, a + N[m], b + sojourn[m])
cp <- cp[-1]
}
else{
m <- 1
cp <- NULL
Q <- NULL
#param <- rnorm(1, mean(obs), sd=sigma/n)
Xbar <- mean(obs)
KappaN <- Kappa0 + n
MuN <- (Kappa0*Mu0 + sum(obs))/KappaN
NvN <- Nv0 + n
SigmaN <- (Nv0*Sigma0^2 + sum((obs - Xbar)^2) + n*Kappa0*(Mu0 - Xbar)^2/KappaN)/NvN  
sigma2 <- rinvchisq(1, df=NvN, scale=sqrt(SigmaN))
Mu <- rnorm(1, mean=MuN, sd=sqrt(SigmaN/KappaN))      
logL <- logL + sum(dnorm(obs, mean=Mu, sd=sqrt(sigma2), log=T))
}
CP <- c(0, CP, eventtimes[n])
cat("\n Number:", m-1, "\n")
cat("Locations:", cp, "\n")
print(CP)
print(Q) 
return(list(Number=m, cp=cp, CP=CP, Q=Q, filtering=filtering, Sigma=sigma2, Mu=Mu, logL=logL))
}