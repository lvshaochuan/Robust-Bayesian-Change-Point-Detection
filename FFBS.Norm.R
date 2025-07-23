FFBS <- function(Uniformtimes, JumpProb, eventtimes, obs){
n <- length(eventtimes)
k <- length(Uniformtimes)
filtering <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
filtering[1, ] <- c(1, rep(0,k))
#marginalLiks <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
changeprobs <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
survivals <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
JumpProb <- rbind(JumpProb, JumpProb[k,])
Uniformtimes <- c(Uniformtimes, ending)
for (i in 2:(k+1)){
     for (j in 1:i){
      if (j < i){
    cusumprob <- 0
## caculate the cumulative distribution function of arriving at changepoint
   for (t in j:(i-1)){
   changeprob <- 1
   for (r in j:t)
   changeprob <- changeprob*JumpProb[r,1]
   changeprob <- changeprob*JumpProb[r+1,2]
   cusumprob <- cusumprob + changeprob
   }

  survivala <- 1-cusumprob
  survivalb <- 1-cusumprob+changeprob 
  marginalLika <- MarginallogLik.Norm(eventtimes, obs, c(starting, 
  Uniformtimes)[j], Uniformtimes[i])     
  marginalLikb <- MarginallogLik.Norm(eventtimes, obs, c(starting, 
  Uniformtimes)[j], Uniformtimes[i-1])      
  changeprobs[i,j] <- changeprob
  survivals[i,j] <- survivalb
  filtering[i,j] <- filtering[i-1,j]*(survivala/survivalb)*exp(marginalLika-marginalLikb)       
       }
       
  if (j ==i){
  marginalLik <- MarginallogLik.Norm(eventtimes, obs, Uniformtimes[i-1], Uniformtimes[i])
  filtering[i,i] <- exp(marginalLik) * sum(filtering[i-1,1:(j-1)]*
  changeprobs[i,1:(j-1)]/survivals[i,1:(j-1)])  
       }     
  }
  filtering[i,] <- filtering[i,]/sum(filtering[i,])
 }

cpn <- cps <- NULL
cpn <- min(sample(1:(k+1), size=3, prob=filtering[k+1,], replace=T))
if (cpn != (k+1)) cps <- cpn
while (cpn > 1){
P <- rep(0, k+1)
index <- cpn-1
if(index==1) break
P[1:(index-1)] <- filtering[index,1:(index-1)]*changeprobs[index,1:(index-1)]/survivals[index,1:(index-1)]
P <- P/sum(P)
cpn <- min(sample(1:(index-1), size=3, prob=P[1:(index-1)], replace=T))
cps <- c(cpn, cps)
}

############################## Segment Parameter Update ########################

CP <- c(Uniformtimes[cps][-1])
m <- length(CP) + 1
cp <- NULL
j <- 1
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
#param <- rep(NA, m)
for (i in 1:(m-1)){
Q[i , i] <- -rgamma(1, alpha + 1, zeta+sojourn[i])
Q[i , (i+1)] <- -Q[i , i]
#param[i] <- mean(obs[cp[i]:(cp[i+1]-1)])
#param[i] <- rnorm(1, mean(obs[cp[i]:(cp[i+1]-1)]), sd=sigma/N[i])
##Lambda[i] <- rgamma(1, a + N[i], b + sojourn[i])
}
#param[m] <- mean(obs[cp[m]:n])
#param[m] <- rnorm(1, mean(obs[cp[m]:n]), sd=sigma/N[m])
##Lambda[m] <- rgamma(1, a + N[m], b + sojourn[m])
}
else{
m <- 1
cp <- NULL
Q <- NULL
#param <- rnorm(1, mean(obs), sd=sigma/n)
}
CP <- c(0, CP, eventtimes[n])
cat("\n Number:", m-1, "\n")
cat("Locations:", cp, "\n")
print(CP)
print(Q) 
return(list(Number=m, cp=cp, CP=CP, Q=Q, filtering=filtering))
}