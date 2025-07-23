require(ggdist)
library(LaplacesDemon)
################################################################################
logLik.t <- function(eventtimes, obs, DF, sim_mu, starts, ends){
z <- 0
i <- 1
Sigma <- Var
###########  <=
while(eventtimes[i] <= starts) i <- i+1
########### <
while (eventtimes[i] < ends){
z <- z + dstudent_t(obs[i], df=DF, mu=sim_mu, log=TRUE)
i <- i+1
if (i > length(eventtimes)) break
  }
return(z)
}

FFBS.t <- function(Uniformtimes, JumpProb, eventtimes, obs, param, DF, Tau, Ui, Alpha){
n <- length(obs)
k <- length(Uniformtimes)
filtering <- matrix(rep(0, (k+1)*(k+1)), ncol=k+1)
filtering[1,] <- c(1, rep(0,k))

JumpProb <- rbind(JumpProb, JumpProb[k,])
Uniformtimes <- c(Uniformtimes, ending)

for (i in 2:(k+1)){
   for (j in 1:(i-1))  filtering[i,j+1] = filtering[i-1, j+1]*JumpProb[i,1] + filtering[i-1, j]*JumpProb[i,2]
   filtering[i,1] = filtering[i-1,1]*JumpProb[i,1]
   ff <- logLik.t(eventtimes, obs, DF=3, param, Uniformtimes[i-1], Uniformtimes[i])
#  ff <- logLik.Norm(eventtimes, obs, param, Sigma, Uniformtimes[i-1], Uniformtimes[i])
   filtering[i,] = log(filtering[i,]) + ff
   filtering[i,] = exp(filtering[i,] - max(filtering[i,]))
   filtering[i,] <- filtering[i,]/sum(filtering[i,])
   }
######  Sampling latent states   
S <- rep(1, k+1)
S[k+1] = min(sample(k+1, size=5, prob=filtering[k+1,], replace=T))
if(S[k+1] == 1) S[1:k] = 1

for (i in k:2){
  Prob <- rep(0, 2)
  active = c(S[i+1]-1, S[i+1])
  scaling <- filtering[i, active[1]] + filtering[i, active[2]] 
  Prob[2] = filtering[i, active[2]]/scaling * JumpProb[i+1, 1]
  Prob[1] = filtering[i, active[1]]/scaling * JumpProb[i+1, 2]
  S[i] = min(sample(active, size=5, prob=Prob, replace=T))
  if(S[i]==1){
              S[1:(i-1)] = 1
              break
             }
  }
cat(S, "\n")

CP <- NULL
CP <- Uniformtimes[which(diff(S)==1)]
m <- length(CP)
cp <- NULL
j <- 1

mu <- rep(NA, m+1)
#Ui <- rep(NA, n)

if(m > 0){
for (i in 1:n){
     if(eventtimes[i] >= CP[j]){
     cp <- c(cp, i)
     j <- j + 1
     if(j > length(CP)) break
     }
     }
       
N <- c(cp[1]-1, diff(cp), (n+1)-cp[m])
sojourn <- c(CP[1], diff(CP), eventtimes[n]-CP[m])
cp <- c(1, cp)
Q <- matrix(rep(0, (m+1)*(m+1)), nrow=m+1)
#theta <- rep(NA, m+1)
#precision <- rep(NA, m+1)

temp0 <- rep(NA, n)

for(i in 1:m){
Q[i, i] <- -rgamma(1, alpha + 1, zeta + sojourn[i])
Q[i, (i+1)] <- -Q[i,i]
#precision[i] <- 1/sigmaprior^2 + N[i]/Sigma^2
#theta[i] = rnorm(1, (muprior/sigmaprior^2 + N[i]*mean(obs[cp[i]:(cp[i+1]-1)])/Sigma^2)/precision[i], sd=1/precision[i])

temp1 <- temp2 <- 0

for (j in cp[i] : (cp[i+1] - 1)){
temp1 <- temp1 + 1/Ui[j]
temp2 <- temp2 + obs[j]/Ui[j]
}
mu[i] <- rnorm(1, temp2/temp1, Alpha/temp1)
for (j in cp[i] : (cp[i+1] - 1)){
Ui[j] <- rinvchisq(1, df=DF + 1, scale=(DF*Tau + (obs[j] - mu[i])^2/Alpha )/(DF + 1))
temp0[j] <- (obs[j] - mu[i])^2/Ui[j]
}
}

temp1 <- temp2 <- 0
for (j in cp[m+1] : n){
temp1 <- temp1 + 1/Ui[j]
temp2 <- temp2 + obs[j]/Ui[j]
}
mu[m+1] <- rnorm(1, temp2/temp1, Alpha/temp1)

for (j in cp[m+1] : n){
Ui[j] <- rinvchisq(1, df=DF + 1, scale=(DF*Tau + (obs[j] - mu[m+1])^2/Alpha )/(DF + 1))
temp0[j] <- (obs[j] - mu[m+1])^2/Ui[j]
}


} else{
       cp <- NULL
       Q <- NULL
#      theta <- rnorm(1, mean(obs), sd=Sigma/n)
       mu <- rnorm(1, sum(obs/Ui)/sum(1/Ui), Alpha/sum(1/Ui))
       for (k in 1:n) Ui <- rinvchisq(1, df=DF + 1, scale=(DF*Tau + (obs[k] - mu)^2/Alpha )/(DF + 1))
       temp0 <- (obs-mu)^2/Ui
      }
Tau <- rgamma(1, shape = n*DF/2, rate = DF*sum(1/Ui)/2)
Alpha <- rinvchisq(1, df=n, scale= mean(temp0))
      
CP <- c(0, CP, eventtimes[n])
cat("\n Number:", m, "\n")
cat("Locations:", cp, "\n")
#cat("Parameters:", theta, "\n")  
return(list(Number=m+1, cp=cp, CP=CP, Q=Q, states= S, filtering=filtering, param=mu, Tau=Tau, Ui=Ui, Alpha=Alpha)) 

}
