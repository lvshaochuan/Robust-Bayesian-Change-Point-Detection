MarginallogLik.Norm2 <- function(eventtimes, obs, starts, ends, Mu0, Kappa0, Nv0, Sigma0){
N <- 0
Cumsum <- 0
Squaresum <- 0

z <- 0
i <- 1
while(eventtimes[i] < starts) i <- i+1
while (eventtimes[i] <= ends){
  N <- N + 1
  Cumsum <- Cumsum + obs[i]
  Squaresum <- Squaresum + obs[i]^2
  i <- i+1
  if (i > length(eventtimes)) break
  }
Xbar <- Cumsum/N
SSum <- Squaresum - N*Xbar^2

KappaN <- Kappa0 + N
MuN <- (Kappa0*Mu0 + Cumsum)/KappaN
NvN <- Nv0 + N
SigmaN <- (Nv0*Sigma0^2 + SSum + N*Kappa0*(Mu0 - Xbar)^2/KappaN)/NvN                                   
z <- loggamma(NvN/2) - loggamma(Nv0/2) + 0.5*(log(Kappa0) - log(KappaN)) - N/2*log(pi) 
z <- z + 0.5*Nv0*log(Nv0*Sigma0^2) - 0.5*NvN*log(NvN*SigmaN)                              
return(list(z=z, N=N, MuN=MuN, SigmaN=SigmaN))                                 
}

loggamma <- function(u){
z <- 0
while (u>1) {
  z <- log(u-1) + z
  u <- u-1
  }
z <- z + log(gamma(u))
return(z)
}
