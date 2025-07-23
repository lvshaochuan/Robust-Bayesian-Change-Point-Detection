MarginallogLik.Norm <- function(eventtimes, obs, starts, ends){
N <- 0
Cumsum <- 0
Squaresum <- 0
logNormConst <- 0
logkernel <- 0
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
logNormConst <- (1-N)*log(sigma)-0.5*N*log(2*pi)-0.5*log(N*sigmaprior^2+sigma^2)
logkernel <-  -0.5*Squaresum/(sigma^2)-0.5*muprior^2/sigmaprior^2
logkernel <- logkernel + 0.5*(Cumsum/sigma^2+muprior/sigmaprior^2)^2/(N/sigma^2+1/sigmaprior^2)
z <- logNormConst + logkernel   
return(z)
}