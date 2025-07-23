MarginallogLik.Exp <- function(evtimes, obs, initial, last){
  Num <- 0
  Cumsum <- 0
  n <- length(evtimes)
  i <- 1
  if ((evtimes[i] < initial) && (evtimes[i] < evtimes[n])){
  while ((evtimes[i] < initial) && (evtimes[i] < evtimes[n])) i <- i+1
  }
  while (evtimes[i] <= last) {
  Num <- Num + 1
  Cumsum <- Cumsum + obs[i]
  i <- i + 1
  if (i > length(evtimes)) break
  }
  if (Num > 0) logL <- a*log(b)-loggamma(a) + loggamma(a+Num) - (a+Num)*log(b+Cumsum)
  else logL <- 0
  return(logL)
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

#loggamma <- function(u){
#if (u > 1)
#return(sum(log(1:(u-1))))
#}