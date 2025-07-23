  CP <- res$CP
  Q <- res$Q
  theta <- res$theta
  changepoints <- Changepoints <- intensity <- jumprates <- CPnum <- cutoffs <- NULL
################################################################################
################################################################################
################################################################################
for (j in 1:iter){
obs <- NULL
Ind <- NUll
cuts <- NULL
cpt <- c(1, z$cp, n)
delta <- 0.1 

for (j in 1:(length(cpt) - 1)){
mag <- magnitude[cpt[j] : (cpt[j+1] - 1)]
minmag <- min(mag, na.rm = TRUE)
counts <- hist(mag, plot = F, breaks = seq(minmag, 
        (max(mag, na.rm = TRUE) + delta), delta), include.lowest = TRUE, right = FALSE)
x <- counts$breaks
y <- 1 - cumsum(counts$density) * delta
y <- c(1, y)
ind <- counts$counts >= 10
x <- x[ind]
y <- y[ind]
ind <- which(is.na(x))
if (length(ind) != 0){
x <- x[-ind]
y <- y[-ind]
}
        
        aa <- (y > 1e-12)
        x <- x[aa]
        y <- log10(y[aa])
        temp <- lqs(x, y)
        bvalue <- -temp$coefficients[2]
        Lambda <- log(10) * bvalue
for (i in 1:10){
     P <-  (exp(-Lambda)-exp(-2*Lambda))/(exp(-Lambda*(i-1)*delta)-exp(-Lambda*i*delta) + exp(-Lambda)-exp(-2*Lambda)) 
     cat(P, "\n")
     Prob <- pnbinom(q=counts$counts[i], prob=P, size=sum(counts$counts[11:20]))
     cat(Prob, "\n")
     if (Prob >= 0.01) {
          cutting <- i
          break
          }
      }
tag <- which(mag >= 0.1 * cutting)
Ind <- c(Ind, cpt[i] + tag)
obs <- c(obs, mag[tag] - 0.1 * cutting)
cuts <- c(cuts, cutting)
   } 
cat("\n cutoff:", cuts, "\n")
######################   #######################
  if (res$Number != 1){  
  rates <- -diag(Q)
  rates <- rates[-length(rates)]
  rates <- c(rates, 1/(ending-starting))
  flag <- FALSE
  while (flag==FALSE){
  A <- Uniform(CP, rates)
  Uniformtimes <- A$U
  index <- which(diff(c(starting,Uniformtimes,ending)) <= max(diff(eventtimes[Ind])))
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
  index <- which(diff(c(starting,Uniformtimes,ending)) <= max(diff(eventtimes[Ind])))
  Uniformtimes <- Uniformtimes[-index]
  JumpProb <- A$JumpProb[-index, ]
  if ((max(diff(Uniformtimes)) < 5) && (length(Uniformtimes) <= 300) && (length(Uniformtimes) >=3) ) flag <- TRUE
  }
  }
################################################################################ 
  evtimes <- eventtimes[Ind]
  z <- res <- FFBS(Uniformtimes, JumpProb, evtimes, obs)
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
#  theta <- res$theta
  if (j > burnin){
  changepoints <- c(changepoints, list(cp))
  Changepoints <- c(Changepoints, list(CP))
  jumprates <- c(jumprates, list(Q))
  CPnum <- c(CPnum, num)
  intensity <- c(intensity, list(theta))
  cutoffs <- c(cutoffs, list(cuts))
        }
  
################################################################################
}