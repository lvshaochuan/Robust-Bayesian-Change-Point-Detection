Uniform <- function(jumptimes, rates, starts, ends, pointer, rate){
##### Simulate Uniformisation times #####
initial <- starts
Poissontimes <- NULL
tag <- NULL
pointer <- c(pointer, ends)
K <- length(pointer)
i <- j <- 1
while(initial <= ends){
if (initial < pointer[j]) {
Poissontime <-  rexp(1, rate[i])
initial <- initial + Poissontime
if(initial <= pointer[j]) {
   Poissontimes <- c(Poissontimes, initial)
   tag <- c(tag, i)
   } else{
          if (initial >= ends) break
          else {
          i <- i + 1
          initial <- pointer[j]
          j <- j + 1
          }
          }
   } 
 }
##### Tagging jumptimes #####
N <- length(jumptimes) - 2
tag0 <- rep(NA, N)

if(N > 0){
for (i in 1:N)
tag0[i] <- sum((jumptimes[i + 1] - c(starts, pointer[-K])) >= 0)
}

Poissontimes <- c(Poissontimes, jumptimes[-c(1,length(jumptimes))])
tag <- c(tag, tag0)
temp <- order(Poissontimes)
tag <- tag[temp]
Poissontimes <- Poissontimes[temp] 
#####  Transition Probability Matrix  #####
n <- length(Poissontimes)
m <- length(tag0)
mark <- rep(1/n, n)
delimit <- NULL
if(m != 0){
for(i in 1:m){
  j = which(Poissontimes == jumptimes[i+1])
  delimit <- c(delimit, j)
  }
mark[1:delimit[1]-1] = rates[1] 
if(m >= 2){
   for(i in 2:m) mark[delimit[i-1]:(delimit[i]-1)] = rates[i]}
#mark[delimit[m]:n] = rates[m + 1]
# The last segment has no more than 1 change point 
}
Rate <- rep(NA, n)
for (i in 1:length(rate)){
Rate[tag == i] = rate[i]
}
JumpProb <- matrix(rep(NA, 2*n), ncol=2)
for(i in 1:n){
JumpProb[i,2] = mark[i]/(Rate[i] + mark[i])
JumpProb[i,1] = 1- JumpProb[i,2]
}
return(list(U=Poissontimes, JumpProb=JumpProb, tag=tag))
}
