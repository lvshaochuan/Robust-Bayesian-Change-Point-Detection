long <- c(170,175,180,180,173,170)
lat <- c(-43,-36,-36,-38,-45,-43)

Polygon Subset from the NZ Catalogue

Number of Events Selected = 23032 

      Polygon Vertices:
              longitude latitude
            1       170      -43
            2       175      -36
            3       180      -36
            4       180      -38
            5       173      -45
            6       170      -43

      Depth Range:  min   = 45 
                    max   = Inf 

  Magnitude Range:  min   = 3.5 
                    max   = Inf 

    Time Interval:  start = 01Jan1955 00:00:00.00 
                    end   = Inf 

 Number: 4 
Locations: 497 2610 6037 7549                                           i
[1] 0.0000000 0.4960393 2.6093418 6.0365900 7.5489823 8.0200000
           [,1]       [,2]       [,3]       [,4]      [,5]
[1,] -0.2905832  0.2905832  0.0000000  0.0000000 0.0000000
[2,]  0.0000000 -0.6040577  0.6040577  0.0000000 0.0000000
[3,]  0.0000000  0.0000000 -0.4939123  0.4939123 0.0000000
[4,]  0.0000000  0.0000000  0.0000000 -0.8199915 0.8199915
[5,]  0.0000000  0.0000000  0.0000000  0.0000000 0.0000000
 Number: 3 
Locations: 383 5985 7667 
[1] 0.0000000 0.3823864 5.9842443 7.6668157 8.0200000
           [,1]       [,2]       [,3]    [,4]
[1,] -0.4556743  0.4556743  0.0000000 0.00000
[2,]  0.0000000 -0.1589591  0.1589591 0.00000
[3,]  0.0000000  0.0000000 -0.4704100 0.47041
[4,]  0.0000000  0.0000000  0.0000000 0.00000
Iteration 995
################################################################################
magnitude <- events[,"magnitude"] - 3.5
################################################################################
  goodpoints <- rep(0, length(magnitude))
  require(MASS)
  CP <- res$CP
  Q <- res$Q
  theta <- res$theta
  changepoints <- Changepoints <- intensity <- jumprates <- CPnum <-  threshold <- NULL
  ################################################################################
  ##############################################################################
  ##############################################################################
  for (j in 1:iter){

  cat("\n iteration", j, "\n")

################################################################################
obs <- NULL
Ind <- NULL
cuts <- NULL
cpt <- z$cp
cpt <- c(0, z$cp, length(magnitude))
delta <- 0.1

for (i in 1:(length(cpt) - 1)){
mag <- magnitude[(cpt[i] + 1) : cpt[i+1]]
minmag <- min(mag, na.rm = TRUE)
counts <- hist(mag, plot = FALSE, breaks = seq(minmag, 
        (max(mag, na.rm = TRUE) + delta), delta), include.lowest = TRUE, 
        right = FALSE)
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
#        ylab <- expression(paste(log[10], "(Proportion of Events with ", 
#            Magnitude >= ~m, ")"))  
            temp <- lqs(x, y)
            bvalue <- -temp$coefficients[2]
#            abline(y[1] + minmag * bvalue, -bvalue, lty = 3)
MAD <- median(abs(temp$residuals - median(temp$residuals)))
IND <- which(abs(temp$residuals) > 2*MAD)
IND <- IND[IND <= median(1:length(x))]
cat("\n outliers:", IND, "\n")            
cutting <- max(c(0,IND))
if (cutting >= 8) cutting <- 0
tag <- which(mag >= 0.1 * cutting)
Ind <- c(Ind, cpt[i] + tag)
obs <- c(obs, mag[tag] - 0.1 * cutting)
cuts <- c(cuts, cutting)
   }
goodpoints[Ind] <- goodpoints[Ind] + 1  
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
  Lambda <- res$Lambda
  if (j > burnin){
  changepoints <- c(changepoints, list(cp))
  Changepoints <- c(Changepoints, list(CP))
  jumprates <- c(jumprates, list(Q))
  CPnum <- c(CPnum, num)
  intensity <- c(intensity, list(Lambda))
  threshold <- c(threshold, list(cuts))
        }
  
################################################################################
  } 