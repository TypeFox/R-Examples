value<- function(obs, pred = NULL, baseline = NULL,
                 cl = seq(0.05, 0.95, 0.05), plot = TRUE, all = FALSE,
                 thresholds = seq(0.05, 0.95, 0.05), ylim = c(-0.05,
                 1), xlim = c(0,1), ...){
####
# obs
# pred = NULL
# baseline = NULL
# cl = seq(0.05, 0.95, 0.05)
# plot = TRUE
# all = FALSE
#                 thresholds = seq(0.05, 0.95, 0.05), ylim = c(-0.05,
#                 1), xlim = c(0,1
####
  
if(!is.null(pred) ){
id   <- is.finite(obs) & is.finite(pred)
obs  <- obs[id]
pred <- pred[id] }else{
obs<- obs[is.finite(obs)]
}


### 2 by 2 contingency table

if(is.null(pred) & length(obs) ==4 ){
  print(" Assume data entered as c(n11, n01, n10, n00) Obs*Forecast")
n <- sum(obs)
a <- obs[1]
b <- obs[2]
c <- obs[3]
d <- obs[4]

F <- b/(b+d) ## CONDITIONAL FALSE ALARM RATE
H <- a/(a+c) ## CONDITIONAL HIT RATE

  if(is.null(baseline)){s <- (a+c)/n; baseline.tf <- FALSE}else
{s<- baseline; baseline.tf <- TRUE}

cl <- sort(c(cl, s) )## always add s to list of cl to ensure peak is used.

V1<- (1-F) - s/(1-s)*(1-cl)/cl *(1-H)
V2 <- H - (1-s)/s*cl/(1-cl)*F
V <- numeric(length(cl) )
V[cl < s]  <- V1[cl < s]
V[cl >= s] <- V2[cl >= s]

V<- matrix(V, ncol = 1)
Vmax       <- H - F ## Pierce skill score
positive   <- c(c/(c+d), a/(a+b) ) ## range of positive skill score
type <- "binary"
} else { ## prediction a vector of binary forecasts
#####################################################
#####################################################
if(prod( unique(pred)%in%c(0,1) )) { # prediction are a vector binary
  
if(is.null(baseline)){s <- mean(obs); baseline.tf <- FALSE}else
{s<- baseline; baseline.tf <- TRUE}

cl <- sort(c(cl, s) )## always add s to list of cl

F <- numeric()
H <- numeric()
Vmax <- numeric()
V <- matrix(nrow = length(cl), ncol = 1 )

A <- table(data.frame(obs = obs,pred = pred) )
a <- A[2,2]
b <- A[1,2]
c <- A[2,1]
d <- A[1,1]
n <- a + b + c + d
F[1] <- b/(b+d) ## FALSE ALARM RATE
H[1] <- a/(a+c) ## HIT RATE


V1<- (1-F[1]) - s/(1-s)*(1-cl)/cl *(1-H[1])
V2 <- H[1] - (1-s)/s*cl/(1-cl)*F[1]
VV <- numeric(length(cl) )
VV[cl < s] <- V1[cl < s]
VV[cl >= s]<- V2[cl >= s]
V[,1]<- VV
Vmax[1]     <- H[1] - F[1] #
positive <- c(c/(c+d), a/(a+b) ) ## range of positive skill score
} else {## close binary vector, open probabilistic forecast

## check?
if(max(pred)>1 | min(pred)<0 ) {
    ## predictions are a vector of probabilities
    stop("Predictions outside [0,1] range.  \n I am a bit confused. \n")
}

## make matrix with prediction for all thresholds

if(is.null(baseline)) {
    s <- mean(obs)
    baseline.tf <- FALSE
} else {
    s<- baseline
    baseline.tf <- TRUE
}

cl <- sort(c(cl, s) )## always add s to list of cl
  NCOL <- length(thresholds)
  PRED<- matrix(NA, nrow = length(pred), ncol = NCOL )
  for(i in 1:NCOL) PRED[,i] <- pred > thresholds[i]
  

F    <- numeric()
H    <- numeric()
Vmax <- numeric()
V    <- matrix(nrow = length(cl), ncol = ncol(PRED) )
n    <- length(pred)

for(i in 1:ncol(PRED)){
    ## try inserted for when table is 2X2
    A    <- table(data.frame(obs, PRED[,i]) )
    a    <- try(A[2,2], silent = TRUE )
    b    <- try(A[1,2], silent = TRUE)
    c    <- try(A[2,1], silent = TRUE)
    d    <- try(A[1,1], silent = TRUE)
    
    
    if(class(a) == "try-error") a<- NA
    if(class(b) == "try-error") b<- NA
    if(class(c) == "try-error") c<- NA
    if(class(d) == "try-error") d<- NA
    
    F[i] <- b/(b+d) ## FALSE ALARM RATE
    H[i] <- a/(a+c) ## HIT RATE
    
    
    V1<- (1-F[i]) - s/(1-s)*(1-cl)/cl *(1-H[i])
    V2 <- H[i] - (1-s)/s*cl/(1-cl)*F[i]
    VV <- numeric(length(cl) )
    VV[cl < s] <- V1[cl < s]
    VV[cl >= s]<- V2[cl >= s]
    V[,i]<- VV
    Vmax[i]     <- H[i] - F[i] #
} #close for 1:ncol(PRED)
#V.ind <- V
#V <- apply(V, 1, max)  ### outer envelope
} ## close probablistic option.

} # Close second if else stmt.

if(plot){

    if(!all) V <- apply(V, 1, max)
    matplot(cl, V, type = "l", ylim = ylim, xlim = xlim, ... )
    if(all) lines(cl,apply(V, 1, max), lwd = 2)
    abline(h=0)
    abline(v = s, lty = 2, lwd = 0.4)
}  ## close if plot
 
(aa<- list(vmax = Vmax, V = V, F= F, H = H, cl = cl , s = s,
           n = n) )

invisible(aa)

}## close function

