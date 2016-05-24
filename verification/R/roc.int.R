####################################################### internal roc function
roc.int <- function(x, pred, thres,  binormal ){ # internal function
  # that returns plot points
thres <- c(thres[1] - 1, thres )
  
H     <- numeric()
F     <- numeric()
n        <- length(x)
n.thres  <- length(thres)  # number of unique thresholds, thres =  
lng      <- length(x)            # 

a         <- x > 0 # event happened
a.sum     <- sum(a) # n*1
a.not.sum <- sum(!a) # n*0

for(i in 1:(n.thres) ){
  
  b      <- pred >   thres[i] # predict if value is greater than
                              # decision thresholds
                              # add point c(1, 1)
 # browser()
  H[i]<- sum(  b * a  )/ a.sum ## hit rate
  F[i]<- sum(  b  * (!a)  )/ a.not.sum  ## False alarm rate
 
} ## close for loop 1:n.thres

    if(binormal){
  
       zH <- c( qnorm( H ) )  # NA are for top and bottom value
       zF <- c( qnorm(F) )# NA are for top and bottom value
       } else {
         zH <- rep(NA, n.thres)
         zF <- rep(NA, n.thres)
       } ## close if binormal

return(cbind(thres, H, F, zH, zF))
       }    # close roc function

##############################################################
