## ranked probability score
rps <- function(obs, pred, baseline = NULL){

## pred is a matrix.  each column the prob of a given outcome
## obs is a vector with the number of the column that occured.  

nr <- nrow(pred)


################
id   <- is.finite(obs) & is.finite(apply(pred, 1, sum) )
obs  <- obs[id]
pred <- matrix(pred[id,], nrow = nr )
########


OBS <- matrix(0, nrow = length(obs), ncol = ncol(pred) )

## a loop, but it seems quick enough
for(i in 1:nrow(OBS) ){
OBS[i,obs[i] ]  <- 1  
}

OBS2 <- OBS


for(i in 1:ncol(OBS) ){
OBS2[,i] <- apply(matrix(OBS[, 1:i], nrow = nr), 1, sum)
}

PRED <- OBS
for(i in 1:ncol(pred) ){
PRED[,i] <- apply(matrix(pred[, 1:i], nrow = nr), 1, sum)
}



RPS <- mean( apply( ( PRED - OBS2)^2,1, sum)  )/ ( ncol(pred) -1 )
####

 if(is.null(baseline)){
xxx <- apply(OBS, 2, sum)/ nrow(OBS) ## avg occurrence of each event
   pred.climo <- matrix( xxx, nrow = nrow(OBS), ncol = ncol(OBS), byrow = TRUE)
} else{ ## climo should be a vector of the climo probs for each cat.
 pred.climo <- matrix( baseline, nrow = nrow(OBS), ncol = ncol(OBS), byrow = TRUE) 
}


PRED.climo <- OBS
for(i in 1:ncol(pred.climo) ){
PRED.climo[,i] <- apply(matrix(pred.climo[, 1:i], nrow = nr), 1, sum)
}


RPS.climo <- mean( apply( ( PRED.climo - OBS2)^2,1, sum)  )/ ( ncol(PRED.climo) -1 )

RPSS <- 1 - RPS/RPS.climo

return(list(rps = RPS, rpss = RPSS, rps.clim = RPS.climo ) )
  
}
