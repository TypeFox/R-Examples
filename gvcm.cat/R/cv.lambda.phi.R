cv.lambda <-
function(
x, 
y, 
weights, 
family,
control, 
acoefs, # a.coefs element
lambda, # TRUE or vector for cv
phis,
weight, # weight.const
start,
offset,
L.index, T.index, 
indices, 
...
)

{

# lambda vector
if(!is.logical(lambda)){
  cross <- cv.vectors(x, y, weights, family, control, acoefs, lambda=lambda, phis,
           weight, start, offset, L.index, T.index, indices)
  lambda <- cross$lambda
  lambdas <- cross$lambdas
  scores <- cross$score
  scores.sd <- cross$score.sd
  colnames(scores) <- as.character(lambdas)
  coefs <- cross$coefs
}

# lambda==TRUE
if(is.logical(lambda)){

  # definitions
  p <- ncol(x)
  basis <- exp(1/4*log(control$lambda.upper-control$lambda.lower))

  scores <- matrix(nrow=1, ncol=0)
  scores.sd <- matrix(nrow=1, ncol=0)
  coefs <- matrix(nrow=nrow(acoefs$A), ncol=0)
  lambdas <- c()
  closer <- if (basis >1) {round(control$lambda.lower + c(.1, basis^(0:4)), 2)} else 
#  closer <- if (basis >1) {round(control$lambda.lower + c(basis^(0:4)), 2)} else 
         {seq(control$lambda.lower, control$lambda.upper, by=control$lambda.accuracy)}
  i <- 1

  while (max(closer)-min(closer)>=control$lambda.accuracy && i <11){

      closer <- closer[which(closer<=control$lambda.upper)]
      closer <- unique(closer)
      if (any(closer %in% lambdas)) closer <- closer[-which(closer %in% lambdas)]
      
      if (length(closer)>0){
          cross <- cv.vectors(x, y, weights, family, control, acoefs, lambda=closer, phis,
                   weight, start, offset, L.index, T.index, indices)

          scores <- cbind(scores, cross$score)
          scores.sd <- cbind(scores.sd, cross$score.sd)
          lambdas <- c(lambdas, cross$lambdas)
          lambda <- max(lambdas[which(scores==min(scores))])[1]
          coefs <- cbind(coefs, cross$coefs)
      }

      j <- log(lambda-control$lambda.lower)/log(basis)
      closer <- round(control$lambda.lower + basis^(c((j-2^(-i)), (j+2^(-i)))),digits=2)
      i <- i+1

  }

##
  if (control$lambda.lower==0) {
      closer <- seq(.05, .5, length.out=10)
      if (any(closer %in% lambdas)) closer <- closer[-which(closer %in% lambdas)]
  
      if (length(closer)>0){
          cross <- cv.vectors(x, y, weights, family, control, acoefs, lambda=closer, phis,
                   weight, start, offset, L.index, T.index, indices)
  
          scores <- cbind(scores, cross$score)
          scores.sd <- cbind(scores.sd, cross$score.sd)
          lambdas <- c(lambdas, cross$lambdas)
          lambda <- max(lambdas[which(scores==min(scores))])[1]
          coefs <- cbind(coefs, cross$coefs)
      }
  }
##

  coefs  <- coefs[,order(as.numeric(colnames(coefs)))]
  scores <- scores[,order(as.numeric(colnames(scores)))]
  scores.sd <- scores.sd[,order(as.numeric(colnames(scores.sd)))]

  if (control$tuning.criterion=="1SE") {
       which.lambda <- which(lambdas==lambda)
       lambda.sd <- scores.sd[which.lambda]
       lambda.can <- which( lambdas[which(scores <= (scores[which.lambda] + lambda.sd))] >= lambda )
       lambda <-  max(lambdas[lambda.can]) # if (length(lambda.can)>0) max(lambdas[lambda.can]) else lambda
      }
 
}
    
return(list(lambda=lambda, score=scores, coefs=coefs, score.sd=scores.sd))
}

