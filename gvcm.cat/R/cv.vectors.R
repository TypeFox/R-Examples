cv.vectors <-
function(
x, 
y, 
weights, 
family,
control, 
acoefs, # acoefs object
lambda, # vector for cv
#phi, # vector for cv
phis, # grouped.fused, vs, spl, 
weight, # weight.const
start,
offset,
L.index, T.index, 
indices, 
...
)

{
n <- nrow(x)
which.a <- acoefs$which.a
acoefs <- acoefs$A

losses <- matrix(ncol=length(lambda), nrow=1)
losses.sd <-  matrix(ncol=length(lambda), nrow=1) 
colnames(losses) <- colnames(losses.sd) <- as.character(lambda) #####
coefs <- matrix(ncol=length(lambda), nrow=nrow(acoefs))
colnames(coefs) <- as.character(lambda)

# if "GCV", "UBRE"
if (control$tuning.criterion %in% c("GCV", "UBRE")) {

    if(control$tuning.criterion == "GCV"){
        crit <- function(n, dev, rank, sc=1) {n * dev/(n-rank)^2} # dev + rank von fit!!      
        }
    if(control$tuning.criterion == "UBRE"){
        crit <- function(n, dev, rank, sc=1) {dev/n + 1*rank*sc/n - sc}
        }
        
    if(control$cv.refit == FALSE){  
       evalcv <- function(x, coefficients, control, y, weights, dev, rank) {
            output <- list(deviance=dev, rank=rank)
            return(output)
            }   
       } else {
       evalcv <- function(x, coefficients, control, y, weights, dev, rank){
            reductionC <- reduce(coefficients, indices, control$assured.intercept, control$accuracy)$C
            x.reduced <- as.matrix(x %*% reductionC)     
            try(refitted <- glm.fit(x.reduced, y, weights, family=family,
                intercept = FALSE))
            dev.refitted <- refitted$deviance
            rank.refitted <- refitted$rank        
            output <- list(deviance=dev.refitted, rank=rank.refitted)
            return(output)
            }
       }
    
    for (i in 1:length(lambda)) {
            suppressWarnings(model <- gvcmcatfit(x, y, weights=weights, family,
                  control, acoefs, lambda=lambda[i], phis=phis, weight, 
                  which.a, start = start, offset = offset))
            eval.cv <- evalcv(x, model$coefficients, control, y, weights, model$deviance, model$rank)
            losses[1,i] <- crit(n, eval.cv$deviance, eval.cv$rank, sc=1)
            coefs[,i] <- model$coefficients
    }    
}
   
# if K-fold/deviance
if(control$tuning.criterion %in% c("deviance", "1SE")){

if(family$family == "binomial"){
  l <- function(y,mudach,weights=weights){sum((y*log(mudach) + (1-y)*log(1-mudach))*weights + lchoose(weights,y*weights))}
  d <- function(y,mudach,weights=weights){e <- matrix(0, nrow=length(y), ncol=1)
   rein.binaere <- if (any(y==0) || any(y==1)) c(which(y==0),which(y==1)) else -(1:length(y))
   e[rein.binaere,] <- (-2*log(1-abs(y-mudach))*weights)[rein.binaere]
   e[-rein.binaere,] <- (weights*(y*log(y/mudach)+(1-y)*log((1-y)/(1-mudach))))[-rein.binaere]
   return(e)}
  }
if(family$family == "gaussian"){
  l <- function(y,mudach,weights=weights){-1/2 * sum(weights*(y - mudach)^2) - n*log(sqrt(2*pi))}
  d <- function(y,mudach,weights=weights){weights*(y-mudach)^2}   #
  }
if(family$family == "poisson"){
  l <- function(y,mudach,weights=weights){sum(weights*((y*log(mudach)) - mudach - lgamma(y+1)))}
  d <- function(y,mudach,weights=weights){e <- matrix(0, nrow=length(y), ncol=1)
   e[which(y==0),] <-(2*mudach*weights)[which(y==0)]
   e[which(y!=0),] <- (2*weights*((y*log(y/mudach))+mudach-y))[which(y!=0)]
   return(e)}
  }
if(family$family == "Gamma"){
  l <- function(y,mudach,weights=weights){sum(weights*(log(mudach) - y/mudach))}
  d <- function(y,mudach,weights=weights){2*(-log(y) + log(mudach) + y/mudach -1)*weights}
  }
dev.res <- function(y,mudach,weights) {(y-mudach)/abs(y-mudach) * sqrt(d(y, mudach, weights)) }
dev <- function(y,mudach,weights) {sum(d(y=y,mudach=mudach, weights))}

    crit <- dev
    if(control$cv.refit == FALSE){  
       evalcv. <- function(x.tr, x.te, coefficients, control, training.y, training.weights) {
            output <- list(coefficients=coefficients, test.x=x.te)
            return(output)
            }   
       } else {
       evalcv. <- function(x.tr, x.te, coefficients, control, training.y, training.weights){
            reductionC <- reduce(coefficients, indices, control$assured.intercept, control$accuracy)$C
            x.reduced.tr <- as.matrix(x.tr %*% reductionC)     
            x.reduced.te <- as.matrix(x.te %*% reductionC) 
            try(beta.refitted <- glm.fit(x.reduced.tr,training.y, training.weights, family=family,
                intercept = FALSE)$coefficients)
            output <- list(coefficients=beta.refitted, test.x=x.reduced.te)
            return(output)
           }
       }

    for (i in 1:length(lambda)) {

        loss <- c()
        for (j in 1:control$K){
  
            training.x <- x[L.index[[j]],]  
            training.y <- y[L.index[[j]]]
            training.weights <- weights[L.index[[j]]]
            test.x <- x[T.index[[j]],]
            test.y <- y[T.index[[j]]]
            test.weights <- weights[T.index[[j]]]
            
            if(control$standardize){
               if (sum(x[,1])==n || sum(x[,1])==sum(weights)) {
                     scaling <- c(1, sqrt(colSums(t((t(training.x[,-1]) - colSums(diag(training.weights)%*%training.x[,-1])/sum(training.weights))^2*training.weights))/(sum(training.weights)-1)))
                   } else {
                     scaling <- sqrt(colSums(t((t(training.x)      - colSums(diag(training.weights)%*%training.x)     /sum(training.weights))^2*training.weights))/(sum(training.weights)-1))
                   }
               training.x <- scale(training.x, center = FALSE, scale = scaling)
               test.x <- scale(test.x, center = FALSE, scale = scaling)
            }

            suppressWarnings(model <- gvcmcatfit(training.x, training.y, weights=training.weights, family,
                  control, acoefs, lambda=lambda[i], phis=phis, weight, 
                  which.a, start = start, offset = offset))
    
            eval.cv <- evalcv.(training.x, test.x, model$coefficients, control, training.y, training.weights)
            test.mudach <- family$linkinv(eval.cv$test.x %*% eval.cv$coefficients)                            
            
            loss <- c(loss, crit(test.y, test.mudach, test.weights))
        }

        losses[1,i] <- sum(loss) / control$K
        losses.sd[1,i] <- sd(loss) * (control$K-1) / control$K  #######
    }
    coefs <- NA
} # end deviance

# Ergebnis
 opt <- max(lambda[which(losses==min(losses))])[1] # max(lambda[(which(losses==min(losses))-1)%/%dim(losses)[1]+1])[1]
      
return(list(
lambda=opt, 
lambdas=lambda,
score=losses, 
coefs=coefs, 
score.sd = losses.sd
))

}

