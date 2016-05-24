################################################
#### AUTHOR:     Arnost Komarek             ####
####             (2005)                     ####
####                                        ####
#### FILE:       bayessurvreg.design.R      ####
####                                        ####
#### FUNCTIONS:  bayessurvreg.design        ####
################################################

### ======================================
### bayessurvreg.design
### ======================================
## Subfunction common for all 'bayessurvreg' functions
##  - extract a design information
##
## 18/03/2004
##
bayessurvreg.design <- function(m, formula, random, data, transform, dtransform)
{

   tempF <- c("", "formula", "data", "subset", "na.action")
   mF <- m[match(tempF, names(m), nomatch=0)]
   mF[[1]] <- as.name("model.frame")
   special <- c("cluster")
   TermsF <- if(missing(data)) terms(formula, special)
             else              terms(formula, special, data = data)
   mF$formula <- TermsF
   mF <- eval(mF, parent.frame())
 
   ## Response matrix
   Y <- model.extract(mF, "response")
   Yinit <- Y
   if (!inherits(Y, "Surv"))
      stop("Response must be a survival object. ")
   type <- attr(Y, "type")
   if (type == 'counting') stop ("Invalid survival type ('counting' is not implemented). ")
   n <- nrow(Y)
   nY <- ncol(Y)
   
   ## Cluster indicators   
   cluster <- attr(TermsF, "specials")$cluster
   dropx <- NULL
   if (length(cluster)) {
     tempc <- untangle.specials(TermsF, "cluster", 1:10)
     ord <- attr(TermsF, "order")[tempc$terms]
     if (any(ord > 1))
         stop("Cluster can not be used in an interaction")
     cluster <- strata(mF[, tempc$vars], shortlabel = TRUE)
     dropx <- tempc$terms
   }
   if (length(dropx)){
     newTermsF <- TermsF[-dropx]
     attr(newTermsF, "intercept") <- attr(TermsF, "intercept")  ## I do not why but the command on the previous row
                                                                 ## sets attr(newTermsF, "intercept") always to 1,
                                                                 ## irrespective of what attr(TermsF, "intercept") was...
   }
   else{
     newTermsF <- TermsF
   }      

   ## Design matrix for both fixed and random effects X
   ## (finally, always without the intercept)
     ## Temporarily, include always at least the intercept
     ##  (to get nice rownames and initial estimates)  
   attr(newTermsF, "intercept") <- 1
   Xinit <- model.matrix(newTermsF, mF)
   rnamesX <- row.names(Xinit)
   cnamesX <- colnames(Xinit)
   
     ## Finally, intercept will be always removed
   nX <- ncol(Xinit) - 1
   cnamesX <- cnamesX[-1]
   if (nX){
     X <- Xinit
     X <- X[,-1]                                      ## removal of the intercept
     attr(X, "assign") <- attr(Xinit, "assign")[-1]   ## removal of the intercept
     attr(X, "contrasts") <- attr(Xinit, "contrasts")
   }
   else{
     X <- NULL
   }     
   indb <- if (nX) rep(-1, nX)                          ## initially, all effects are only fixed
           else    0
      
   ## Design matrix for random effects 
   randomInt <- FALSE
   nrandom <- 0
   if (!missing(random)){
     if (!length(cluster)) stop ("You have to indicate clusters when you wnat to include some random effects. ")
     tempR <- c("", "random", "data", "subset", "na.action")
     mR <- m[match(tempR, names(m), nomatch=0)]
     mR[[1]] <- as.name("model.frame")
     names(mR)[2] <- "formula"
     TermsR <- if(missing(data)) terms(random)
               else              terms(random, data = data)
     lTR <- length(attr(TermsR, "variables"))
     if (lTR == 1 & !attr(TermsR, "intercept")){        ## do nothing, in reality no random terms in the model
       names.random <- character(0)
     }
     else{
       if (lTR == 1 & attr(TermsR, "intercept")){        ## the only random term is the intercept
         randomInt <- TRUE
         names.random <- character(0)
       }
       else{
         mR$formula <- TermsR
         mR <- eval(mR, parent.frame())
         if (attr(TermsR, "intercept")){
           randomInt <- TRUE
           names.random <- colnames(model.matrix(TermsR, mR))[-1]
           ## attr(TermsR, "intercept") <- 0     ## remove it from the design matrix of random effects
         }
         else{
           names.random <- colnames(model.matrix(TermsR, mR))
         }           
       }
       nrandom <- 1*randomInt + length(names.random)
       if (sum(names.random %in% cnamesX) != nrandom - 1*randomInt) stop("Each random effect has to have also its fixed counterpart.")
       find.indeces <- function(all.eff){     ## Find indeces of random effects in a design matrix to be passed to C++
         where <- names.random %in% all.eff
         if (!sum(where)) return (-1)
         if (sum(where) > 1) stop("Error, contact the author.")
         index <- (1:length(names.random))[where]
         if (!randomInt) index <- index - 1
         return(index)
       }
       if (nX) indb <- as.numeric(apply(matrix(cnamesX, ncol = 1), 1, find.indeces))
     }       
   }
   else{
     names.random <- character(0)
   }     
   if (randomInt) names.random <- c("(Intercept)", names.random)
   nfixed <- nX - (nrandom - 1*randomInt)
   
   ## Give indeces of factors in the design matrix, it was used to define MH blocks in the earlier version of this program
   n.factors <- 0
   n.in.factors <- NULL
   factors <- NULL
   if (nX){
     temp <- attr(X, "assign")
     if (length(temp) == 1) factors <- 0
     else{
       factors <- numeric(length(temp))
       n.in.factors <- numeric(0)
       temp  <- temp - c(0, temp[1:(length(temp)-1)])
       i <- length(temp)
       while (i >= 1){
         if (temp[i] == 0){
           n.factors <- n.factors + 1
           factors[i] <- n.factors
           n.in.factor <- 1
           while (temp[i-1] == 0){
             i <- i - 1
             factors[i] <- n.factors
             n.in.factor <- n.in.factor + 1
           }
           i <- i - 1
           factors[i] <- n.factors
           n.in.factors <- c(n.in.factor + 1, n.in.factors)
         }
         else
           n.in.factors <- c(1, n.in.factors)         
         i <- i - 1
       }         
     }
     if (length(temp) != nX) stop("Something is wrong, contact the author.")
   }

   ## Sort everything according to the cluster indicator
   ##   and find the numbers of observations per cluster
   if (length(cluster)){
     ordering <- order(cluster)
     Y <- Y[ordering, ]
     cluster <- cluster[ordering]
     rnamesX <- rnamesX[ordering]
     if (nX){
       namesX <- cnamesX
       if (nX == 1) X <- matrix(X[ordering], ncol = 1)
       else         X <- as.matrix(X[ordering, ])
       colnames(X) <- namesX
     }
     ncluster <- length(attr(cluster, "levels"))
     helpf <- function(cl){return(sum(cluster %in% attr(cluster, "levels")[cl]))}
     nwithin <- apply(matrix(1:ncluster, ncol = 1), 1, "helpf")
   }
   else{
     if (nX){
       namesX <- cnamesX
       if (nX == 1) X <- matrix(X, ncol=1)
       else         X <-  as.matrix(X[,])
       colnames(X) <- namesX
     }
     cluster <- 1:n
     ncluster <- n
     nwithin <- rep(1, n)
   }     
   
   ## Transform the response
   if (type == 'interval') {
      if (any(Y[,3]==3)) Y <- cbind(eval(call("transform", Y[,1:2])), Y[,3])
      else               Y <- cbind(eval(call("transform", Y[,1])), Y[,3])
   }
   else if (type == 'left'){
          Y <- cbind(eval(call("transform", Y[,1])), 2-Y[,2])   ## change 0 indicator into 2 indicating left censoring
        }
        else  ## type = 'right' or 'interval2'
           Y <- cbind(eval(call("transform", Y[,1])), Y[,2])

   if (!all(is.finite(Y))) stop("Invalid survival times for this distribution (infinity on log-scale not allowed). ")

   design <- list(n = n, ncluster = ncluster, nwithin = nwithin, nY = nY, nX = nX,
                  nfixed = nfixed, nrandom = nrandom, randomInt = randomInt,
                  Y = Y, X = X, Yinit = Yinit, Xinit = Xinit, cluster = cluster, indb = indb,
                  rnamesX = rnamesX, names.random = names.random, factors = factors, n.factors = n.factors, n.in.factors = n.in.factors)

   return(design)
}  
