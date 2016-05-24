########################################################################
# COR MAIN FUNCTION
########################################################################
svyPVcor <- function(formula, 
                     design, 
                     placeholder = 1:10){
  # prepare args
  x.name <- as.character(formula)[2]
  y.name <- as.character(formula)[3] 
  if(length(grep("\\.\\.", x.name))==1){
    x.name  <- lapply(placeholder,
                      function(x) gsub("\\.\\.", x, x.name))
  }
  if(length(grep("\\.\\.", y.name))==1){
    y.name  <- lapply(placeholder,
                      function(x) gsub("\\.\\.", x, y.name))
  }
  x.name.4.df <- unlist(x.name)
  y.name.4.df <- unlist(y.name)
  Y <- as.data.frame(design$variables[y.name.4.df])
  X <- as.data.frame(design$variables[x.name.4.df])
  # create NA-dummy to count cases/weights & for survey update
  NA.dummy.matrix <- cbind(X,Y)
  NA.dummy <- apply(NA.dummy.matrix, 
                    1, 
                    function(x){
                      any(is.na(x))}
  )
  n <- sum(design$pweights[NA.dummy==FALSE])
  N <- sum(NA.dummy==FALSE)
  updated.design <- update(design, NA.dummy = NA.dummy)
  # computation with two vectors as input
  if(max(dim(as.matrix(X))[2],
         dim(as.matrix(Y))[2]) == 1){
    result <- cor.aux.func(X[NA.dummy == FALSE,], 
                           Y[NA.dummy == FALSE,], 
                           subset(updated.design, 
                                  NA.dummy == FALSE)
    )  
    finalres <- rbind(t(as.data.frame(result)), N, n) 
    rownames(finalres) <- c("COR", "SE", "number.of.cases",
                            "sum.of.weights")
    t(finalres)
  }
  else{
    # computation with one vectors and one matrix as input
    if(dim(as.matrix(X))[2] == 1 &
         dim(as.matrix(Y))[2] > 1){
      X.d <- matrix(data = unlist(rep(X, dim(as.matrix(Y))[2])), 
                    ncol = dim(as.matrix(Y))[2])  
      PVs.d <- lapply(1:dim(as.matrix(X.d))[2], 
                      function(n){
                        inp1 <- X.d[NA.dummy == FALSE,n]
                        inp2 <- Y[NA.dummy == FALSE,n]
                        cor.aux.func(inp1, 
                                     inp2, 
                                     subset(updated.design, 
                                            NA.dummy == FALSE))}
      )
    }  
    if(dim(as.matrix(X))[2] > 1 &
         dim(as.matrix(Y))[2] == 1){
      Y.d <- matrix(data = unlist(rep(Y, dim(as.matrix(X))[2])), 
                    ncol = dim(as.matrix(X))[2])  
      PVs.d <- lapply(1:dim(as.matrix(Y.d))[2], 
                      function(n){
                        inp1 <- X[NA.dummy == FALSE,n]
                        inp2 <- Y.d[NA.dummy == FALSE,n]
                        cor.aux.func(inp1, 
                                     inp2, 
                                     subset(updated.design, 
                                            NA.dummy == FALSE))}
      )
    }  
    # computation with two matrices as input
    if(dim(as.matrix(X))[2] > 1 &
         dim(as.matrix(Y))[2] > 1){
      PVs.d <- lapply(1:dim(as.matrix(X))[2], function(n){
        inp1 <- X[NA.dummy == FALSE,n]
        inp2 <- Y[NA.dummy == FALSE,n]
        cor.aux.func(inp1, inp2, subset(updated.design, 
                                        NA.dummy == FALSE))}
      ) 
    }
    # dataframe with coefficients
    unlPV.with.lab <- as.data.frame(PVs.d)
    # point estimators + SD
    ndummy1 <- seq(1,2*max(dim(as.matrix(X))[2], 
                           dim(as.matrix(Y))[2]),2)
    ndummy2 <- 1+seq(1,2*max(dim(as.matrix(X))[2],
                             dim(as.matrix(Y))[2]),2)
    PV.means <- unlist(unlPV.with.lab[ndummy1])
    PV.ses.d <- unlist(unlPV.with.lab[ndummy2])
    PV.ses <- PV.ses.d^2
    Means <- mean(PV.means)
    sampling.v <- mean(PV.ses)
    imputation.v <- var(PV.means) 
    tot.v <- sqrt(sampling.v + ((1+1/length(placeholder)) * 
                                  imputation.v))
    # compilation of final result object
    result <- rbind(Means,tot.v,N,n)
    rownames(result) <- c("COR", "SE", "number.of.cases",
                          "sum.of.weights")
    # return final result object
    t(result)
  }
}
