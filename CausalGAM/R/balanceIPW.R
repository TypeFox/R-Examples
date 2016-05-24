

balance.IPW <- function(pscore.formula, pscore.family,
                        treatment.var,
                        outcome.var,
                        data=NULL,
                        divby0.action=c("fail", "truncate", "discard"),
                        divby0.tol=1e-8,
                        nboot=501,
                        suppress.warnings=TRUE, ...){
  
  
  if (is.null(data)){
    stop("'data' must be specified in call to 'balance.IPW'\n")
  }

  if (nboot < 0){
    nboot <- 0
  }
  
  call <- match.call()
  arg.list <- as.list(call)
  
  
  divby0.action <- match.arg(divby0.action)

  if (max(is.na(data)) > 0){
    stop("'data' contains NAs. Remove NAs and call balance.IPW again.\n")
  }


  treatment.vec <- data[, treatment.var]
  
  if (is.factor(treatment.vec)){
    treatment.values <- levels(treatment.vec)
  }
  else{
    treatment.values <- sort(unique(treatment.vec))
  }


  non.numeric.indic <- rep(FALSE, ncol(data))
  for (i in 1:ncol(data)){
    if (!is.numeric(data[,i])){
      data[,i] <- -999
      non.numeric.indic[i] <- TRUE
    }
  }



  
  treated.data <- data[treatment.vec==treatment.values[2],]
  control.data <- data[treatment.vec==treatment.values[1],]

  
  if (suppress.warnings){
    gam.ps <- suppressWarnings(gam(pscore.formula, family=pscore.family,
                  data=data, na.action="na.fail", ...))
  }
  else{
    gam.ps <- gam(pscore.formula, family=pscore.family,
                  data=data, na.action="na.fail", ...)
  }


  
  pscore.probs <- predict(gam.ps, newdata=data, type="response", ...)


  

  
  pscores.pre <- pscore.probs
    
  truncated.indic <- rep(FALSE, nrow(data))
  discarded.indic <- rep(FALSE, nrow(data))
  if(min(pscore.probs) <= divby0.tol || max(pscore.probs) >= (1-divby0.tol)){
    if (divby0.action == "fail"){
      stop("\nCannot compute IPW estimate because some\nprobabilities of treatment are numerically 0 and/or 1\n\n")
    }
    if (divby0.action == "truncate"){
      truncated.indic[pscore.probs <= divby0.tol] <- TRUE
      pscore.probs[pscore.probs <= divby0.tol] <- divby0.tol
      truncated.indic[pscore.probs >= (1-divby0.tol)] <- TRUE
      pscore.probs[pscore.probs >= (1-divby0.tol)] <- (1-divby0.tol)
    }
    if (divby0.action == "discard"){
      discarded.indic[pscore.probs <= divby0.tol] <- TRUE
      discarded.indic[pscore.probs >= (1-divby0.tol)] <- TRUE
      pscore.probs <- pscore.probs[!discarded.indic]
      treatment.vec <- treatment.vec[!discarded.indic]
      data <- data[!discarded.indic,]
      treated.data <- data[treatment.vec==treatment.values[2],]
      control.data <- data[treatment.vec==treatment.values[1],]
    }    
  }





  ## number of observations
  n <- length(treatment.vec)

  control.indic <- treatment.vec == treatment.values[1]
  treated.indic <- treatment.vec == treatment.values[2]


  norm1 <- 1 / sum(treated.indic/pscore.probs)
  norm0 <- 1 / sum((1-treated.indic)/(1-pscore.probs))

  data1 <- data
  data1[treatment.vec==treatment.values[1],] <- 0
  
  data0 <- data
  data0[treatment.vec==treatment.values[2],] <- 0
  
  for (i in 1:ncol(data1)){
    data1[,i] <- data1[,i] / pscore.probs
    data0[,i] <- data0[,i] / (1-pscore.probs)
  }

  
  drop.indic <- (colnames(data) == treatment.var) | (colnames(data) == outcome.var)
  data1 <- data1[,!drop.indic]
  data0 <- data0[,!drop.indic]
  treated.data <- treated.data[,!drop.indic]
  control.data <- control.data[,!drop.indic]
  non.numeric.indic <- non.numeric.indic[!drop.indic]
  
  w.mean.1 <- norm1 * apply(data1, 2, sum)
  w.mean.0 <- norm0 * apply(data0, 2, sum)
  obs.mean.1 <- colMeans(treated.data)
  obs.mean.0 <- colMeans(control.data)
  for (i in 1:ncol(data1)){
    if (non.numeric.indic[i]){
      w.mean.1[i] <- NA
      w.mean.0[i] <- NA
      obs.mean.1[i] <- NA
      obs.mean.0[i] <- NA
    }
  }



  


  ## calculate bootstrap SEs
  w.diff.SE <- NULL
  if (nboot > 0){

    bs.mat <- matrix(NA, nboot, ncol(data1))
    for (biter in 1:nboot){
      boot.inds <- sample(1:n, n, replace=TRUE)
      data.bs <- data[boot.inds,]

      bs.out <- balance.IPW(pscore.formula=pscore.formula,
                            pscore.family=pscore.family,
                            treatment.var=treatment.var,
                            outcome.var=outcome.var,
                            data=data.bs,
                            divby0.action=divby0.action,
                            divby0.tol=divby0.tol,
                            nboot=0,
                            suppress.warnings=TRUE, ...)



      bs.mat[biter,] <- bs.out$weighted.mean.treated - bs.out$weighted.mean.control
    }

    w.diff.SE <- apply(bs.mat, 2, sd)
    for (i in 1:ncol(data1)){
      if (non.numeric.indic[i]){
        w.diff.SE[i] <- NA
      }
    }
    
  } ## end nboot > 0
  

  return(structure(list(obs.mean.control=obs.mean.0,
                        obs.mean.treated=obs.mean.1,
                        weighted.mean.control=w.mean.0,
                        weighted.mean.treated=w.mean.1,
                        weighted.diff.SE=w.diff.SE),
         class="balance"))
         
}





"print.balance" <- function(x, ...){

  w.diff <- x$weighted.mean.treated - x$weighted.mean.control

  mat <- cbind(x$obs.mean.treated,
               x$obs.mean.control,
               x$weighted.mean.treated,
               x$weighted.mean.control,
               w.diff,
               w.diff/x$weighted.diff.SE)
  colnames(mat) <- c("obs.mean.t", "obs.mean.c",
                     "w.mean.t", "w.mean.c", "w.mean.diff", "z")
  mat <- round(mat, 3)
  print(mat)

}











