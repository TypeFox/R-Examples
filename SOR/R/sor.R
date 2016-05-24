sor <- function(y.formula,
                w1.formula,
                w2.formula = ~1,
                id,
                family = "binomial",
                y0 = 0,
                hfunc = identity, 
                support = c(0,1),
                pi1.pi0.ratio = 1,  
                data = parent.frame(),
                init.beta=NULL,
                init.sig.2 = 1,
                est.var = TRUE,
                CORSTR="independence"){
  call <- match.call()
  
  fam <- as.character(charmatch(family, c("gaussian", "normal", "poisson", "binomial")))
  if(fam=="2") fam <- "1"
  
  if(typeof(data) == "environment"){
    id = id
    DAT.ods <- data.frame(model.frame(y.formula), model.frame(w1.formula), model.frame(w2.formula))
  }
  else{
    DAT.ods <- data 
    subj.col <- which(colnames(data) == call$id)  
    id <- data[,subj.col]
  }
  ny <- length(support)  
  DAT.ods$id       <- id
  DAT.ods$offset.z <- log(pi1.pi0.ratio) 
  DAT.ods$pi.ratio <- pi1.pi0.ratio
  
  yframe <- model.frame(y.formula, DAT.ods)
  resp <- model.response(yframe)
  w1frame <- model.frame(w1.formula, DAT.ods)
  ref <- model.response(w1frame)
  w2frame <- model.frame(w2.formula, DAT.ods)
  
  if(!all( is.element(ref, c(0,1)))){stop("Response in w1.formula must be binary")}
  if(!all(ref==model.response(w2frame)) && !is.null(model.response(w2frame))){stop("There should be no response for w2.formula")}
  
  if(max(resp) > max(support) | min(resp) < min(support)) warning("Some values of response are outside support")
  
  switch(fam,
          #### NORMAL ####
          "1"= {
            results <- norm.SOR(y.formula, w1.formula, w2.formula, y0, hfunc, id, support, pi1.pi0.ratio, DAT.ods, init.beta, init.sig.2, est.var, CORSTR)                
            },
          #### POISSONS ####
          "3" = {
            if(!all( is.wholenumber(resp))){stop("Response in y.formula must be positive integers if family is Poisson")}
            if(init.sig.2 != 1){warning("Initial variance specification ignored")}
            results <- pois.SOR(y.formula, w1.formula, w2.formula, y0, hfunc, id, support, pi1.pi0.ratio, DAT.ods, init.beta,  CORSTR)
           },
          #### BINOMIAL ####
           "4" = {
             if(!all( is.element(resp, c(0,1)))){stop("Response in y.formula must be binary if family is binomial")}
             if(y0 != 0){warning("Setting y0 to 0")}
             results <- binom.SOR(y.formula, w1.formula, w2.formula, DAT.ods, init.beta,  CORSTR)
           }
         )
  class(results) <- "sor"
  results$coefnames <- colnames(model.matrix(y.formula, DAT.ods))
  results$call <- call
  results$family <- c("normal","normal", "Poisson", "binomial")[as.numeric(fam)]
  return(results)
}

print.sor <- function(x, ...){
  coefdf <- data.frame(x$coefs)
  rownames(coefdf) <- x$coefnames
  colnames(coefdf) <- ""
  print(x$call)
  cat("\n", "Coefficients:", "\n")
  print(t(coefdf))
  cat("\n Reference Distribution: ", x$family, "\n")
}

summary.sor <- function(object, ...)  {
  Coefs <- matrix(0,nrow=length(object$coefs),ncol=4)
  Coefs[,1] <- c(object$coefs)
  Coefs[,2] <- object$se.coefs
  Coefs[,3] <- Coefs[,1]/Coefs[,2]
  Coefs[,4] <- round(2*pnorm(abs(Coefs[,3]), lower.tail=F), digits=8)
  colnames(Coefs) <- c("Estimates","Std. Err.", "wald", "p")
  
  summ <- list(coefs = Coefs[,1], se.coefs = Coefs[,2],  wald.test = Coefs[,3], p = Coefs[,4],
                coefnames = object$coefnames, family=object$family, call=object$call)
  
  class(summ) <- 'summary.sor'
  return(summ)
}

print.summary.sor <- function(x, ...){
  Coefs <- matrix(0,nrow=length(x$coefnames),ncol=4)
  rownames(Coefs) <- c(x$coefnames)
  colnames(Coefs) <- c("Estimates","Std. Err.", "wald", "p")
  Coefs[,1] <- x$coefs
  Coefs[,2] <- x$se.coefs
  Coefs[,3] <- x$wald.test
  Coefs[,4] <- x$p
  
  print( x$call)
  cat("\n")
  print(Coefs)
  cat("\n Reference Distribution: ", x$family, "\n")
}

