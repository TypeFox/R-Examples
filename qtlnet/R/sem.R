#################################################################################################
qdg.sem <- function(qdgObject, cross)
{
#################################################################################
  score.sem.models <- function(cross,pheno.names,all.solutions,steptol,addcov=NULL) {
    n.sol <- length(all.solutions[[1]])
    mypheno <- cross$pheno[,pheno.names]
    np <- length(mypheno[1,])
    n.paths <- nrow(all.solutions[[1]][[1]])
    semBIC <- rep(NA,n.sol)
    path.coeffs <- matrix(NA,n.paths,n.sol)
    if(!is.null(addcov)){
      addcov <- paste("cross$pheno$",addcov,sep="")
      myresid <- matrix(0,nind(cross),np)
      for(i in 1:np){
        fm <- lm(as.formula(paste("mypheno[,i] ~ ", paste(addcov, collapse = "+"))))
        myresid[,i] <- fm$resid
      }
      mycov <- cov(myresid)
      for(i in 1:n.sol){
        ramMatrix <- create.sem.model(DG=all.solutions[[1]][[i]],pheno.names=pheno.names)
        mysem <- try(sem(ramMatrix, S = mycov, N = nind(cross), var.names = pheno.names,
                         steptol = steptol, analytic.gradient = FALSE,
                         param.names = paste("Param", seq(nrow(ramMatrix)), sep = "")), silent = TRUE)
        if(class(mysem)[1] != "try-error"){
          aux.summary <- try(summary(mysem),silent=TRUE)
          if(class(aux.summary)[1] != "try-error"){ 
            semBIC[i] <- aux.summary$BIC
            path.coeffs[,i] <- include.path.coefficients(sem.summary=aux.summary,output=all.solutions[[1]][[i]])
          }
        }
      }
    }
    else {
      mycov <- cov(mypheno)
      for(i in 1:n.sol){
        ramMatrix <- create.sem.model(DG=all.solutions[[1]][[i]],pheno.names=pheno.names)	
        mysem <- try(sem(ramMatrix, S = mycov, N = nind(cross), var.names = pheno.names,
                         steptol = steptol, analytic.gradient = FALSE,
                         param.names = paste("Param", seq(nrow(ramMatrix)), sep = "")), silent = TRUE)
        if(class(mysem)[1] != "try-error"){
          aux.summary <- try(summary(mysem),silent=TRUE)
          if(class(aux.summary)[1] != "try-error"){ 
            semBIC[i] <- aux.summary$BIC
            path.coeffs[,i] <- include.path.coefficients(sem.summary=aux.summary,output=all.solutions[[1]][[i]])
          } 
        }
      }
    }
    
    ## Drop solutions that did not work with sem().
    tmp <- !is.na(semBIC)
    if(!any(tmp)) {
      stop("No qdg solutions could be fit with sem().")
    }
    if(any(!tmp)) {
      warning(paste(sum(!tmp), "qdg solutions could not be fit with sem() and were dropped."))
      semBIC <- semBIC[tmp]
      path.coeffs <- path.coeffs[, tmp, drop = FALSE]
      n.sol <- sum(tmp)
      dropped <- which(!tmp)
    }
    else
      dropped <- NULL
    
    output <- data.frame(cbind(semBIC,approx.posterior(semBIC)))
    names(output) <- c("sem.BIC","posterior prob")
    row.names(output) <- paste("model.",1:n.sol,sep="")
    ## if there are ties, returns the first.
    best <- which(output[,2] == max(output[,2]))[1]	
    list(output,path.coeffs[,best], dropped)
  }
#########################################################
  include.path.coefficients <- function(sem.summary,output) {
    ne <- length(output[,1])
    mypathcoef <- rep(NA,ne)
    aux <- sem.summary$coeff
    aux <- aux[1:ne,]
    for(i in 1:ne){
      if(output[i,2] == "---->") aux1 <- paste(output[i,3], output[i,1], sep=" <--- ")
      if(output[i,2] == "<----") aux1 <- paste(output[i,1], output[i,3], sep=" <--- ")
      aux2 <- match(aux1,aux[,5])
      mypathcoef[i] <- aux[aux2,1]
    }
    mypathcoef
  }
############################################
  create.sem.model <- function(DG,pheno.names) {
    n <- length(DG[,1])
    myvector <- c()
    for(i in 1:n){
      aux1 <- which(DG[i,1]==pheno.names)
      aux2 <- which(DG[i,3]==pheno.names)
      if(DG[i,2] == "---->"){
        aux.vector <- c(1,aux2,aux1,i,NA)
      }
      else{aux.vector <- c(1,aux1,aux2,i,NA)}
      myvector <- c(myvector,aux.vector)
    }
    for(i in 1:length(pheno.names)){
      aux.vector <- c(2,i,i,n+i,NA)
      myvector <- c(myvector,aux.vector)
    }
    matrix(myvector,ncol=5,byrow=TRUE)
  }
##################################
  approx.posterior <- function(bics) {
    aux <- min(bics)
    round(exp(-0.5*(bics-aux))/sum(exp(-0.5*(bics-aux))),6)
  }

#################################################
  ss <- score.sem.models(cross = cross,
                         pheno.names = qdgObject$phenotype.names,
                         all.solutions = qdgObject$Solutions,
                         steptol = 1 / 100000,
                         addcov = qdgObject$addcov)
  best <- which(ss[[1]][,1] == min(ss[[1]][,1]))
  mylist <- list(best, ss[[1]], ss[[2]])
  names(mylist) <- c("best.SEM","BIC.SEM","path.coeffs")
  mylist$Solutions <- qdgObject$Solutions
  mylist$marker.names <- qdgObject$marker.names
  mylist$phenotype.names <- qdgObject$phenotype.names
  mylist$dropped <- ss[[3]]
  class(mylist) <- c("qdg.sem", "qdg", "list")
  
  mylist
}

summary.qdg.sem <- function(object, ...)
{
  cat("\nBest SEM solution:\n")
  print(object$Solution$solution[[object$best.SEM]])
  bic.sem <- object$BIC.SEM[object$best.SEM, "sem.BIC"]
  cat("\nBIC:\n")
  print(c(sem = bic.sem))
  cat("\nBest SEM solution is solution number:\n")
  print(object$best.SEM)
  if(!is.null(object$dropped)) {
    cat(length(object$dropped), "qdg.sem solution were dropped; sem() failed for graphs",
        paste(object$dropped, collapse = ","))
  }
  invisible()
}

print.qdg.sem <- function(x, ...) summary(x, ...)
