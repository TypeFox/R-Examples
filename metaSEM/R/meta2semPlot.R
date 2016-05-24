meta2semPlot <- function(object, manNames=NULL, latNames=NULL, labels=c("labels", "RAM"), ...) {

    if (!requireNamespace("semPlot", quietly=TRUE))    
        stop("\"semPlot\" package is required for this function.")
    
    if (class(object)=="wls") {
    ## No visible binding for global variable Note in R CMD check    
    ## A <- mxEval(Amatrix, object$mx.fit)
    ## S <- mxEval(Smatrix, object$mx.fit)
    ## F <- mxEval(Fmatrix, object$mx.fit)
       A <- object$mx.fit@matrices$Amatrix$values
       S <- object$mx.fit@matrices$Smatrix$values
       ## S may be calculated when diag.constraints=FALSE
       if (is.null(S)) S <- object$mx.fit$algebras$Smatrix$result      
       F <- object$mx.fit@matrices$Fmatrix$values
       Id <- diag(nrow(S))
       ObsCovs <- object$Cov
    ## ImpCovs <- mxEval(Fmatrix%*%solve(Id-Amatrix)%*%S%*%solve(t(Id-Amatrix))%*%t(Fmatrix), 
    ##                   object$mx.fit)
    ## ImpCovs <- F%*%solve(Id-A)%*%S%*%solve(t(Id-A))%*%t(F)
       ImpCovs <- object$mx.fit@algebras$impliedS$result    
        
       if (is.null(manNames)) {
               manNames <- colnames(ObsCovs)
       }    
       ## If colnames(ObsCovs) is still NULL
       if (is.null(manNames)) {
           manNames <- paste("X", seq(1, nrow(F), 1), sep="")
       }

       dimnames(ImpCovs) <- dimnames(ObsCovs) <- list(manNames, manNames)
            
       if (is.null(latNames)) {
         no.latent <- ncol(F) - nrow(F)
         if (no.latent > 0) latNames <- paste("L", seq(1,no.latent), sep="")
       }
    
       out <- semPlot::ramModel(A=A, S=S, F=F, manNames=manNames, latNames=latNames, ObsCovs=ObsCovs, 
                                ImpCovs=ImpCovs, ...)
    
       labels <- match.arg(labels)    
       if (labels=="labels") {
         A.labels <- object$mx.fit@matrices$Amatrix$labels
         S.labels <- object$mx.fit@matrices$Smatrix$labels
         labels <- c(c(A.labels), c(S.labels))
         row.pars <- as.numeric(row.names(out@Pars))
         out@Pars$label <- labels[row.pars]
         ## replace NA labels with estimates
         na.index <- is.na(out@Pars$label)
         out@Pars$label[na.index] <- sprintf("%.2f", out@Pars$est[na.index])
       }   
  } else {
    stop("\"object\" must be an object of class \"wls\".")
    
  }
    
  out
}
