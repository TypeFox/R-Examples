predict.plsRcoxmodel <- function(object,newdata,comps=object$computed_nt,type=c("lp", "risk", "expected", "terms", "scores"),se.fit=FALSE,weights,methodNA="adaptative", verbose=TRUE,...)
{
  if (!inherits(object, "plsRcoxmodel")) 
    stop("Primary argument much be a plsRcoxmodel object")
  if(missing(type)){type="lp"}
  if (!(type %in% c("lp", "risk", "expected", "terms", "scores")))
    stop("Invalid type specification")
  if (comps>object$computed_nt){stop("Cannot predict using more components than extracted.")}
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    nrtt <- nrow(object$tt)
    if (type=="lp"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "lp",se.fit=se.fit))
    }
    if (type=="risk"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "risk",se.fit=se.fit))
    }
    if (type=="expected"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "expected",se.fit=se.fit))
    }
    if (type=="terms"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
    ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "terms",se.fit=se.fit))
    }  
    if (type=="scores"){return(object$tt[,1:comps])}   
  } else {
    nrnd <- nrow(newdata)
    if(any(apply(is.na(newdata),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of newdata is completely filled with missing data\n"); stop()}
    if(any(is.na(newdata))) {na.miss.newdata <- TRUE} else {na.miss.newdata <- FALSE}
    if(!is.null(object$XplanFormula)){
    mf0 <- match.call(expand.dots = FALSE)
    m0 <- match(c("weights"), names(mf0), 0L)
    mf0 <- mf0[c(1L, m0)]
    mf0$data <- newdata
    mf0$formula <- object$XplanFormula
    mf0$drop.unused.levels <- TRUE
    mf0[[1L]] <- as.name("model.frame")
    mf0 <- eval(mf0, parent.frame())
    mt0 <- attr(mf0, "terms")
    attr(mt0,"intercept")<-0L
    newdata.frame <- if (!is.empty.model(mt0)) model.matrix(mt0, mf0, contrasts)[,-1]
    else matrix(, nrnd, 0L)
    weights <- as.vector(model.weights(mf0))
    } else {newdata.frame <- newdata}
    newdata.scaled <- sweep(sweep(newdata.frame, 2, attr(object$ExpliX,"scaled:center")), 2 ,attr(object$ExpliX,"scaled:scale"), "/")
  newdataNA <- !is.na(newdata)
  newdata.scaledwotNA <- as.matrix(newdata.scaled)
  newdata.scaledwotNA[!newdataNA] <- 0
  ttpredY <- NULL
  if (methodNA=="adaptative") {
    for(ii in 1:nrnd){
      if (all(newdataNA[ii,])){
        ttpredY <- rbind(ttpredY, c(newdata.scaledwotNA[ii,]%*%object$wwetoile[,1:comps],rep(0,object$computed_nt-comps))) 
      }
      else {
        if(verbose){cat("Missing value in row ",ii,".\n")}
        ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
      }}}
  if(methodNA=="missingdata") {
    if(verbose){cat("Prediction as if missing values in every row.\n")}
    for (ii in 1:nrnd) {  
      ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
    }
  }
colnames(ttpredY) <- NULL
ttpred<-data.frame(tt=ttpredY)
    if (type=="lp"){return(predict(object$FinalModel,newdata=ttpred,type = "lp",se.fit=se.fit))
    }
    if (type=="risk"){return(predict(object$FinalModel,newdata=ttpred,type = "risk",se.fit=se.fit))
    }
    if (type=="expected"){return(predict(object$FinalModel,newdata=ttpred,type = "expected",se.fit=se.fit))
    }
    if (type=="terms"){return(predict(object$FinalModel,newdata=ttpred,type = "terms",se.fit=se.fit))
    }  

    if (type=="scores"){colnames(ttpredY) <- paste("Comp_",1:object$computed_nt,sep="");
                      return(ttpredY[,1:comps,drop=FALSE])}
}
}

