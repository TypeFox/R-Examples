predict.plsRmodel <- function(object,newdata,comps=object$computed_nt,type=c("response","scores"),weights,methodNA="adaptative",verbose=TRUE,...)
{
    if (!inherits(object, "plsRmodel")) 
        stop("Primary argument much be a plsRmodel object")
    if (!(type %in% c("response","scores"))) 
      stop("Invalid type specification")
    if (comps>object$computed_nt){stop("Cannot predict using more components than extracted.")}
    type <- match.arg(type)
    if (missing(newdata) || is.null(newdata)) {
    nrtt <- nrow(object$tt)
    if (type=="response"){return(attr(object$RepY,"scaled:center")+attr(object$RepY,"scaled:scale")*object$tt%*%object$CoeffC)}
    if (type=="scores"){return(object$tt[,1:comps])}   
    } else {
    nrnd <- nrow(newdata)
    if(any(apply(is.na(newdata),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of newdata is completely filled with missing data\n"); stop()}
    if(any(is.na(newdata))) {na.miss.newdata <- TRUE} else {na.miss.newdata <- FALSE}
    if(!is.null(object$call$formula)){
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("subset", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$data <- newdata
    mf$formula <- object$call$formula[-2]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- na.pass    
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
#    attr(mt,"intercept")<-0L    
    newdata.frame <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)[,-1]
    else matrix(, nrnd, 0L)
    weights <- as.vector(model.weights(mf))
    } else {newdata.frame <- newdata}
    newdata.scaled <- sweep(sweep(newdata.frame, 2, attr(object$ExpliX,"scaled:center")), 2 ,attr(object$ExpliX,"scaled:scale"), "/")
    newdataNA <- !is.na(newdata)
    newdata.scaledwotNA <- as.matrix(newdata.scaled)
    newdata.scaledwotNA[!newdataNA] <- 0
    tt <- NULL
    if (methodNA=="adaptative") {
    for(ii in 1:nrnd){
    if (all(newdataNA[ii,])){
    tt <- rbind(tt, c(newdata.scaledwotNA[ii,]%*%object$wwetoile[,1:comps],rep(0,object$computed_nt-comps))) 
    }
    else {
      if(verbose){cat("Missing value in row ",ii,".\n")}
      tt <- rbind(tt, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
    }}}
    if(methodNA=="missingdata") {
      if(verbose){cat("Prediction as if missing values in every row.\n")}
    for (ii in 1:nrnd) {  
      tt <- rbind(tt, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
    }
    }
    colnames(tt) <- paste("Comp_",1:object$computed_nt,sep="")
    if (type=="response"){return(attr(object$RepY,"scaled:center")+attr(object$RepY,"scaled:scale")*tt%*%object$CoeffC)}
    if (type=="scores"){return(tt[,1:comps,drop=FALSE])}      
}
}
