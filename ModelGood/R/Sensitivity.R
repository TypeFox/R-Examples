##' Compute sensitivity, specificity and predictive values
##'
##' Confidence intervals are obtained with \code{binom.test} 
##' @aliases Specificity NPV PPV Diagnose
##' @usage Sensitivity(x,event,cutoff,comparison=">=",...)
##'        Specificity(x,event,cutoff,comparison=">=",...)
##'        NPV(x,event,cutoff,comparison=">=",...)
##'        PPV(x,event,cutoff,comparison=">=",...)
##' @title Compute sensitivity, specificity and predictive values
##' @param x Either a binary 0,1 variable, or a numeric marker which is cut into binary. 
##' @param event Binary response variable. Either a 0,1 variable where 1 means 'event', or a factor where the second level means 'event'.
##' @param cutoff When x is a numeric marker, it is compared to this cutoff to obtain a binary test.
##' @param comparison How x is to be compared to the cutoff value
##' @param ... passed on to \code{binom.test}
##' @return list with Sensitivity, Specificity, NPV, PPV and confidence interval
##' @seealso binom.test
##' @examples
##'
##' set.seed(17)
##' x <- rnorm(10)
##' y <- rbinom(10,1,0.4)
##' Sensitivity(x,y,0.3)
##' Specificity(x,y,0.3)
##' PPV(x,y,0.3)
##' NPV(x,y,0.3)
##' 
##' Diagnose(x,y,0.3)
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
Sensitivity <- function(x,event,cutoff,comparison=">=",...){
stopifnot(length(x)==length(event))
    if (length(unique(x))<=2){
        if (is.factor(x))
            Test <- x==levels(x)[1]
        else
            Test <- x
    } else{
        if (missing(cutoff)) stop("No cutoff given")
        Test <- switch(comparison,"<="=I(x<=cutoff),"<"=I(x<cutoff),">="=I(x>=cutoff),">"=I(x>cutoff))
        Test <- 1*Test
    }
    ue <- unique(event)
    stopifnot(length(ue)==2)
    if (is.factor(event))
        Event <- 1*(event==levels(event)[2])
    else
        Event <- 1*(event==max(event))
    ## print(table(Test,Event))
    bt <- binom.test(x=sum(Test==1 & Event==1),n=sum(Event==1),...)
    out <- list(Sensitivity=as.numeric(bt$estimate),conf.int=bt$conf.int)
    class(out) <- "Sensitivity"
    out
}

##' @S3method print Sensitivity
print.Sensitivity <- function(x,...){
    cat(names(x)[1],
        ": ",
        round(100*x[[1]],1),
        " (CI_95:",
        paste("[",paste(round(100*x$conf.int,1),collapse=","),"])",sep=""),
        "\n",
        sep="")
}

##' @export
Specificity <- function(x,event,cutoff,comparison=">=",...){
    stopifnot(length(x)==length(event))
    if (length(unique(x))<=2){
        if (is.factor(x))
            Test <- x==levels(x)[1]
        else
            Test <- x
    } else{
        if (missing(cutoff)) stop("No cutoff given")
        Test <- switch(comparison,"<="=I(x<=cutoff),"<"=I(x<cutoff),">="=I(x>=cutoff),">"=I(x>cutoff))
        Test <- 1*Test
    }
    ue <- unique(event)
    stopifnot(length(ue)==2)
    if (is.factor(event))
        Event <- 1*(event==levels(event)[2])
    else
        Event <- 1*(event==max(event))
    ## print(table(Test,Event))
    bt <- binom.test(x=sum(Test==0 & Event==0),n=sum(Event==0),...)
    out <- list(Specificity=as.numeric(bt$estimate),conf.int=bt$conf.int)
    class(out) <- "Sensitivity"
    out
}

##' @export
PPV <- function(x,event,cutoff,comparison=">=",...){
    stopifnot(length(x)==length(event))
    if (length(unique(x))<=2){
        if (is.factor(x))
            Test <- x==levels(x)[1]
        else
            Test <- x
    } else{
        if (missing(cutoff)) stop("No cutoff given")
        Test <- switch(comparison,"<="=I(x<=cutoff),"<"=I(x<cutoff),">="=I(x>=cutoff),">"=I(x>cutoff))
        Test <- 1*Test
    }
    ue <- unique(event)
    stopifnot(length(ue)==2)
    if (is.factor(event))
        Event <- 1*(event==levels(event)[2])
    else
        Event <- 1*(event==max(event))
    ## print(table(Test,Event))
    bt <- binom.test(x=sum(Test==1 & Event==1),n=sum(Test==1),...)
    out <- list("Positive predictive value"=as.numeric(bt$estimate),conf.int=bt$conf.int)
    class(out) <- "Sensitivity"
    out
}
##' @export
NPV <- function(x,event,cutoff,comparison=">=",...){
    stopifnot(length(x)==length(event))
    if (length(unique(x))<=2){
        if (is.factor(x))
            Test <- x==levels(x)[1]
        else
            Test <- x
    } else{
        if (missing(cutoff)) stop("No cutoff given")
        Test <- switch(comparison,"<="=I(x<=cutoff),"<"=I(x<cutoff),">="=I(x>=cutoff),">"=I(x>cutoff))
        Test <- 1*Test
    }
    ue <- unique(event)
    stopifnot(length(ue)==2)
    if (is.factor(event))
        Event <- 1*(event==levels(event)[2])
    else
        Event <- 1*(event==max(event))
    ## print(table(Test,Event))
    bt <- binom.test(x=sum(Test==0 & Event==0),n=sum(Test==0),...)
    out <- list("Negative predictive value"=as.numeric(bt$estimate),conf.int=bt$conf.int)
    class(out) <- "Sensitivity"
    out
}

##' @export
Diagnose <- function(x,event,cutoff,comparison=">=",...){
    stopifnot(length(x)==length(event))
    if (length(unique(x))<=2){
        if (is.factor(x))
            Test <- x==levels(x)[1]
        else
            Test <- x
    } else{
        if (missing(cutoff)) stop("No cutoff given")
        Test <- switch(comparison,"<="=I(x<=cutoff),"<"=I(x<cutoff),">="=I(x>=cutoff),">"=I(x>cutoff))
        Test <- 1*Test
    }
    ue <- unique(event)
    stopifnot(length(ue)==2)
    if (is.factor(event))
        Event <- 1*(event==levels(event)[2])
    else
        Event <- 1*(event==max(event))
    sens <- binom.test(x=sum(Test==1 & Event==1),n=sum(Event==1),...)
    spec <- binom.test(x=sum(Test==0 & Event==0),n=sum(Event==0),...)
    ppv <- binom.test(x=sum(Test==1 & Event==1),n=sum(Test==1),...)
    npv <- binom.test(x=sum(Test==0 & Event==0),n=sum(Test==0),...)
    out <- list(table=table(test=Test, event=Event),
                diag=list(Sensitivity=list(Sensitivity=sens$estimate,conf.int=sens$conf.int),
                    Specificity=list(Specificity=spec$estimate,conf.int=spec$conf.int),
                    PPV=list(PPV=ppv$estimate,conf.int=ppv$conf.int),
                    NPV=list(NPV=npv$estimate,conf.int=npv$conf.int)))
    class(out) <- "Diagnose"
    out
}


##' @S3method print Diagnose
print.Diagnose <- function(x,...){
    cat("2x2 table:\n")
    print(x$table)
    dd <- do.call("rbind",lapply(x$diag,function(d){
        c(names(d)[1],
          round(100*d[[1]],1),
          paste("[",paste(round(100*d$conf.int,1),collapse=","),"]",sep=""))
    }))
    colnames(dd) <- c("Parameter","Estimate","CI.95")
    rownames(dd) <- NULL
    cat("\nDiagnostic parameters:\n")
    print(dd,quote=FALSE,rownames=FALSE)
    invisible(dd)
}
