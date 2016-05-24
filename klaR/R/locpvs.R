pvs <- function(x, ...) 
{
  UseMethod("pvs")
}


pvs.default <- function(x, grouping, prior=NULL, method="lda", vs.method=c("ks.test","stepclass","greedy.wilks"), 
                 niveau=0.05, fold=10, impr=0.1, direct="backward", out=FALSE, ...) { 

    cl <- match.call()
    cl[[1]] <- as.name("pvs")
    vs.method <- match.arg(vs.method)
    classes <- levels(grouping)
    classcombins <- combn(classes, 2)
    if (is.null(prior)) { prior <- table(grouping)/sum(table(grouping)) }
    names(prior) <- classes

    if (dim(x)[2] == 1) { stop("If only one variable is used, there is no need for pairwise variable selection!") }
    if (sum(!sapply(x,is.numeric)) > 0) { stop("Implemented variable selection methods are only for numerical data.") }
    if (!is.null(prior)) { 
    if (round(sum(prior)) != 1) stop("Sum of priors must equal 1.")
    if (length(prior) != length(levels(grouping))) stop("Prior must be set in harmony with grouping.") 
    }
    if (length(classes) <= 2) { stop("pvs with 2 or less classes makes no sense.") }
    if (sum(!(table(grouping) > 2)) > 0) { stop("There should be at least 5 observations per group to perform pairwise variable selection.") }
    if ( ((method == "svm") || (method == "randomForest")) && (vs.method == "stepclass") ) { 
    stop("The chosen classifier method can't be used together with selection method stepclass.") 
    }


    
## for a given pair of classes (classvec) test whether variables differ,  
## estimate distribution of resulting subspace for both classes:

    calc.classcompare <- function(classvec, pval=niveau, ...) {

        class1 <- which(grouping==classvec[1])
        class2 <- which(grouping==classvec[2])     
        apriori <- as.vector(prior[classvec]/sum(prior[classvec]))   
     
        switch(vs.method,
        "ks.test" = {
                suppressWarnings(pvalues <- 
            apply(x, 2, function(vec) ks.test(vec[class1], vec[class2], alternative="two.sided", exact=TRUE)$p.value))
                relevant <- pvalues <= pval
                if (!sum(relevant)) { relevant <- pvalues <= min(pvalues) } ## if all variables are discarded, build model of those with minimal p.values...
            subgrouping <- factor(grouping[c(class1, class2)],labels=classvec) 
            subdata <- as.matrix(x[c(class1, class2), relevant])
            colnames(subdata) <- colnames(x)[relevant]
        if (method=="svm") {    
            names(apriori) <- classvec
            model <- try(do.call(method, list(subdata, subgrouping, type = "C-classification", 
                                              class.weights=apriori, probability = TRUE, ...)))
            } else {
                if (method=="randomForest") {
                    if(is.null(colnames(subdata)))  warning("If using randomForest at least up to version 4.5-25 for prediction colnames are required!")
                    model <- try(do.call(method, list(subdata, subgrouping, classwt=apriori, ...)))
                    } 
                else {
                    model <- try(do.call(method, list(subdata, subgrouping, prior=apriori, ...)))
                    }
                }
            if (inherits(model, "try-error")) { stop("Method ", sQuote(method), " resulted in the error reported above.") }
            return(list(classpair=classvec,  subspace=list(variables=which(relevant), pvalues=pvalues), model=model))     
        },

        "stepclass" = {            
                which.relevant <- stepclass(x=x[c(class1, class2), ], grouping=grouping[c(class1, class2)], method=method, 
                        prior=apriori, fold=fold, improve=impr, direction=direct, output=out, ...)   
                details <- which.relevant$process
                which.relevant <- which.relevant$model$nr
                relevant <- !logical(ncol(x)) ## include all variables (avoids no variable will be selected)
                if (length(which.relevant)) { relevant[-which.relevant] <- FALSE } ## if any variable selected put it into the model

                subgrouping <- factor(grouping[c(class1, class2)],labels=classvec) 
                subdata <- as.matrix(x[c(class1, class2), relevant])
                colnames(subdata) <- colnames(x)[relevant]
                model <- try(do.call(method, list(subdata, subgrouping, prior=apriori, ...)))
                if (inherits(model, "try-error")) { stop("Method ", sQuote(method), " resulted in the error reported above.") }
                return(list(classpair=classvec,  subspace=list(variables=which(relevant), details=details), model=model))     
            },

        "greedy.wilks" = {
                which.relevant <- greedy.wilks(X=x[c(class1, class2), ], grouping=grouping[c(class1, class2)], 
                                       method=method, prior=apriori, niveau=niveau)
                results <- which.relevant$results
                which.relevant <- as.numeric(which.relevant[[1]]$vars) ## get indexes of chosen variables
                relevant <- !logical(ncol(x)) ## include all variables (avoids no variable will be selected)
                if (length(which.relevant)) { relevant[-which.relevant] <- FALSE } # if any variable selected put it into the model
        
                subgrouping <- factor(grouping[c(class1, class2)],labels=classvec) 
                subdata <- as.matrix(x[c(class1, class2), relevant])
                colnames(subdata) <- colnames(x)[relevant]

        if (method=="svm") { 
            names(apriori) <- classvec
            model <- try(do.call(method, list(subdata, subgrouping, type = "C-classification", 
                                              class.weights=apriori, probability = TRUE, ...)))
            } else {
                if (method=="randomForest") {
                    if(is.null(colnames(subdata)))  warning("If using randomForest at least up to version 4.5-25 for prediction colnames are required!")
                    model <- try(do.call(method, list(subdata, subgrouping, classwt=apriori, ...)))
                    }    
                else {
                    model <- try(do.call(method, list(subdata, subgrouping, prior=apriori, ...)))
                    }
                }
            if (inherits(model, "try-error")) { stop("Method ", sQuote(method), " resulted in the error reported above.") }
            return(list(classpair=classvec,  subspace=list(variables=which(relevant), results=results), model=model))     
           } 
        ) ## end switch

    } ## end calc.classcompare
        
    models <- apply(classcombins, 2, calc.classcompare, ...)
    
    result <- list(classes=classes, prior=prior, method=method, vs.method=vs.method, submodels=models, call=cl)    
    class(result) <- "pvs"
    
    return(result)

} ## end pvs.default



pvs.formula <- function(formula, data = NULL,...)
{
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval.parent(m$data))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) 
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0) 
        x <- x[, -xint, drop = FALSE]
    res <- pvs(x, grouping,...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1]] <- as.name("pvs")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        res$na.action <- attr(m, "na.action")
    res
}




predict.pvs <- function(object, newdata, quick = FALSE, detail = FALSE, ...) {
 
    if (!inherits(object, "pvs")) { stop("Object is not of class", " 'pvs'", ".") }
    if (!is.null(terms <- object$terms)) {
        if (missing(newdata)) {
            newdata <- model.frame(object)
            } 
        else {        
            newdata <- model.frame(as.formula(delete.response(terms)), newdata, na.action = function(x) x, xlev = object$xlevels)
            }
        x <- model.matrix(delete.response(terms), newdata, contrasts = object$contrasts)
        xint <- match("(Intercept)", colnames(x), nomatch = 0)
        if (xint > 0) { x <- x[, -xint, drop = FALSE] }
        } 
    else {     
        if (missing(newdata)) {    
            if (!is.null(sub <- object$call$subset)) {
                newdata <- eval.parent(parse(text = paste(deparse(object$call$x, backtick = TRUE), 
                                           "[", deparse(sub, backtick = TRUE), ",]")))
                } 
            else {    
                newdata <- eval.parent(object$call$x)
                }    
            if (!is.null(nas <- object$call$na.action)) { newdata <- eval(call(nas, newdata)) }
            }
        if (is.null(dim(newdata))) { dim(newdata) <- c(1, length(newdata)) }

    x <- as.matrix(newdata)
    }

    m <- diag(0, length(object$classes))
    colnames(m) <- rownames(m) <- object$classes 
    comp.result <- vector("list",dim(x)[1])

    for (i in 1:dim(x)[1]) { comp.result[[i]] <- m }

## for all submodels (i.e. each pair of variables) predict class probabilities for each observation:

    for (i in object$submodels) {
        if (object$method=="svm") {
            pairclass <- try(predict(i$model, as.matrix(x[,i$subspace$variables]), probability=TRUE, ...), silent=TRUE)
            if (is.factor(pairclass)) { pairclass <- attr(pairclass,"probabilities") }
            } 
        else {     
            if (object$method=="randomForest") {
                pairclass <- try(predict(i$model, as.matrix(x[,i$subspace$variables]), type="prob", ...), silent=TRUE)
                } 
            else {
                pairclass <- try(predict(i$model, as.matrix(x[,i$subspace$variables]), ...), silent=TRUE)
                if (is.list(pairclass)) { pairclass <- pairclass$posterior }
                }
            }
    
        coln <- colnames(pairclass)
    
        for (k in 1:dim(x)[1]) {
            for (j in 1:2) { comp.result[[k]][coln[j], coln[3-j]] <- pairclass[k,j] } ##  compares every pair of classes for each observation: rows indicate 'winner',  coloumns 'loser'
            }
        } ## end loop submodels

    if (quick) {    
        posterior <- sapply(comp.result, rowMeans) ## fast, but questionable computation of posterior probabilities
        posterior <- as.matrix(apply(posterior,2, function(x) return(x/sum(x))))
        } 
    else {
        posterior <- sapply(comp.result, pfromp) ## computation of posterior probabilities according to Hastie, Tibshirani
        }

    rownames(posterior) <- object$classes
    classes <- apply(posterior, 2, function(x) return(names(x)[which.max(x)]))

    if (detail) { 
        result <- list(class=factor(classes), posterior=t(posterior), details=comp.result)    
        } 
    else {
        result <- list(class=factor(classes), posterior=t(posterior))
        }

    return(result)
} ## end predict.pvs




#Function for estimating class probabilities from probabilities of pairwise KLOGREG
#see: Hastie,  Tibshirani: Classification by pairwise coupling. In Jordan,  Kearns,  Solla,  editors, 
#Advances in Neural Information Processing Systems,  volume 10. The MIT Press,  1998.
#Author: Marcos Marin-Galiano,  Dept. Of Statistics,  University Of Dortmund,  Germany
#        modified by Gero Szepannek
#email: marcos.marin-galiano@uni-dortmund.de.
#input: matr is the matrix of r_ij = p_i / (p_i + p_j). matr must have the same number of
#rows and columns,  matnum is the matrix of the n_ij

pfromp <- function(matr, tolerance=1.E-4, matnum=matrix(rep(1, dim(matr)[1]^2), nrow=dim(matr)[1]))
{
    classnum<-dim(matr)[1]
    #initial set for the matrix of mu´s.
    matmue<-matrix(0, nrow=classnum, ncol=classnum)
    #initial estimate for class probabilities = laplace probabilities
    pestnew<-rep(1/classnum, classnum)
    wmatr<-diag(matnum%*%t(matr))


    #algorithm loop
    repeat
        {pest<-pestnew

        #computation of new matrix mu and the new estimate for p
        for (i in 2:classnum){
                temp <- pest[i]/(pest[i]+pest[1:(i-1)])
                matmue[i, 1:(i-1)] <- temp
                matmue[1:(i-1), i] <- 1 - temp
        }
        ##normalization of p
        divisor <- diag(matnum%*%t(matmue))
        pestnew<-pest*wmatr/divisor
        if(any(divisor==0)) pestnew[divisor==0] <- 0
        pestnew <- pestnew/sum(pestnew)
        #breaking rule
        if (sum(abs(pestnew-pest))<tolerance) break
    }
    return(pestnew)
}





print.pvs <- function(x, ...){
    cat("Used classifier: ", x$method, "\n\n")
    dummy <- x$vs.method
    if(is.null(dummy)) dummy <- "Kolmogorov Smirnov - test"
    cat("Used variable selection: ", dummy, "\n\n")
    cat("Pairwise subspaces: \n")
    lapply(x$submodels, function(x) cat("classes: ", x$classpair, "\t variable subset: ", x$subspace$variables, "\n"))
    invisible(x)
}





locpvs <- function(x, subclasses, subclass.labels, prior=NULL, method="lda", vs.method=c("ks.test","stepclass","greedy.wilks"), 
                   niveau=0.05, fold=10, impr=0.1, direct="backward", out=FALSE, ...) {
    ## subclass.labels must be a matrix with 2 coloumns: col.1: subclass, col.2: according upper class
    
    vs.method <- match.arg(vs.method)
    result <- pvs(x = x, grouping = subclasses, prior = prior, method = method, vs.method = vs.method, niveau  = niveau, 
          fold = fold, impr = impr, direct = direct, out = out, ...) ## call pvs for the subgroups
    
    result <- list(result,subclass.labels) ## keep the subgroup-group-relationship in the model
    names(result) <- c("pvs.result","subclass.labels") 
    class(result) <- "locpvs"

    return(result)   

}



predict.locpvs <- function(object,newdata,quick=FALSE,return.subclass.prediction=TRUE, ...) {

    prediction <- predict(object$pvs.result,newdata,quick=quick, ...)$posterior
  
## assign to each subclass its upperclass:
  
    for (n in 1:dim(prediction)[2]) {     
        colnames(prediction)[n] <- object$subclass.labels[which(object$subclass.labels[,1]==colnames(prediction)[n]),2]
        }
    classes <- colnames(prediction)

## sum over all subclasses for each class:
    add.subclass.posteriors <- function(x,klassen=classes) { return(by(x,klassen,sum)) }
    
    posteriors <- t(apply(prediction,1,add.subclass.posteriors))
    classes <- as.factor(apply(posteriors, 1, function(x) { return(names(x)[which.max(x)]) } ))
    

    if (return.subclass.prediction) { 
        result <- list(class=factor(classes), posterior=posteriors, subclass.posteriors=prediction)
        } 
    else {
        result <- list(class=factor(classes), posterior=posteriors)
        }

    return(result)
}
