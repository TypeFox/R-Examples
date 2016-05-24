#################################################################
## Workflows
#################################################################
## This module contains functions related to workflows
#################################################################




# =====================================================
# This function generates a named list of objects of class
# learner.
# It is used for easily generating a list of variants of
# a learning system that is usable by the experimentalComparison()
# function.
# If you give only the learning system name it will generate
# a learner wit the default parameters of that system.
# The names of the components are generated automatically
# and have the form <sys>-v<x> where <x> is an increasing
# integer. In the case of using defaults the name is
# <sys>-defaults
# =====================================================
# Luis Torgo, Ago 2013
# =====================================================
# Example call:
# 
# 
#

workflowVariants <- function(wf,...,varsRootName,as.is=NULL) {
    vars <- list(...)
    
    ## if no ID was provided then it is one of the standard workflows
    if (missing(wf)) {
        if ("type" %in% names(vars)) {
            n <- paste(vars[["learner"]],vars[["type"]],sep='.')
            wf <- "timeseriesWF"
        } else {
            n <- vars[["learner"]]
            wf <- "standardWF"
        }
    ## using a user-defined workflow
    } else {
        n <- wf <- wf
    }
    
    if (missing(varsRootName)) varsRootName <- n

    default.novar <- c('evaluator.pars.stats','evaluator.pars.allCls','evaluator.pars.benMtrx')
    allnovar <- c(as.is,default.novar)

    ## unfolding the parameters hidden inside the special parameters


    if (!length(vars)) {
        ##vars <- c(Workflow(wf,list()))
        vars <- c(Workflow(wf))
        names(vars)[1] <- paste(varsRootName,'.v1',sep='')
        vars[[1]]@name <- names(vars)[1]
        return(vars)
    }
  
    ## the special parameters that are list and thus need to be unfolded
    islist <- sapply(vars,function(v) is.list(v))
    spec <- names(vars)[islist]
    ##spec <- c('learner.pars','predictor.pars','evaluator.pars')  
    newVars <- list()
    for(i in 1:length(vars))
        if (names(vars)[i] %in% spec)
            newVars <- c(newVars,unlist(vars[i],recursive=F))
        else
            newVars <- c(newVars,vars[i])
    vars <- newVars

    ## the parameters not involved in variants generation
    ## their names:
    allnovar <- c(allnovar,names(vars)[which(sapply(vars,length)==1)],
                  names(vars)[as.numeric(unlist(sapply(as.is,function(g) grep(g,names(vars)))))])
    ## their positions in vars:
    toExcl <- which(names(vars) %in% allnovar)
    ## their number:
    nExcl <- length(toExcl)
    
    ## checking how many variants per parameter and generate the grid
    nvarsEach <- rep(1,length(vars))
    varying <- if (nExcl) (1:length(vars))[-toExcl] else 1:length(vars)
    if (length(varying)) nvarsEach[varying] <- sapply(vars[varying],length)
    idxsEach <- lapply(nvarsEach,function(x) 1:x)
    theVars <- expand.grid(idxsEach)
    
    ## now go for generating the different variants
    vs <- vector("list",nrow(theVars))
    for(v in 1:nrow(theVars)) {
        ## start
        varPars <- list()
        for(k in 1:ncol(theVars)) {
            if (nExcl & (k %in% toExcl))
                varPars <- c(varPars,vars[k])
            else {
                x <- vars[k]
                x[[1]] <- x[[1]][theVars[v,k]]
                varPars <- c(varPars,x)
            }
        }
        specParsPos <- grep(paste(paste('^',spec,'[:.:]',sep=''),collapse='|'),
                            names(varPars))
        normal <- 1:length(varPars)
        if (length(specParsPos)) normal <- normal[-specParsPos]
        normalPars <- if (length(normal)) varPars[normal] else NULL
        
        finalVars <- list()
        for(i in 1:length(spec)) {
            pos <- grep(paste('^',spec[i],'[:.:]',sep=''),names(varPars))
            if (length(pos)) {
                x <- list()
                for(j in pos) {
                    x <- c(x,list(varPars[[j]]))
                    names(x)[length(x)] <- gsub(paste('^',spec[i],'[:.:]',sep=''),'',names(varPars)[j])
                }
                finalVars <- c(finalVars,list(x))
                names(finalVars)[length(finalVars)] <- spec[i]
            }
        }
        finalVars <- c(finalVars,normalPars)
        ##vs[[v]] <- Workflow(wf,finalVars)
        vs[[v]] <- do.call("Workflow",c(list(wf),finalVars))
    }
    
    for(i in 1:length(vs)) 
        names(vs)[i] <- vs[[i]]@name <- if (wf=="standardWF" || wf=="timeseriesWF") vs[[i]]@pars$learner else varsRootName

    cs <- table(names(vs))
    tochg <- names(cs)[which(cs>1)]
    for(n in tochg) {
        pos <- which(names(vs)==n)
        for(i in seq_along(pos))
            names(vs)[pos[i]] <- vs[[pos[i]]]@name <- paste0(n,".v",i)
    }
    vs
}


# =====================================================
# This function obtains the parameter settings associated
# to a certain variant name in the context of the variants
# of an experimental comparison
# =====================================================
# Luis Torgo, Ago 2013
# =====================================================
# Example Call:
# > getWorkflow('cv.nnet-v6',cvResults)
#
# Note: The result of this can then be "run" as follows,
# > runWorkflow(getWorkflow('cv.nnet-v6',cvResults),
#               medv~.,Boston[1:100,],Boston[-(1:100),])
#
getWorkflow <- function(var,obj)
    obj[[1]][[which(names(obj[[1]]) == var)]]@workflow



# =====================================================
# Function that can be used to call a workflow function
# whose information is stored in an object of class workflow.
# =====================================================
# Luis Torgo, Jul 2013
# =====================================================
# Example run:
# l  <- workflow('nnet',pars=list(size=4,linout=T))
# runWorkflow(l,medv ~ ., Boston)
#
runWorkflow <- function(l,...) {
  if (!inherits(l,'Workflow')) stop(l,' is not of class "Workflow".')
  do.call(l@func,c(list(...),l@pars))
}



## --------------------------------------------------------------
## Standard Workflows
## --------------------------------------------------------------


## =====================================================================
## A function implementing a typical workflow for predictive tasks.
## ---------------------------------------------------------------------
## L. Torgo (Ago, 2013)
##
standardWF <- function(form,train,test,
                       learner,learner.pars=NULL,
                       predictor='predict',predictor.pars=NULL,
                       pre=NULL,pre.pars=NULL,
                       post=NULL,post.pars=NULL,
                       .fullOutput=FALSE)
{
    .fullRes <- if (.fullOutput) list() else NULL

    ## Data pre-processing stage
    if (!is.null(pre)) {
        preprocRes <- do.call("standardPRE",c(list(form,train,test,steps=pre),pre.pars))
        train <- preprocRes$train
        test <- preprocRes$test
        if (.fullOutput) .fullRes$preprocessing <- preprocRes
    }

    ## Learning and prediction stage
    tm <- Sys.time()
    if (is.null(predictor)) {  ## there is no separate predict stage (e.g. kNN)
        ps <- do.call(learner,c(list(form,train,test),learner.pars))
        t.tr <- t.ts <- as.numeric(Sys.time() - tm,units="secs")
        if (.fullOutput) .fullRes$modeling <- ps
    } else {
        m <- do.call(learner,c(list(form,train),learner.pars))
        t.tr <- as.numeric(Sys.time() - tm,units="secs")
        tm <- Sys.time()
        ps <- do.call(predictor,c(list(m,test),predictor.pars))
        t.ts <- as.numeric(Sys.time() - tm,units="secs")
        if (.fullOutput) .fullRes$modeling <- if (is.null(post)) m else list(model=m,initPreds=ps)
    }

    ## Checking for strange things (like regression methods not returning a vector, e.g. earth)
    if (is.numeric(train[,as.character(form[[2]])]) && !is.null(dim(ps))) ps <- ps[,1]
        
    trues <- responseValues(form,test)
    ## Checking for learners that do not ouput as many predictions as test cases!
    ## (e.g. SVM from e1071!)
    if (NROW(ps) != length(trues)) {
        warning("standardWF:: less predictions than test cases, filling with NAs.")
        if (is.null(dim(ps))) {
            t <- trues
            t[] <- NA
            t[names(ps)] <- ps
            ps <- t
        } else {
            t <- matrix(NA,nrow=length(trues),ncol=ncol(ps),dimnames=list(names(trues),colnames(ps)))
            t[names(ps),] <- ps
            ps <- t
        }
    }
    

    ## Data post-processing stage
    if (!is.null(post)) {
        ps <- do.call("standardPOST",c(list(form,train,test,ps,steps=post),post.pars))
        if (.fullOutput) .fullRes$postprocessing <- ps
    }

    ## giving correct names to the predictions
    if (is.null(dim(ps))) names(ps) <- rownames(test)
    else rownames(ps) <- rownames(test)

    ## the final return object (a list) from the workflow
    res <- list(trues=trues,preds=ps,times=c(trainT=t.tr,testT=t.ts))
    if (.fullOutput) res <- c(res,.fullRes) 
    res
}
 

## =====================================================================
## The workhorse function implementing sliding and growing window worflows
## for time series prediction tasks
## ---------------------------------------------------------------------
## L. Torgo (Jan, 2013)
##
timeseriesWF <- function(form,train,test,
                         learner,learner.pars=NULL,
                         type='slide',relearn.step=1,
                         predictor='predict',predictor.pars=NULL,
                         pre=NULL,pre.pars=NULL,
                         post=NULL,post.pars=NULL,
                         .fullOutput=FALSE,verbose=FALSE)
{
    .fullRes <- if (.fullOutput) list() else NULL

    ## Data pre-processing stage
    if (!is.null(pre)) {
        preprocRes <- do.call("standardPRE",c(list(form,train,test,steps=pre),pre.pars))
        train <- preprocRes$train
        test <- preprocRes$test
        if (.fullOutput) .fullRes$preprocessing <- preprocRes
    }
   
    ## Learning and prediction stage
    data <- rbind(train,test)
    n <- NROW(data)
    train.size <- NROW(train)
    sts <- seq(train.size+1,n,by=relearn.step) # relearn times

    preds <- vector()
    if (.fullOutput && !is.null(predictor)) models <- list()

    t.tr <- t.ts <- 0
    
    for(s in sts) { # the learn+test iterations to obtain all predictions

        ## the train and test sets to use on this iteration
        tr <- if (type=='slide') data[(s-train.size):(s-1),] else data[1:(s-1),]
        ts <- data[s:min((s+relearn.step-1),n),]
        
        if (verbose) cat('*')

        tm <- Sys.time()
        if (is.null(predictor)) {
            ps <- do.call(learner,c(list(form,tr,ts),learner.pars))
            t.tr <- t.ts <- t.tr + as.numeric(Sys.time() - tm,units="secs")
            if (.fullOutput) models <- c(models,list(start=s,ps=ps))
        } else {
            m <- do.call(learner,c(list(form,tr),learner.pars))
            t.tr <- t.tr + as.numeric(Sys.time() - tm,units="secs")
            tm <- Sys.time()
            ps <- do.call(predictor,c(list(m,ts),predictor.pars))
            t.ts <- t.ts + as.numeric(Sys.time() - tm,units="secs")
            if (.fullOutput) models <- c(models,list(start=s,model=m,preds=ps))
        }
        ## Checking for strange things (like regression methods not returning a vector, e.g. earth)
        if (is.numeric(tr[,as.character(form[[2]])]) && !is.null(dim(ps))) ps <- ps[,1]
        
        preds <- c(preds,ps)
    }
    if (verbose) cat('\n')
    if (.fullOutput) .fullRes$modeling <- models

    
    trues <- responseValues(form,test)
    ## Checking for learners that do not ouput as many predictions as test cases!
    ## (e.g. SVM from e1071!)
    if (length(preds) != length(trues)) {
        warning("timeseriesWF:: less predictions than test cases, filling with NAs.")
        t <- trues
        t[] <- NA
        t[names(preds)] <- preds
        preds <- t
    }    

    ## Data post-processing stage
    if (!is.null(post)) {
        preds <- do.call("standardPOST",
                         c(list(form,train,test,preds,steps=post),post.pars))
        if (.fullOutput) .fullRes$postprocessing <- preds
    }

    ## giving correct names to the predictions
    names(preds) <- rownames(test)
    
    ## the final return object (a list) from the workflow
    res <- list(trues=trues,preds=preds,times=c(trainT=t.tr,testT=t.ts))
    if (.fullOutput) res <- c(res,.fullRes) 
    res
}




## =====================================================================
## A function implementing some typical pre-processing steps/functions
## ---------------------------------------------------------------------
## L. Torgo, Oct 2014
##
standardPRE <- function(form,train,test,steps,...) {

    tgtVar <- deparse(form[[2]])
    tgtCol <- which(colnames(train)==tgtVar)
    allPreds <- setdiff(colnames(train),tgtVar)
    
    for(s in steps) {
        if (s == "scale") {
            numPreds <- allPreds[sapply(allPreds,function(p) is.numeric(train[[p]]))]
            scaledTrain <- scale(train[,numPreds],...)
            train[,numPreds] <- scaledTrain
            test[,numPreds] <- scale(test[,numPreds],
                                     center=attr(scaledTrain,"scaled:center"),
                                     scale=attr(scaledTrain,"scaled:scale"))
        } else if (s == "centralImp") {
            for (i in allPreds) {
                cval <- if (is.numeric(train[[i]])) median(train[[i]],na.rm=TRUE) else { x <- as.factor(train[[i]]) ; levels(x)[which.max(table(x))] }
                if (any(idx <- is.na(train[[i]]))) train[[i]][idx] <- cval
                if (any(idx <- is.na(test[[i]]))) test[[i]][idx] <- cval
            }
        } else if (s == "knnImp") {
            pars <- list(...)
            if (!("k" %in% names(pars))) pars$k <- 10 # default nr. neigh.
            train[,-tgtCol] <- knnImp(train[,-tgtCol],k=pars$k)
            test[,-tgtCol] <- knnImp(test[,-tgtCol],k=pars$k,distData=train[,-tgtCol])
        } else if (s == "na.omit") {
            train <- na.omit(train)
            test <- na.omit(test)
        } else if (s == "undersample") {
            if (is.numeric(train[,tgtVar])) stop("Undersampling is currently only available for classification tasks. Check http://www.dcc.fc.up.pt/~ltorgo/ExpertSystems/ for approaches applicable to regression.",call.=FALSE)
            pars <- list(...)
            clDistr <- table(train[,tgtVar])
            minCl <- which.min(clDistr)
            minClName <- names(minCl)
            nMin <- min(clDistr)
            minExs <- which(train[,tgtVar] == minClName)
            if (!("perc.under" %in% names(pars))) pars$perc.under <- 1 # default under %
            selMaj <- sample((1:NROW(train))[-minExs],
                             as.integer(pars$perc.under*nMin),
                             replace=TRUE)
            train <- train[c(minExs,selMaj),]
        } else if (s == "smote") {
            if (is.numeric(train[,tgtVar])) stop("SMOTE is currently only available for classification tasks. Check http://www.dcc.fc.up.pt/~ltorgo/ExpertSystems/ for approaches applicable to regression.",call.=FALSE)
            train <- smote(form,train,...)
        } else {
            user.pre <- do.call(s,c(list(form,train,test),list(...)))
            train <- user.pre$train
            test  <- user.pre$test
        }
    }

    list(train=train,test=test)
}



# =====================================================
# Luis Torgo, Mar 2009, Nov 2011
# =====================================================
knnImp <- function(data,k=10,scale=TRUE,distData=NULL) {
    
    n <- nrow(data)  
    if (!is.null(distData)) {
        distInit <- n+1
        data <- rbind(data,distData)
    } else distInit <- 1
    N <- nrow(data)
    
    ncol <- ncol(data)
    nomAttrs <- rep(F,ncol)
    for(i in seq(ncol)) nomAttrs[i] <- is.factor(data[,i])
    nomAttrs <- which(nomAttrs)
    hasNom <- length(nomAttrs)
    contAttrs <- setdiff(seq(ncol),nomAttrs)
    
    dm <- data
    if (scale) dm[,contAttrs] <- scale(dm[,contAttrs])
    if (hasNom)
        for(i in nomAttrs) dm[,i] <- as.integer(dm[,i])
    
    dm <- as.matrix(dm)
    
    nas <- which(!complete.cases(dm))
    if (!is.null(distData)) tgt.nas <- nas[nas <= n]
    else tgt.nas <- nas
    
    if (length(tgt.nas) == 0)
        warning("No case has missing values. Stopping as there is nothing to do.")
    
    xcomplete <- dm[setdiff(distInit:N,nas),]
    if (nrow(xcomplete) < k)
        stop("Not sufficient complete cases for computing neighbors.",call.=FALSE)
    
    for (i in tgt.nas) {
        
        tgtAs <- which(is.na(dm[i,]))
        
        dist <- scale(xcomplete,dm[i,],FALSE)
        
        xnom <- setdiff(nomAttrs,tgtAs)
        if (length(xnom)) dist[,xnom] <-ifelse(dist[,xnom]>0,1,dist[,xnom])
        
        dist <- dist[,-tgtAs]
        dist <- sqrt(drop(dist^2 %*% rep(1,ncol(dist))))
        ks <- order(dist)[seq(k)]
        for(j in tgtAs) {
            vals <- data[setdiff(distInit:N,nas),j][ks]
            data[i,j] <- if (is.numeric(vals)) median(vals,na.rm=TRUE) else { x <- as.factor(vals) ; levels(x)[which.max(table(x))] }
        }
    }
    data[1:n,]
}


## =====================================================================
## A function implementing some typical pre-processing steps/functions
## ---------------------------------------------------------------------
## L. Torgo, Oct 2014
##
standardPOST <- function(form,train,test,preds,steps,...) {

    tgtVar <- deparse(form[[2]])
    allPreds <- setdiff(colnames(train),tgtVar)
    
    for(s in steps) {
        ## -----------
        if (s == "na2central") {
            if (!is.vector(preds)) stop("standardPOST:: 'na2central' is only applicable to predictions that are vectors of values.",call.=FALSE)
            if (any(idx <- is.na(preds))) {
                cval <- if (is.numeric(train[[tgtVar]])) median(train[[tgtVar]],na.rm=TRUE) else { x <- as.factor(train[[tgtVar]]) ; levels(x)[which.max(table(x))] }
                preds[idx] <- cval
            }
            
        ## -----------
        } else if (s == "onlyPos") {
            if (!is.vector(preds)) stop("standardPOST:: 'onlyPos' is only applicable to predictions that are vectors of values.",call.=FALSE)
            if (! is.numeric(train[[tgtVar]])) stop("standardPOST:: 'onlyPos' is only applicable to numeric predictions.",call.=FALSE)
            if (any(idx <- preds < 0)) preds[idx] <- 0

        ## -----------
        } else if (s == "cast2int") {
            if (!is.vector(preds)) stop("standardPOST:: 'cast2int' is only applicable to predictions that are vectors of values.",call.=FALSE)
            if (! is.numeric(train[[tgtVar]])) stop("standardPOST:: 'cast2int' is only applicable to numeric predictions.",call.=FALSE)
            pars <- list(...)
            if (any(idx <- preds < pars$infLim)) preds[idx] <- pars$infLim
            if (any(idx <- preds > pars$supLim)) preds[idx] <- pars$supLim

        ## -----------    
        } else if (s == "maxutil") {
            pars <- list(...)
            if (is.numeric(train[,tgtVar])) stop("'maxutil'' is only available for classification tasks.",call.=FALSE)
            if (!("cb.matrix" %in% names(pars))) stop("'maxutil' requires that you specify a cost-benefit matrix.",call.=FALSE)
            if (is.null(dim(preds))) stop("'maxutil' requires that the classifier outputs a matrix of probabilities as predictions.",call.=FALSE)
            if (ncol(preds) != ncol(pars$cb.matrix)) stop("Error in 'maxutil': predictions do not contain as many class probabilities as there are classes in the cost-benefit matrix.",call.=FALSE)
            ps <- apply(preds,1,function(ps) which.max(apply(pars$cb.matrix,2,function(cs) sum(ps*cs))))
            preds <- levels(train[,tgtVar])[ps]
            
        ## -----------    
        } else {
            preds <- do.call(s,c(list(form,train,test,preds),...))
        }
    }

    preds
}
 
