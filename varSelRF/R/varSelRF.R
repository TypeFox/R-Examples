#### watch out use of "..." inside Boot function

#####   in bootSimplifyRF maybe allow allBootRuns = NULL so that object is small.



varSelRF <- function(xdata, Class,
                     c.sd = 1,
                     mtryFactor = 1,
                     ntree = 5000,
                     ntreeIterat = 2000,                     
                     vars.drop.num = NULL,
                     vars.drop.frac = 0.2,
                     whole.range = TRUE,
                     recompute.var.imp = FALSE,
                     verbose = FALSE,
                     returnFirstForest = TRUE,
                     fitted.rf = NULL,
                     keep.forest = FALSE) {

    if(!is.factor(Class))
        stop("Class should be a factor")
    if( (is.null(vars.drop.num) & is.null(vars.drop.frac)) |
       (!is.null(vars.drop.num) & !is.null(vars.drop.frac)))
        stop("One (and only one) of vars.drop.frac and vars.drop.num must be NULL and the other set")
  
    max.num.steps <- dim(xdata)[2]
    num.subjects <- dim(xdata)[1]

    if(is.null(colnames(xdata)))
        colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep ="")
    
    ##oversize the vectors; will prune later.
    n.vars <- vars <- OOB.rf <- OOB.sd <- rep(NA, max.num.steps)

    oobError <- function(rf) { ## should not be exported in the namespace.
        ## The out of bag error
        ooo <- rf$confusion[, -dim(rf$confusion)[2]]
        s.ooo <- sum(ooo)
        diag(ooo) <- 0
        sum(ooo)/s.ooo
    }

    
    if(!is.null(fitted.rf)) {
        if(ncol(fitted.rf$importance) < 2)
            stop("The fitted rf was not fitted with importance = TRUE")
        n.ntree <- fitted.rf$ntree
        mtry <- fitted.rf$mtry
        n.mtryFactor <- mtry/sqrt(ncol(xdata))
        if((n.ntree != ntree) | (n.mtryFactor != mtryFactor))
            warning("Using as ntree and mtry the parameters obtained from fitted.rf",
                    immediate.= TRUE)
        ntree <- n.ntree
        mtryFactor <- n.mtryFactor
        rm(n.ntree, n.mtryFactor)
        rf <- fitted.rf
    } else {
        mtry <- floor(sqrt(ncol(xdata)) * mtryFactor)
        rf <- randomForest(x = xdata, y = Class,
                           ntree = ntree, mtry = mtry,
                           importance = TRUE,
                           keep.forest = keep.forest)
    }
    
    if(returnFirstForest)
        FirstForest <- rf
    else
        FirstForest <- NULL
    m.iterated.ob.error <- m.initial.ob.error <- oobError(rf)
    sd.iterated.ob.error <- sd.initial.ob.error <-
        sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) *
             (1/num.subjects))

    if(verbose) {
        print(paste("Initial OOB error: mean = ",
                    round(m.initial.ob.error, 4),
                    "; sd = ", round(sd.initial.ob.error, 4), sep = ""))
    }

    importances <- importance(rf, type = 1, scale = FALSE)
    selected.vars <- order(importances, decreasing = TRUE)
    ordered.importances <- importances[selected.vars]
    
    initialImportances <- importances
    initialOrderedImportances <- ordered.importances
    
    j <- 1
    n.vars[j] <- dim(xdata)[2] 
    vars[j] <- paste(colnames(xdata), collapse = " + ")
    OOB.rf[j] <- m.iterated.ob.error
    OOB.sd[j] <- sd.iterated.ob.error

    var.simplify <- TRUE
    
    while(var.simplify) {
        if (verbose){
          print("gc inside loop of varSelRF")
          print(gc())
        } else {
          gc()
        }
        
        last.rf <- rf 
        last.vars <- selected.vars
        previous.m.error <- m.iterated.ob.error
        previous.sd.error <- sd.iterated.ob.error

        ## If this is left as only
        ## "if(length(selected.vars) <= 2) var.simplify <- FALSE"
        ## under the if((length(selected.vars) < 2) | (any(selected.vars < 1))),
        ## as it used to be, then we fit a 2 model var, which might be
        ## better or as good as others, but as we never re-enter,
        ## we cannot return it, even if we see it in the history.

        ## Alternatively, we cannot just convert
        ## "if((length(selected.vars) < 2) | (any(selected.vars < 1))"
        ## to the <= 2, as we then would re-enter many times because
        ## of the way selected.vars <- selected.vars[1:2] when
        ## num.vars < (vars.drop + 2)

        ## This way, we enter just to set last.rf, last.vars and
        ## we bail out
        
        if(length(selected.vars) <= 2) {
          var.simplify <- FALSE
          break
        }

        
        if(recompute.var.imp & (j > 1)) {
            importances <- importance(rf, type = 1, scale = FALSE)
            tmp.order <- order(importances, decreasing = TRUE)
            selected.vars <- selected.vars[tmp.order]
            ordered.importances <- importances[tmp.order]
        }
        
        num.vars <- length(selected.vars)
  
        if(is.null(vars.drop.num))
            vars.drop <- round(num.vars * vars.drop.frac)
        else vars.drop <- vars.drop.num
        
        if(num.vars >= (vars.drop + 2)) {
            selected.vars <- selected.vars[1: (num.vars - vars.drop)]
            ordered.importances <- ordered.importances[1: (num.vars - vars.drop)]
        }
        else {
            selected.vars <- selected.vars[1:2]
            ordered.importances <- ordered.importances[1:2]
        }
        
        ## couldn't we eliminate the following?
        if((length(selected.vars) < 2) | (any(selected.vars < 1))) {
            var.simplify <- FALSE
            break
        }
        
        
        
        mtry <- floor(sqrt(length(selected.vars)) * mtryFactor)
        if(mtry > length(selected.vars)) mtry <- length(selected.vars)
        
        if(recompute.var.imp) 
            rf <- randomForest(x = xdata[, selected.vars],
                               y = Class, importance= TRUE,
                               ntree = ntree, mtry = mtry,
                               keep.forest = keep.forest)
        else
            rf <- randomForest(x = xdata[, selected.vars],
                               y = Class, importance= FALSE,
                               ntree = ntreeIterat, mtry = mtry,
                               keep.forest = keep.forest)
        
        m.iterated.ob.error <- oobError(rf)
        sd.iterated.ob.error <-
            sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) *
                 (1/num.subjects))
        
        if(verbose) {
            print(paste("..... iteration ", j, "; OOB error: mean = ",
                        round(m.iterated.ob.error, 4),
                        "; sd = ", round(sd.iterated.ob.error, 4),
                        "; num. vars = ", length(selected.vars), 
                        sep = ""))
        }
        j <- j + 1
        
        
        n.vars[j] <- length(selected.vars)
        vars[j] <- paste(colnames(xdata)[selected.vars],
                         collapse = " + ")
        OOB.rf[j] <- m.iterated.ob.error
        OOB.sd[j] <- sd.iterated.ob.error


        if(!whole.range &
           (
            (m.iterated.ob.error >
             (m.initial.ob.error + c.sd*sd.initial.ob.error))
            |
            (m.iterated.ob.error >
             (previous.m.error + c.sd*previous.sd.error)))
           )
            var.simplify <- FALSE
    }

    if (!whole.range) {
        if(!is.null(colnames(xdata)))
            selected.vars <- sort(colnames(xdata)[last.vars])
        else
            selected.vars <- last.vars

        out <- list(selec.history = data.frame(
                    Number.Variables = n.vars,
                    Vars.in.Forest = vars,
                    OOB = OOB.rf,
                    sd.OOB = OOB.sd)[1:j,],
                    rf.model = last.rf,
                    selected.vars = selected.vars,
                    selected.model =  paste(selected.vars, collapse = " + "),
                    best.model.nvars = length(selected.vars),
                    initialImportances = initialImportances,
                    initialOrderedImportances = initialOrderedImportances,
                    ntree = ntree,
                    ntreeIterat = ntreeIterat,
                    mtryFactor = mtryFactor,
#                    mtry = mtry,
                    firstForest = FirstForest)
        class(out) <- "varSelRF"
        return(out)
    }
    else {
      ## Prune the too long vectors created at begin.
      ## not needed above, because we select the 1:j rows
      ## of the return matrix selec.history.
        n.vars <- n.vars[1:j]
        vars <- vars[1:j]
        OOB.rf<- OOB.rf[1:j]
        OOB.sd <- OOB.sd[1:j]
        min.oob.ci <- min(OOB.rf) + c.sd * OOB.sd[which.min(OOB.rf)]
        best.pos <-
            which(OOB.rf <= min.oob.ci)[which.min(n.vars[which(OOB.rf <= min.oob.ci)])]

        selected.vars <- sort(unlist(strsplit(vars[best.pos],
                                              " + ", fixed = TRUE)))
        out <- list(selec.history = data.frame(
                    Number.Variables = n.vars,
                    Vars.in.Forest = vars,
                    OOB = OOB.rf,
                    sd.OOB = OOB.sd),
                    rf.model = NA,
                    selected.vars = selected.vars,
                    selected.model = paste(selected.vars, collapse = " + "),
                    best.model.nvars = n.vars[best.pos],
                    initialImportances = initialImportances,
                    initialOrderedImportances = initialOrderedImportances,
                    ntree = ntree,
                    ntreeIterat = ntreeIterat,
                    mtryFactor = mtryFactor,
                    ##mtry = mtry,
                    firstForest = FirstForest)
        class(out) <- "varSelRF"
        return(out)
    }
}


print.varSelRF <- function(x, ...) {
    cat("\nBackwards elimination on random forest; ")
  cat(paste("ntree = ", x$ntree,";  mtryFactor = ",
            x$mtryFactor, "\n"), sep ="")
  cat("\n Selected variables:\n")
  print(x$selected.vars)
  cat("\n Number of selected variables:", x$best.model.nvars, "\n\n")
}

#summary.varSelRF <- print.varSelRF

plot.varSelRF <- function(x, nvar = NULL, which = c(1, 2), ...) {
    if (length(which) == 2 && dev.interactive()) {
        op <- par(ask = TRUE, las = 1)
    } else {
        op <- par(las = 1)
    }
    
    on.exit(par(op))
    
    if(is.null(nvar))
        nvar <- min(30,
                    length(x$initialOrderedImportances))
    
    show <- c(FALSE, FALSE)
    show[which] <- TRUE

    if (show[1]){
      dotchart(rev(x$initialOrderedImportances[1:nvar]),
          main = "Initial importances",
          xlab = "Importances (unscaled)")
    }
    
    if (show[2]){
      ylim <- c(0, max(0.50, x$selec.history$OOB))
      plot(x$selec.history$Number.Variables,
          x$selec.history$OOB, type = "b",
          xlab = "Number of variables used",
          ylab = "OOB error", log = "x",
          ylim = ylim,
          ...)
      
      ##  if(max(x$selec.history$Number.Variables) > 300)
      ##      axis(1, at = c(1, 2, 3, 5, 8, 15, 25, 50, 75, 150, 200, 300),
      ##           labels = c(1, 2, 3, 5, 8, 15, 25, 50, 75, 150, 200, 300))
      
      lines(x$selec.history$Number.Variables,
          x$selec.history$OOB +
              2 * x$selec.history$sd.OOB, lty = 2)
      lines(x$selec.history$Number.Variables,
          x$selec.history$OOB -
              2 * x$selec.history$sd.OOB, lty = 2)
    }
}
  
## We could also write a varSelRFCV zz: move al TODO

varSelRFBoot <- function(xdata, Class,
                         c.sd = 1,
                         mtryFactor = 1,
                         ntree = 5000,
                         ntreeIterat = 2000,
                         vars.drop.frac = 0.2,
                         bootnumber = 200,
                         whole.range = TRUE,
                         recompute.var.imp = FALSE,
                         usingCluster = TRUE,
                         TheCluster = NULL,
                         srf = NULL,
                         verbose = TRUE,
                         ...) {
    
    ## beware there is a lot of data copying... pass by reference, or
    ## minimize copying or something.
    

    if(is.null(colnames(xdata)))
        colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep ="")
    
    if(!is.null(srf)) { ## we are passing a simplified rf object.
        if(class(srf) != "varSelRF")
            stop("srf must be the results of a run of varSelRF")
        n.ntree <- srf$ntree
        n.ntreeIterat <- srf$ntreeIterat
        n.mtryFactor <- srf$mtryFactor
        if((n.ntree != ntree) | (n.mtryFactor != mtryFactor) |
           (n.ntreeIterat != ntreeIterat))
        warning("Using as ntree and mtryFactor the parameters obtained from srf",
                immediate.= TRUE)
        ntree <- n.ntree
        mtryFactor <- n.mtryFactor
        rm(n.ntree, n.mtryFactor)
        all.data.run <- srf
    } else { ## we are simplifying the random forest
        all.data.run <- varSelRF(Class = Class,
                             xdata = xdata,
                             c.sd = c.sd,
                             mtryFactor = mtryFactor,
                             ntree = ntree,
                             ntreeIterat = ntreeIterat,
                             vars.drop.frac = vars.drop.frac,
                             whole.range = whole.range,
                             recompute.var.imp = recompute.var.imp)
###                             ...) 
    }
    columns.data <- which(colnames(xdata) %in%
                          all.data.run$selected.vars)
    
    all.data.rf.mtry <- floor(mtryFactor * sqrt(length(columns.data)))
    if(all.data.rf.mtry > length(columns.data))
        all.data.rf.mtry <- length(columns.data)
    all.data.rf.predict <- randomForest(y = Class,
                                        x = xdata[, columns.data],
                                        ntree = all.data.run$ntree,
                                        mtry = all.data.rf.mtry,
                                        xtest = xdata[, columns.data],
                                        ytest = Class,
                                        keep.forest = FALSE)
    
    full.pred <- all.data.rf.predict$test$predicted
    
    all.data.selected.vars <- all.data.run$selected.vars
    ##    all.data.selected.model <- all.data.run$selected.model
    all.data.best.model.nvars <- all.data.run$best.model.nvars
    
    N <- length(Class)
    solution.sizes <- rep(NA, bootnumber)
    overlap.with.full <- rep(NA, bootnumber)
    vars.in.solutions <- vector()
##    solutions <- rep(NA, bootnumber)
    
    bootTrainTest <- function(dummy,
        xdataTheCluster,
        ClassTheCluster,
        c.sd,
        mtryFactor,
        ntree, ntreeIterat,
        whole.range, recompute.var.imp,
        ...) {
        N <- length(ClassTheCluster)
        sample.again <- TRUE
        while(sample.again) {
            bootsample <- unlist(tapply(1:N,  ClassTheCluster,
                                        function(x)
                                        sample(x, size = length(x),
                                               replace = TRUE)))
            ## sure, this isn't the fastest, but will do for now.
            nobootsample <- setdiff(1:N, bootsample)
            if(!length(nobootsample))
                sample.again <- TRUE
            else sample.again <- FALSE
        }
        ## this is an ugly hack to prevent nobootsamples
        ## of size 0.
        
        train.data <- xdataTheCluster[bootsample, , drop = FALSE]
        test.data <- xdataTheCluster[nobootsample, , drop = FALSE]
        train.class <- ClassTheCluster[bootsample]
        test.class <- ClassTheCluster[nobootsample]

        boot.run <- varSelRF(Class = train.class,
                             xdata = train.data,
                             c.sd = c.sd,
                             mtryFactor = mtryFactor,
                             ntree = ntree,
                             ntreeIterat = ntreeIterat,
                             whole.range = whole.range,
                             recompute.var.imp = recompute.var.imp,
                             vars.drop.frac = vars.drop.frac)
###                             ...)
        output.cl <- list()
        output.cl$best.model.nvars <- boot.run$best.model.nvars
        output.cl$selected.model <- boot.run$selected.model
        output.cl$selected.vars <- boot.run$selected.vars
        output.cl$nobootsample <- nobootsample
        output.cl$bootsample <- bootsample
        output.cl$initialImportances <- boot.run$initialImportances
        output.cl$initialOrderedImportances <- boot.run$initialOrderedImportances
        output.cl$selec.history <- boot.run$selec.history

        boot.col.data <- which(colnames(xdataTheCluster) %in%
                               boot.run$selected.vars)
        run.test.mtry <- floor(mtryFactor * sqrt(length(boot.col.data)))
        if(run.test.mtry > length(boot.col.data))
            run.test.mtry <- length(boot.col.data)

        boot.run.test <- randomForest(y  = train.class,
                                      x = train.data[, boot.col.data,
                                      drop = FALSE],
                                       ntree = boot.run$ntree,
                                      mtry = run.test.mtry,
                                      keep.forest = FALSE,
                                      xtest = test.data[, boot.col.data,
                                      drop = FALSE])$test
##                                      ytest = test.class)$test
        output.cl$class.pred.array <- boot.run.test$predicted
        output.cl$prob.pred.array <- boot.run.test$votes
        rm(boot.run.test)
        gc()
        return(output.cl)
    } ## bootTrainTest; this function is defined to be used in the cluster.
    ## we use it too if not in the cluster. There is a bit too much
    ## copying, because we define a pair of new objects
    ## xdataTheCluster and ClassTheCluster that need not be defined
    ## if we just copied the code. But that is ugly and hard to follow,
    ## and the non-cluster is slow for other reasons.
    
    if(usingCluster) {
        if (verbose){
          print("gc inside varSelRFBoot papply")
          print(gc())
        } else {
          gc()
        }
        
        cat("\n      Running bootstrap iterations using cluster (can take a while)\n")
        boot.runs <- clusterApplyLB(TheCluster,
                                    1:bootnumber,
                                    bootTrainTest,
                                    xdataTheCluster = xdata,
                                    ClassTheCluster = Class,
                                    c.sd = c.sd,
                                    mtryFactor = mtryFactor,
                                    ntree = ntree,
                                    ntreeIterat = ntreeIterat,
                                    whole.range = whole.range,
                                    recompute.var.imp = recompute.var.imp)
        ## clean up before leaving
        clusterEvalQ(TheCluster,
                     rm(list = c("xdataTheCluster", "ClassTheCluster")))

    } else { ## Not using Cluster
        boot.runs <- list()
        xdataTheCluster <- xdata
        ClassTheCluster <- Class
        cat("\n      Running bootstrap iterations")
        for(nboot in 1:bootnumber) {
            cat(".")
            boot.runs[[nboot]] <-
                bootTrainTest(nboot,
                              xdataTheCluster = xdata,
                              ClassTheCluster = Class,
                              c.sd = c.sd,
                              mtryFactor = mtryFactor,
                              ntree = ntree,
                              ntreeIterat = ntreeIterat,
                              whole.range = whole.range,
                              recompute.var.imp = recompute.var.imp)
        }
        cat("\n")
    }
        

    solutions <- unlist(lapply(boot.runs, function(z) {
        paste(sort(z$selected.vars), collapse = " + ")}))
    vars.in.solutions <- unlist(lapply(boot.runs,
                                       function(z) z$selected.vars))
    solution.sizes <- unlist(lapply(boot.runs,
                                    function(z) z$best.model.nvars))
    overlap.with.full <-
        unlist(lapply(boot.runs,
                      function(x) {
                          length(intersect(x$selected.vars,
                                           all.data.selected.vars))/
                                               sqrt(all.data.best.model.nvars *
                                                    x$best.model.nvars)})) 

    prob.pred.array <-
        array(NA, dim = c(N,  nlevels(Class), bootnumber),
              dimnames = list(1:N, levels(Class),
              paste("BootstrapReplication.", 1:bootnumber, sep = "")))

    ## to store class predictions as data frame
    class.pred.array <- data.frame(Class)
    
    for(nb in 1:bootnumber) {
        nobootsample <- boot.runs[[nb]]$nobootsample
        class.pred.array[[nb]] <- factor(NA, levels = levels(Class))
        class.pred.array[nobootsample, nb] <-
            boot.runs[[nb]]$class.pred.array
        prob.pred.array[nobootsample, , nb] <-
            boot.runs[[nb]]$prob.pred.array
    }
    names(class.pred.array) <-
        paste("BootstrapReplication.", 1:bootnumber, sep = "")

            
    ############## The .632+ estimate of prediction error ##################
    ##      one:         \hat{Err}^{(1)} in Efron & Tibshirani, 1997, p. 550,
    ##                   or leave-one-out bootstrap error.  
    ##      resubst:     \bar{err}, the apparent error rate, resubstitution rate,
    ##                   or "in sample" error rate.
    ##      full.pred:   the "in sample" prediction from full model; only used
    ##                   to obtain gamma.
    ##      gamma:       gamma in p. 552
    ##      r:           R in p. 552 (bounding in [0, 1]).
    ##      err632:      the .632 error
    ##      errprime:    the one prime; \hat{Err}^{(1)}' 
    ##      err:         the .632+ error
    
    ## This is how I find one
    
    one <- mean(apply(cbind(class.pred.array, Class), 1,
                      function(x) {mean(x[-(bootnumber + 1)] != x[bootnumber + 1],
                                        na.rm = TRUE)}), na.rm = TRUE)
    ## this is equivalent to what Torsten Hothorn does in ipred
            
    resubst <- mean(full.pred != Class)
    ## The following code I take directly from the function
    ## bootest.factor, by Torsten Hothorn, in package ipred.
    err632 <- 0.368 * resubst + 0.632 * one
    gamma <-
        sum(outer(as.numeric(Class),
                  as.numeric(full.pred),
                  function(x, y) ifelse(x == y, 0, 1)))/(length(Class)^2)
    
    r <- (one - resubst)/(gamma - resubst)
    r <- ifelse(one > resubst & gamma > resubst, r, 0)
    if((r > 1) | (r < 0)) { ## just debugging; eliminar m�s adelante
        print(paste("r outside of 0, 1 bounds: one", one,
                    "resubst", resubst, "gamma", gamma))
        if(r > 1) {
            r <- 1
            print("setting r to 1")
        }
        else if(r < 0) {
            r <- 0
            print("setting r to 0")
        }
    }
    
    errprime <- min(one, gamma)
    err <- err632 + (errprime - resubst) *
        (0.368 * 0.632 * r)/(1 - 0.368 * r)
    cat("\n     .632+ prediction error ", round(err, 4), "\n")

    out <- list(number.of.bootsamples = bootnumber, 
                bootstrap.pred.error = err,
                resubstitution.error = resubst,
                leave.one.out.bootstrap.error = one,
                all.data.randomForest = all.data.rf.predict,
                all.data.vars = all.data.selected.vars,
                all.data.run = all.data.run,
                class.predictions = class.pred.array,
                prob.predictions = prob.pred.array,
                number.of.vars = solution.sizes,
                overlap = overlap.with.full,
                all.vars.in.solutions = vars.in.solutions,
                all.solutions = solutions,
                class = Class,
                allBootRuns = boot.runs)
    class(out) <- "varSelRFBoot"
    return(out)
}


print.varSelRFBoot <- function(x, ...) {
    cat("\n\n Variable selection with random forest \n")
    cat(" ------------------------------\n")
##    cat("\n\n randomForest summary \n")
##    print(object$all.data.randomForest)
    cat("\n Variables used \n")
    print(x$all.data.vars)
    cat("\n \n Number of variables used: ", length(x$all.data.vars),
        "\n")
    
    cat("\n\n Bootstrap results\n")
    cat(" ------------------\n")
    cat("\n Bootstrap (.632+) estimate of prediction error: \n",
        " (using", x$number.of.bootsamples, "bootstrap iterations): \n",
        x$bootstrap.pred.error)
    cat("\n\n Number of vars in bootstrapped forests:        \n")
    print(summary(x$number.of.vars))
}


## zz: pass gene names and subject names??
## look at examples from ~/Proyectos/Signatures/Symposum/boot.pamr.knn.dlda.R
## o similar 
## Or not, because I think I am ussing column names
## but look at that code anyway for improvements


#### in summaryBoot: a plot of error rate vs. 
#### number of variables
#### requires saving el "all.data.run" en el Boot, which should ALWAYS be done!!


summary.varSelRFBoot <- function(object,
                                 return.model.freqs = FALSE,
                                 return.class.probs = TRUE,
                                 return.var.freqs.b.models = TRUE,
                                 ...) {
    cat("\n\n Variable selection using all data \n")
    cat(" ------------------------------\n")
###    cat("\n\n randomForest summary \n")
###    print(object$all.data.randomForest)
    cat("\n \n variables used \n")
    print(object$all.data.vars)
    cat("\n \n Number of variables used: ", length(object$all.data.vars),
        "\n")
    cat("\n\n Bootstrap results\n")
    cat(" ------------------\n")
    cat("\n\n Bootstrap (.632+) estimate of prediction error:",
        object$bootstrap.pred.error, " (using",
        object$number.of.bootsamples, "bootstrap iterations).\n")
    cat("\n\n Resubstitution error:                          ",
        object$resubstitution.error, "\n")
    cat("\n\n Leave-one-out bootstrap error:                 ",
        object$leave.one.out.bootstrap.error, "\n")
    nk <- as.vector(table(object$class))
    cat("\n\n Error rate at random:                          ",
        1 - (max(nk)/sum(nk)), "\n")
    cat("\n\n Number of vars in bootstrapped forests:        \n")
    print(summary(object$number.of.vars))
    cat
    cat("\n")
    cat("\n Overlapp of bootstrapped forests with forest from all data\n")
    print(summary(object$overlap))
    if(return.var.freqs.b.models) {
        cat("\n\n Variable freqs. in bootstrapped models \n")
        print(sort(table(object$all.vars.in.solutions),
                   decreasing = TRUE)/object$number.of.bootsamples)
    }
    in.all.data <-
        which(names(table(object$all.vars.in.solutions)) %in% object$all.data.vars)
    cat("\n\n Variable freqs. of variables in forest from all data, and summary \n")
    print(table(object$all.vars.in.solutions)[in.all.data]/object$number.of.bootsamples)
    cat("\n")
    print(summary(table(object$all.vars.in.solutions)[in.all.data]/object$number.of.bootsamples))


    if(return.model.freqs) {
        tmp.table <- sort(table(object$all.solutions),
                          decreasing = TRUE)/object$number.of.bootsamples
        n.tmp.table <- names(tmp.table)
        dim(tmp.table) <- c(dim(tmp.table), 1)
        rownames(tmp.table) <- n.tmp.table
        colnames(tmp.table) <- "Freq."
        cat("\n\n Solutions frequencies in bootstrapped models \n")
        print(tmp.table)
    }
    
    if(return.class.probs)
        cat("\n\n Mean class membership probabilities from out of bag samples\n")
    if(return.class.probs) {
        mean.class.probs <- apply(object$prob.predictions, c(1, 2),
                                  function(x) mean(x, na.rm = TRUE))
        colnames(mean.class.probs) <- levels(object$class)
    }
    if(return.class.probs) return(mean.class.probs)
}



plot.varSelRFBoot <- function(x,
                              oobProb = TRUE,
                              oobProbBoxPlot = FALSE,
                              ErrorNum = TRUE,
                              subject.names = NULL,
                              class.to.plot = NULL,
                              ...) {
    
    if(oobProb | oobProbBoxPlot) {
        mean.class.probs <- apply(x$prob.predictions, c(1, 2),
                                  function(x) mean(x, na.rm = TRUE))
        colnames(mean.class.probs) <- levels(x$class)
        rainbow.col <- rainbow(nlevels(x$class))
	if(dev.interactive()) {
              op <- par(ask = TRUE, las = 1)
    	} else {
               op <- par(las = 1)
    	}
        on.exit(par(op))
####        for(i in 1:ncol(mean.class.probs)) {
        if(is.null(class.to.plot))
            class.to.plot <- 1:ncol(mean.class.probs)
        for(i in class.to.plot) {
            if(oobProbBoxPlot)  {
                boxplot(data.frame(t(x$prob.predictions[ , i, ])),
                        xlab = "Samples",
                        ylab = "Out of bag probability of membership",
                        main = paste("Class", colnames(mean.class.probs)[i]),
                        type = "p",
                        col = rainbow.col[as.numeric(x$class)],
                        pch = 19, ylim = c(0, 1.3),
                        axes = FALSE)
            } else {
###                dotchart(mean.class.probs[, i],
###                         xlab = "Samples",
###                         ylab = "(Average) Out of bag probability of membership",
###                         main = paste("Class", colnames(mean.class.probs)[i]),
###                         type = "p",
###                         col = rainbow.col[as.numeric(x$class)],
###                         pch = 19, ylim = c(0, 1.2),
###                         axes = FALSE)
                plot(mean.class.probs[, i],
                     xlab = "Samples",
                     ylab = "(Average) Out of bag probability of membership",
                     main = paste("Class", colnames(mean.class.probs)[i]),
                     type = "p",
                     col = rainbow.col[as.numeric(x$class)],
                     pch = 19, ylim = c(0, 1.3),
                     axes = FALSE)
            }
            if(!is.null(subject.names))
                text(mean.class.probs[ , i],
                     labels = subject.names, pos = 2, cex = 0.8)
            
            box()
            ##axis(1)
            par(las = 2)
            axis(1, at = seq(1:dim(mean.class.probs)[1]),
                 labels = seq(1:dim(mean.class.probs)[1]))
            segments(x0 = seq(1:dim(mean.class.probs)[1]),
                     y0 = 0,
                     x1 = seq(1:dim(mean.class.probs)[1]),
                     y1 = 1,
                     lty = 2,
                     col = "grey")
            axis(2, at = seq(from = 0, to = 1, by = 0.2))
            legend(y = c(1.29, 1.15), x = c(1, dim(mean.class.probs)[1]),
                   legend = paste("Class", levels(x$class)),
                   col = rainbow.col,
                   pch = 19, bty = "n" )
            abline( h = 1.05, lty = 1)
            abline(h = seq(from = 0, to = 1,
                   length = 1 + nlevels(x$class))[-c(1, nlevels(x$class) + 1)], 
                   lty = 2)
        }
    }
    if(ErrorNum) {
        par(las = 2)
        all.data.errors <- x$all.data.run$selec.history$OOB
        ngenes <- x$all.data.run$selec.history$Number.Variables
        maxplot <-
            max(
                c(unlist(lapply(x$allBootRuns,
                                function(x) max(x$selec.history$OOB))),
                  all.data.errors))
        minplot <-
            min(
                c(unlist(lapply(x$allBootRuns,
                                function(x) min(x$selec.history$OOB))),
                  all.data.errors))
        minplot <- minplot * (1 - 0.1)
        maxplot <- maxplot * (1 + 0.2)
        plot(ngenes, all.data.errors, 
             type = "l", axes = TRUE, xlab = "Number of variables",
             ylab = "OOB Error rate", ylim = c(minplot, maxplot),
             lty = 1,
             log = "x",
             col = "red", lwd = 2,
             main = "OOB Error rate vs. Number of variables in predictor",
             xlim = c(2, max(ngenes)*1.1))
###             plot(num.points.plot:1,
###                  x$allBootRuns[[1]]$other$trained.pam.cv$error,
###                  type = "l", axes = FALSE, xlab = "Number of genes",
###                  ylab = "CV Error rate", ylim = c(minplot, maxplot), lty = 2,
###                  main = "CV Error rate vs. Number of genes in predictor.")
        legend(x = 10, y = maxplot,
               legend = c("Bootstrap samples", "Original sample"),
               lty = c(2, 1), lwd = c(1, 3), col = c("Black", "Red"),
               bty = "n")
##        box()
##        axis(2)
##        axis(1)
        
        if(max(ngenes) > 300)
            axis(1, at = c(2, 3, 5, 8, 15, 20, 25, 35,
                    50, 75, 150, 200, 300),
                 labels = c(2, 3, 5, 8, 15, 10, 25, 35,
                 50, 75, 150, 200, 300))

        for(nb in 1:x$number.of.bootsamples)
            lines(x$allBootRuns[[nb]]$selec.history$Number.Variables,
                  x$allBootRuns[[nb]]$selec.history$OOB, lty = 2)
        lines(ngenes, all.data.errors, col = "red", lwd = 4)
    }
}



randomVarImpsRF <- function(xdata, Class, forest, numrandom = 100,
                            whichImp = "impsUnscaled",
                            usingCluster = TRUE,
                            TheCluster = NULL,
                            ...) {

  if(!all(whichImp %in% c("impsScaled", "impsUnscaled", "impsGini")))
      stop("whichImp contains a non-valid option; should be one or more \n",
           "of impsScaled, impsUnscaled, impsGini")

  
  ontree <- forest$ntree
  omtry <- forest$mtry
  
  nodesize <- 1
  
  if(usingCluster) {
      iRF2.cluster <- function(dummy, xdataTheCluster, ClassTheCluster,
                               ontree, omtry, nodesize, ...) {
        rf <- randomForest(x = xdataTheCluster,
                           y = sample(ClassTheCluster),
                           ntree = ontree,
                           mtry = omtry,
                           nodesize = nodesize,
                           importance = TRUE, keep.forest = FALSE,
                           ...)
        ## If we specify the importance measure, only that var is evaluated.
        ## if we say "ALL", all three.
        impsScaled <- NULL
        impsUnscaled <- NULL
        impsGini <- NULL
        if("impsUnscaled" %in% whichImp) 
          impsUnscaled <- importance(rf, type = 1, scale = FALSE)
        if("impsScaled" %in% whichImp)
          impsScaled <- importance(rf, type = 1, scale = TRUE)
        if("impsGini" %in% whichImp)
          impsGini <- rf$importance[, ncol(rf$importance)]
        
        return(list(impsScaled,
                    impsUnscaled,
                    impsGini))
    }
      
      outCl <- clusterApplyLB(TheCluster,
                              1:numrandom,
                              iRF2.cluster,
                              xdataTheCluster = xdata,
                              ClassTheCluster = Class,
                              ontree = ontree,
                              omtry = omtry,
                              nodesize = nodesize)


      
  } else {
      outCl <- list()
      cat("\n Obtaining random importances ")
      for(nriter in 1:numrandom) {
          cat(".")
        rf <- randomForest(x = xdata,
                           y = sample(Class),
                           ntree = ontree,
                           mtry = omtry,
                           nodesize = nodesize,
                           importance = TRUE,
                           keep.forest = FALSE,
                           ...)
        impsScaled <- NULL
        impsUnscaled <- NULL
        impsGini <- NULL
        if("impsUnscaled" %in% whichImp) 
          impsUnscaled <- importance(rf, type = 1, scale = FALSE)
        if("impsScaled" %in% whichImp)
          impsScaled <- importance(rf, type = 1, scale = TRUE)
        if("impsGini" %in% whichImp)
          impsGini <- rf$importance[, ncol(rf$importance)]
        outCl[[nriter]] <- list(impsScaled,
                                impsUnscaled,
                                impsGini)
      }
      cat("\n")
    } ##</else>
    
  randomVarImps <- list()
  
  if("impsScaled" %in% whichImp) {
    randomVarImps$impsScaled <-
        matrix(unlist(lapply(outCl, function(x) x[[1]])),
               ncol = length(outCl))
    colnames(randomVarImps$impsScaled) <- 1:numrandom
    rownames(randomVarImps$impsScaled) <- rownames(forest$importance)
  }
  if("impsUnscaled" %in% whichImp){
    randomVarImps$impsUnscaled <-
        matrix(unlist(lapply(outCl, function(x) x[[2]])),
               ncol = length(outCl))
    colnames(randomVarImps$impsUnscaled) <- 1:numrandom
    rownames(randomVarImps$impsUnscaled) <- rownames(forest$importance)
  }
  if("impsGini" %in% whichImp) {
    randomVarImps$impsGini <-
        matrix(unlist(lapply(outCl, function(x) x[[3]])),
               ncol = length(outCl))
    colnames(randomVarImps$impsGini) <- 1:numrandom
    rownames(randomVarImps$impsGini) <- rownames(forest$importance)
  }
  
  class(randomVarImps) <- c(class(randomVarImps),
                            "randomVarImpsRF")
  return(randomVarImps)
  ## Hopefully, we are returning a list with only the components selected.zz
}

randomVarImpsRFplot <- function(randomImportances,
                                forest,
                                whichImp = "impsUnscaled",
                                nvars = NULL,
                                show.var.names = FALSE,
                                vars.highlight = NULL,
                                main = NULL,
                                screeRandom = TRUE,
                                lwdBlack = 1.5,
                                lwdRed = 2,
                                lwdLightblue = 1,
                                cexPoint = 1,
                                overlayTrue = FALSE,
                                xlab = NULL,
                                ylab = NULL,
                                ...) {
    
    if(ncol(forest$importance) < 2)
        stop("The fitted rf", deparse(substitute(forest)),
             "was not fitted with importance = TRUE")
    
    randomImportances <-
        switch(whichImp,
               "impsUnscaled" = randomImportances$impsUnscaled,
               "impsScaled" = randomImportances$impsScaled,
               "impsGini" = randomImportances$impsGini,
               )

    originalForestImportance <-
        switch(whichImp,
               "impsUnscaled" =
               importance(forest, type = 1, scale = FALSE),
               "impsScaled" =
               importance(forest, type = 1, scale = TRUE),
               "impsGini" =
               forest$importance[, ncol(forest$importance)]
               )
    if(is.null(originalForestImportance))
        stop("\n Not valid 'whichImp' \n")
    
    if(is.null(xlab)) xlab <- "(Ordered) Variable"
    if(is.null(ylab)) ylab <-
        switch(whichImp,
               "impsUnscaled" = "Importance (unscaled)",
               "impsScaled"   = "Importance (scaled)",
               "impsGini"     = "Importance (Gini)",
               )
               
    nvars <- min(nvars, dim(randomImportances)[1])
    ylim <- range(originalForestImportance, randomImportances)    
    plottingOrder <- order(originalForestImportance,
                           decreasing = TRUE)[1:nvars]

    orderedOriginalImps <- originalForestImportance[plottingOrder]
    
    plot(orderedOriginalImps, type = "n", axes = FALSE, 
         xlab = xlab, ylab = ylab,
         main = main, ylim = ylim,
         ...)
    abline(h = 0, lty = 2, col = "blue")
    axis(2)
    box()
    if(show.var.names) {
        axis(1, labels = names(orderedOriginalImps),
             at = 1:length(orderedOriginalImps))
    } else {
        axis(1)
    }
    
    if(!overlayTrue)
         ###points(x = 1:nvars, orderedOriginalImps, lwd = lwdBlack,
           ###    col = "black", type = "b", cex = cexPoint)
        lines(x = 1:nvars, orderedOriginalImps, lwd = lwdBlack,
              col = "black", type = "b", cex = cexPoint)

    
    if(!is.null(vars.highlight)) {
        if(length(vars.highlight) > nvars) {
            warning("Not all vars. to highlight will be shown; increase nvars\n")
            cat("\n Not all vars. to highlight will be shown; increase nvars\n")
        }
        pos.selected <- which(names(orderedOriginalImps) %in% vars.highlight)
        if(!length(pos.selected)){
            warning("No selected vars. among those to show\n")
            cat("\nNo selected vars. among those to show\n")
        }
        else segments(pos.selected, 0,
                      pos.selected, orderedOriginalImps[pos.selected],
                      col = "blue", lwd = 2)
        if(length(pos.selected) < length(vars.highlight)) {
            warning("Not shown ",
                    length(vars.highlight) - length(pos.selected),
                    " of the 'to-highlight' variables\n")
            cat("Not shown ", length(vars.highlight) - length(pos.selected),
                    " of the 'to-highlight' variables\n")

        }
    }

    column.of.mean <- dim(randomImportances)[2] + 1
    if(screeRandom) {
        randomImportances <- apply(randomImportances, 2,
                                   sort, decreasing = TRUE)
        randomImportances <- randomImportances[1:nvars, ]
        randomImportances <- cbind(randomImportances,
                                   apply(randomImportances, 1, mean))
    } else {
        randomImportances <- randomImportances[plottingOrder, ]
        randomImportances <- cbind(randomImportances,
                                   apply(randomImportances, 1, mean))
    }

    matlines(randomImportances[, -column.of.mean],
             col = "lightblue", lwd = lwdLightblue, lty = 1)
    lines(x = 1:nvars,
          y = randomImportances[, column.of.mean],
          lwd = lwdRed, col = "red")

    if(overlayTrue) {
        points(x = 1:nvars, orderedOriginalImps, lwd = lwdBlack,
               col = "black", type = "b", cex = cexPoint)
        lines(x = 1:nvars, orderedOriginalImps, lwd = lwdBlack,
              col = "black", type = "l")
    }

}



varSelImpSpecRF <- function(forest,
                            xdata = NULL,
                            Class = NULL,
                            randomImps = NULL,
                            threshold = 0.10,
                            numrandom = 20,
                            whichImp = "impsUnscaled",
                            usingCluster = TRUE,
                            TheCluster = NULL,
                            ...) {

    if((is.null(xdata) | is.null(Class)) & is.null(randomImps))
        stop("You must specify a randomVarImpsRF object OR",
             "valid covariates and class objects.\n")

    if((!is.null(xdata) & !is.null(Class)) & !is.null(randomImps))
        warning("Using only the randomVarImpsRF object.",
             "Covariates and class objects discarded.\n", immediate. = TRUE)

    if(length(whichImp) > 1) stop("You can only use one importance measure")

    originalImps <- switch(whichImp,
                   "impsUnscaled" =
                   importance(forest, type = 1, scale = FALSE),
                   "impsScaled" =
                   importance(forest, type = 1, scale = TRUE),
                   "impsGini" =
                   forest$importance[, ncol(forest$importance)]
                   )
    if(is.null(originalImps))
        stop("\n Not valid 'whichImp' \n")

    originalImpsOrder <- order(originalImps,
                               decreasing = TRUE)

    if(is.null(randomImps)) {
        randomImps <-
            switch(whichImp,
                   "impsUnscaled" = randomVarImpsRF(xdata, Class, forest,
                   numrandom = numrandom,
                   whichImp = "impsUnscaled",
                   usingCluster = usingCluster,
                   TheCluster = TheCluster,
                   ...)$impsUnscaled, 
                   "impsUnscaled" = randomVarImpsRF(xdata, Class, forest,
                   numrandom = numrandom,
                   whichImp = "impsScaled",
                   usingCluster = usingCluster,
                   TheCluster = TheCluster,
                   ...)$impsScaled, 
                   "impsUnscaled" = randomVarImpsRF(xdata, Class, forest,
                   numrandom = numrandom,
                   whichImp = "impsGini",
                   usingCluster = usingCluster,
                   TheCluster = TheCluster,
                   ...)$impsGini, 
                   )               
    } else {
        elemento <- match(whichImp, names(randomImps))
        if(is.na(elemento))
            stop("The requested importance was not calculated\n",
                 "for the randomImps", deparse(substitute(randomImps)),"\n",
                 "object.\n")
        cat("\n Using the randomVarImpsRF", deparse(substitute(randomImps)),
            "object. xdata, Class, numrandom ignored.\n")
        randomImps <- randomImps[[elemento]]

    }
     
    randomImps <- apply(randomImps, 2,
                        sort, decreasing = TRUE)
    thresholds <- apply(randomImps, 1,
                        function(x) quantile(x, probs = 1 - threshold, type = 8))

    largest.value <- which(originalImps[originalImpsOrder] <= thresholds)[1]
    if(is.na(largest.value)) {
        selected.vars <- originalImpsOrder
        warning("All variables selected; could signal a problem")
    } else if(largest.value == 1) {
        selected.vars <- NA
        warning("No variables selected; could signal a problem")
    } else selected.vars <- originalImpsOrder[1:(largest.value - 1)]

    return(selected.vars)
}


selProbPlot <- function(object,
                        k = c(20, 100),
                        color = TRUE,
                        legend = FALSE,
                        xlegend = 68,
                        ylegend = 0.93,
                        cexlegend = 1.4,
                        main = NULL,
                        xlab = "Rank of gene",
                        ylab = "Selection probability",
                        pch = 19,
                        ...) {
    ## selection probability plots, such as in Pepe
    ## et al. 2003 (ROC paper).
    if(class(object) != "varSelRFBoot")
        stop("This function only works with objects created\n",
             "with the varSelRFBoot function.\n")
    nboot <- object$number.of.bootsamples
    if(nboot < 100)
        warning("You only used ", nboot,
                " bootstrap samples. Might be too few.",
                immediate. = TRUE)
    
    original.imps <- object$all.data.run$initialImportances
    original.ranks <- rank(-original.imps)
    boot.ranks <- lapply(object$allBootRuns,
                         function(x) {rank(-x$initialImportances)})
    boot.ranks <- matrix(unlist(boot.ranks), ncol = nboot)
    k1 <- apply(boot.ranks, 1, function(z) {sum(z <= k[1])/nboot})
    k2 <- apply(boot.ranks, 1, function(z) {sum(z <= k[2])/nboot})

    if(color) {
        plot(original.ranks[original.ranks <= k[2]],
             k2[original.ranks <= k[2]], xlim = c(1, k[2]),
             col = "red", pch = pch,
             xlab = xlab,
             ylab = ylab,
             main = main,
             ylim = c(0, 1),
             ...)
        points(original.ranks[original.ranks <= k[1]],
               k1[original.ranks <= k[1]],
               col = "blue", pch = pch)
        
        if(legend) legend(x = xlegend, y = ylegend,
                          legend = c(paste("Top", k[2]),
                          paste("Top", k[1])),
                          col = c("red", "blue"),
                          pch = pch, cex = cexlegend)
    } else {
        plot(original.ranks[original.ranks <= k[2]],
             k2[original.ranks <= k[2]], xlim = c(1, k[2]),
             pch = pch,
             xlab = "Rank of gene",
             ylab = "Selection probability",
             main = main, ylim = c(0, 1),
             ...)
        points(original.ranks[original.ranks <= k[1]],
               k1[original.ranks <= k[1]],
               pch = 21)
        if(legend) legend(x = xlegend, y = ylegend,
                          legend = c(paste("Top", k[2]),
                          paste("Top", k[1])),
                          pch = c(19, 21), cex = (cexlegend + 0.2))
    }
}





###################################################################
###################################################################
##########                                           ##############
##########       Custom code used in paper           ##############
##########       (unlikely to be of interest         ##############
##########       for anybody else)                   ##############
##########                                           ##############
###################################################################
###################################################################



    

figureSummary.varSelRFBoot <- function(object) {
## to create data for figures of paper with randomdata
    errorrate <- object$bootstrap.pred.error
    nvused <- length(object$all.data.vars)
    return(c(errorrate, nvused))
}


figureSummary2.varSelRFBoot <- function(object) {
## to create data for figures of paper with randomdata
    errorrate <- object$bootstrap.pred.error
    nvused <- length(object$all.data.vars)
    median.nvars <- median (object$number.of.vars)
    iq1.nvars <- quantile(object$number.of.vars, p = 0.25)
    iq3.nvars <- quantile(object$number.of.vars, p = 0.75)
    string1 <- paste(formatC(round(nvused, 0), width = 6)," (",
                     round(iq1.nvars, 0),", ",
                     round(median.nvars, 0), ", ",
                     round(iq3.nvars, 0),"),  ", 
                     formatC(round(errorrate, 3), width = 4)," \\\\ ",
                     sep = "")
    string1
}




tableSummary.varSelRFBoot <- function(object, name){
### to create data ready for the LaTeX tables of the paper

    neatname <- switch(name,
                        "gl.boot" = "Leukemia    ",
                        "vv.boot" = "Breast 2 cl.",
                        "vv3.boot" ="Breast 3 cl.",
                        "nci.boot" ="NCI 60      ",
                        "ra.boot" = "Adenocar.   ",
                     "brain.boot" = "Brain       ",
                     "colon.boot" = "Colon       ",
                  "lymphoma.boot" = "Lymphoma    ",
                  "prostate.boot" = "Prostate    ",
                     "srbct.boot" = "Srbct       ")
                        
    errorrate <- object$bootstrap.pred.error
    nvused <- length(object$all.data.vars)
    median.nvars <- median (object$number.of.vars)
    iq1.nvars <- quantile(object$number.of.vars, p = 0.25)
    iq3.nvars <- quantile(object$number.of.vars, p = 0.75)
    in.all.data <- which(names(table(object$all.vars.in.solutions)) %in% object$all.data.vars)
    tmp1 <- table(object$all.vars.in.solutions)[in.all.data]/object$number.of.bootsamples
    median.freq <- median(tmp1)
    if(nvused == 2) {
        iq1.freq <- min(tmp1)
        iq3.freq <- max(tmp1)
    } else {
        iq1.freq <- quantile(tmp1, p = 0.25)
        iq3.freq <- quantile(tmp1, p = 0.75)
    }

    if(nvused > 2) 
        string1 <-
            paste(neatname, "&     ",
                  formatC(round(errorrate, 3), width = 4)," &       ",
                  formatC(round(nvused, 0), width = 4)," &       ",
                  formatC(round(median.nvars, 0), width = 4)," (",
                  round(iq1.nvars, 0),", ",
                  round(iq3.nvars, 0),")   &   ",
                  formatC(round(median.freq, 2), width = 4)," (",
                  round(iq1.freq, 2),", ",
                  round(iq3.freq, 2),")\\\\",
                  sep = "") else   string1 <-
                      paste(neatname, "&     ",
                            formatC(round(errorrate, 3), width = 4)," &       ",
                            formatC(round(nvused, 0), width = 4)," &       ",
                            formatC(round(median.nvars, 0), width = 4)," (",
                            round(iq1.nvars, 0),", ",
                            round(iq3.nvars, 0),")   &   ",
                            formatC(round(median.freq, 2), width = 4)," (",
                            round(iq1.freq, 2),", ",
                            round(iq3.freq, 2),")\\footnotemark[1]\\\\",
                            sep = "")
    cat(string1, "\n")
}




###################################################################
###################################################################
##########                                           ##############
##########       Miscellaneous code                  ##############
##########                                           ##############
###################################################################
###################################################################



## rVI <- function(xdata, ydata, forest, numrandom = 40,
##                 whichImp = "impsUnscaled",
##                 ...) {
##     ## from randomVarImpsRF, but simplified to use only
##     ## one tytpe of importance
    
##     ontree <- forest$ntree
##     omtry <- forest$mtry

##     if(usingCluster) {
##         iRF.cluster <- function(dummy, ontree, omtry, ...) {
##             rf <- randomForest(x = xdataTheCluster,
##                                y = sample(ClassTheCluster),
##                                ntree = ontree,
##                                mtry = omtry,
##                                importance = TRUE, keep.forest = FALSE,
##                                ...)
##             imps <- switch(whichImp,
##                            "impsUnscaled" =
##                            importance(rf, type = 1, scale = FALSE),
##                            "impsScaled" =
##                            importance(rf, type = 1, scale = TRUE),
##                            "impsGini" =
##                            rf$importance[, ncol(rf$importance)]
##                            )
##             if(is.null(imps))
##                 stop("\n Not valid 'whichImp' \n")
##             return(imps)
##         }
##         clusterEvalQ(TheCluster,
##                      rm(list = c("xdataTheCluster", "ClassTheCluster")))
##         xdataTheCluster <<- xdata
##         ClassTheCluster <<- ydata
##         clusterExport(TheCluster,
##                       c("xdataTheCluster", "ClassTheCluster"))
##         outCl <- clusterApplyLB(TheCluster,
##                                 1:numrandom,
##                                 iRF.cluster,
##                                 ontree = ontree,
##                                 omtry = omtry)

##     } else {
##         outCl <- list()
##         for(nriter in 1:numrandom) {
##             rf <- randomForest(x = xdata,
##                                y = sample(ydata),
##                                ntree = ontree,
##                                mtry = omtry,
##                                importance = TRUE,
##                                keep.forest = FALSE,
##                                ...)
##             imps <- switch(whichImp,
##                            "impsUnscaled" =
##                            importance(rf, type = 1, scale = FALSE),
##                            "impsScaled" =
##                            importance(rf, type = 1, scale = TRUE),
##                            "impsGini" =
##                            rf$importance[, ncol(rf$importance)]
##                            )
##             if(is.null(imps))
##                 stop("\n Not valid 'whichImp' \n")
##             outCl[[nriter]] <- imps
##         }
##     }
    
##     randomVarImps <- matrix(unlist(outCl), ncol = numrandom)
##     rownames(randomVarImps) <- rownames(forest$importance)
##     colnames(randomVarImps) <- 1:numrandom
##     class(randomVarImps) <- c(class(randomVarImps),
##                               "rVI")
##     return(randomVarImps)
## }


##### old stuff

####nr <- 10; nc <- 20; x <- matrix(rnorm(2* nr*nc), ncol = nc)
####colnames(x) <- paste("v", 1:nc, sep ="")
####Class <- factor(c(rep("A", nr), rep("B", nr)))
####xdata <- x

####library(randomForest)
####rf1 <- randomForest(x = x, y = Class, importance = TRUE,
####                    keep.forest = FALSE, ntree = 2000)

####usingCluster <- TRUE
####if(usingCluster) {
####    library(snow)
####    library(Rmpi)
####    clusterNumberNodes <- 4
####    typeCluster <- "MPI"
####    TheCluster <- makeCluster(clusterNumberNodes,
####                              type = typeCluster)
####    clusterSetupSPRNG(TheCluster)
####    clusterEvalQ(TheCluster, library(randomForest))
####}



####screePlotRF <- function(forest, randomImportances,
####                        nvars = 50,
####                        show.var.names = FALSE,
####                        vars.highlight = NULL,
####                        main = NULL, scale = TRUE) {
####    nvars <- min(nvars, dim(randomImportances)[1])

    
####    if(scale) {
####        original.imps <- importance(forest, type = 1, scale = TRUE)
####    } else {
####        original.imps <- importance(forest, type = 1, scale = FALSE)
####    }
####    ordered.imps <- sort(original.imps, decreasing = TRUE)[1:nvars]

####    plot(ordered.imps, type = "b", axes = FALSE, lwd = 1.5,
####         xlab = "Variable", ylab = "Importance",
####         main = main)
####    abline(h = 0, lty = 2, col = "blue")
####    axis(2)
####    box()
####    if(show.var.names) {
####        axis(1, labels = names(ordered.imps),
####             at = 1:length(ordered.imps))
####    } else {
####        axis(1)
####    }

####    if(!is.null(vars.highlight)) {
####        if(length(vars.highlight) > nvars) {
####            warning("Not all vars. to highlight will be shown; increase nvars\n")
####            cat("\n Not all vars. to highlight will be shown; increase nvars\n")
####        }
####        pos.selected <- which(names(ordered.imps) %in% vars.highlight)
####        if(!length(pos.selected)){
####            warning("No selected vars. among those to show\n")
####            cat("\nNo selected vars. among those to show\n")
####        }
####        else segments(pos.selected, 0,
####                      pos.selected, ordered.imps[pos.selected],
####                      col = "blue", lwd = 2)
####        if(length(pos.selected) < length(vars.highlight)) {
####            warning("Not shown ", length(vars.highlight) - length(pos.selected),
####                    " of the 'to-highlight' variables\n")
####            cat("Not shown ", length(vars.highlight) - length(pos.selected),
####                    " of the 'to-highlight' variables\n")

####        }
####    }

####    randomImportances <- randomImportances[1:nvars, ]
####    column.mean <- dim(randomImportances)[2]
####    matlines(randomImportances[, -column.mean],
####             col = "green", lty = 1)
####    lines(x = 1:nvars,
####          y = randomImportances[, column.mean],
####          lwd = 1.3, col = "red")
####}


####screePlotRF(rf1, ii1, nvars = 20, vars.highligh = paste("v", c(8, 4), sep =""), show.var.names = TRUE)

####i1 <- randomVarImpsRF(xdata, Class, rf1, numrandom = 3)
####screePlotRF(rf1, i1)






### This version allows optimizing mtry and ntree; leave code here,
### but don't use.

## FIXME: uncomment, but fix global function uses. And add examples.

## ExperimentalvarSelRF <- function(xdata, Class, vars.drop.num = NULL,
##                      vars.drop.frac = 0.5,
##                      c.sd = 1,
##                      whole.range = TRUE,
##                      recompute.var.imp = FALSE,
##                      verbose = FALSE,
##                      returnFirstForest = TRUE,
##                      ## next are all for tune2RF
##                      ## but this is a mess; should pass them as
##                      ## control.
##                      tuneMtry = FALSE, tuneNtree = FALSE,
##                      startNtree = 1000, startMtryFactor = 1,
##                      stepFactorMtry = 1.25,
##                      stepFactorNtree = 1.75,
##                      minCorNtree = 0.975,
##                      quantNtree = 0.025,
##                      returnForest = TRUE,
##                      ntreeTry = 2000) {

##     if( (is.null(vars.drop.num) & is.null(vars.drop.frac)) |
##        (!is.null(vars.drop.num) & !is.null(vars.drop.frac)))
##         stop("One (and only one) of vars.drop.frac and vars.drop.num must be NULL and the other set")
    

##     max.num.steps <- dim(xdata)[2]
##     num.subjects <- dim(xdata)[1]

##     if(is.null(colnames(xdata)))
##         colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep ="")
    
##     ##oversize the vectors; will prune later.
##     n.vars <- vars <- OOB.rf <- OOB.sd <- rep(NA, max.num.steps)
    
##     ## First get optimal values for mtry and ntree, and get first run:
    
##     rf <- tune2RF(x = xdata, y = Class,
##                   tuneMtry = tuneMtry, tuneNtree = tuneNtree,
##                   startNtree = startNtree, mtryFactor = mtryFactor,
##                   stepFactorMtry = stepFactorMtry,
##                   stepFactorNtree = stepFactorNtree,
##                   minCorNtree = minCorNtree,
##                   quantNtree = quantNtree,
##                   returnForest = TRUE,
##                   ntreeTry = ntreeTry,
##                   verbose = verbose)

##     ## We'll need it for future runs
##     ntree <- rf$ntree
##     mtry <- rf$mtry

##     if(returnFirstForest)
##         FirstForest <- rf
##     else
##         FirstForest <- NULL
    
## #    rf <- randomForest(x = xdata, y = Class, importance= TRUE,
## #                       ntree = ntree, keep.forest = FALSE)
##     m.iterated.ob.error <- m.initial.ob.error <- oobError(rf)
##     sd.iterated.ob.error <- sd.initial.ob.error <-
##         sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) * (1/num.subjects))

##     if(verbose) {
##         print(paste("Initial OOB error: mean = ", round(m.initial.ob.error, 4),
##                     "; sd = ", round(sd.initial.ob.error, 4), sep = ""))
##     }

## ### previous code (before rF < 4.3.1)
##     #selected.vars <- order(rf$importance[, (ncol(rf$importance) - 1)], decreasing = TRUE)
## #    ordered.importance <- rf$importance[selected.vars,(ncol(rf$importance) -1)]

##     importances <- importance(rf, type = 1, scale = FALSE)
##     selected.vars <- order(importances, decreasing = TRUE)
##     ordered.importances <- importances[selected.vars]
    
##     initialImportances <- importances
##     initialOrderedImportances <- ordered.importances
    
##     j <- 1
##     n.vars[j] <- dim(xdata)[2] 
##     vars[j] <- paste(colnames(xdata), collapse = " + ")
##     OOB.rf[j] <- m.iterated.ob.error
##     OOB.sd[j] <- sd.iterated.ob.error

##     var.simplify <- TRUE
    
##     while(var.simplify) {
##         last.rf <- rf
##         last.vars <- selected.vars
## #        print(paste(".........Number of variables before selection",
## #                    dim(xdata)[2])) ## debug
##         previous.m.error <- m.iterated.ob.error
##         previous.sd.error <- sd.iterated.ob.error

##         if(recompute.var.imp & (j > 1)) {
##             ## need to set indexes as absolute w.r.t. original data
## ####            tmp.order <- order(rf$importance[, (ncol(rf$importance) - 1)],
## ####                                    decreasing = TRUE)
## ####            selected.vars <- selected.vars[tmp.order]
## ####            ordered.importance <- rf$importance[tmp.order, (ncol(rf$importance) -1)]

##             importances <- importance(rf, type = 1, scale = FALSE)
##             tmp.order <- order(importances, decreasing = TRUE)
##             selected.vars <- selected.vars[tmp.order]
##             ordered.importances <- importances[tmp.order]

##         }
        
##         num.vars <- length(selected.vars)
  
## ####        if(any(is.na(ordered.importances))) {
## ####            print("**********  Nas in ordered.importances ******")
## ####            browser()
## ####        }
           
##         if(any(ordered.importances < 0)) {
##             selected.vars <- selected.vars[-which(ordered.importances < 0)]
##             ordered.importances <- ordered.importances[-which(ordered.importances < 0)]
##         } else {
##             if(is.null(vars.drop.num))
##                 vars.drop <- round(num.vars * vars.drop.frac)
##             else vars.drop <- vars.drop.num
                
##             if(num.vars >= (vars.drop + 2)) {
##                 selected.vars <- selected.vars[1: (num.vars - vars.drop)]
##                 ordered.importances <- ordered.importances[1: (num.vars - vars.drop)]
##             }
##             else {
##                 selected.vars <- selected.vars[1:2]
##                 ordered.importances <- ordered.importances[1:2]
##             }
##         }
##         ## couldn't we eliminate the following?
##         if((length(selected.vars) < 2) | (any(selected.vars < 1))) {
##             var.simplify <- FALSE
##             break
##         }

##         if(length(selected.vars) <= 2) var.simplify <- FALSE
      
##         if(recompute.var.imp) 
##             rf <- randomForest(x = xdata[, selected.vars], y = Class, importance= TRUE,
##                                ntree = ntree, mtry = mtry, keep.forest = FALSE)
##         else
##             rf <- randomForest(x = xdata[, selected.vars], y = Class, importance= FALSE,
##                                ntree = ntree, mtry = mtry, keep.forest = FALSE)
        
##         m.iterated.ob.error <- oobError(rf)
##         sd.iterated.ob.error <-
##             sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) * (1/num.subjects))
        
##         if(verbose) {
##             print(paste("..... iteration ", j, "; OOB error: mean = ",
##                         round(m.iterated.ob.error, 4),
##                         "; sd = ", round(sd.iterated.ob.error, 4), sep = ""))
##         }
##         j <- j + 1

        
##         n.vars[j] <- length(selected.vars)
##         vars[j] <- paste(colnames(xdata)[selected.vars],
##                          collapse = " + ")
##         OOB.rf[j] <- m.iterated.ob.error
##         OOB.sd[j] <- sd.iterated.ob.error


##         if(!whole.range &
##            (
##             (m.iterated.ob.error >
##              (m.initial.ob.error + c.sd*sd.initial.ob.error))
##             |
##             (m.iterated.ob.error >
##              (previous.m.error + c.sd*previous.sd.error)))
##            )
##             var.simplify <- FALSE
##     }

##     if (!whole.range) {
##         if(!is.null(colnames(xdata)))
##             selected.vars <- sort(colnames(xdata)[last.vars])
##         else
##             selected.vars <- last.vars

##         out <- list(selec.history = data.frame(
##                     Number.Variables = n.vars,
##                     Vars.in.Forest = vars,
##                     OOB = OOB.rf,
##                     sd.OOB = OOB.sd)[1:j,],
##                     rf.model = last.rf,
##                     selected.vars = selected.vars,
##                     selected.model =  paste(selected.vars, collapse = " + "),
##                     best.model.nvars = length(selected.vars),
##                     initialImportances = initialImportances,
##                     initialOrderedImportances = initialOrderedImportances,
##                     ntree = ntree,
##                     mtry = mtry,
##                     firstForest = FirstForest)
##         class(out) <- "varSelRF"
##         return(out)
##     }
##     else {
##         n.vars <- n.vars[1:j]
##         vars <- vars[1:j]
##         OOB.rf<- OOB.rf[1:j]
##         OOB.sd <- OOB.sd[1:j]
##         ##browser()
##         min.oob.ci <- min(OOB.rf) + c.sd * OOB.sd[which.min(OOB.rf)]
##         best.pos <-
##             which(OOB.rf <= min.oob.ci)[which.min(n.vars[which(OOB.rf <= min.oob.ci)])]

##         selected.vars <- sort(unlist(strsplit(vars[best.pos],
##                                               " + ", fixed = TRUE)))
##         out <- list(selec.history = data.frame(
##                     Number.Variables = n.vars,
##                     Vars.in.Forest = vars,
##                     OOB = OOB.rf,
##                     sd.OOB = OOB.sd),
##                     rf.model = NA,
##                     selected.vars = selected.vars,
##                     selected.model = paste(selected.vars, collapse = " + "),
##                     best.model.nvars = n.vars[best.pos],
##                     initialImportances = initialImportances,
##                     initialOrderedImportances = initialOrderedImportances,
##                     ntree = ntree,
##                     mtry = mtry,
##                     firstForest = FirstForest)
##         class(out) <- "varSelRF"
##         return(out)
##     }
## }



## FIXME: uncomment, but fix global function uses. And add examples.

## tune2RF <- function(x, y, tuneMtry = TRUE, tuneNtree = TRUE,
##                     startNtree = 1000, startMtryFactor = 1,
##                     stepFactorMtry = 1.25,
##                     stepFactorNtree = 1.75,
##                     minCorNtree = 0.9,
##                     quantNtree = 0.025, ## with 4000, look at first 100
##                     returnForest = TRUE,
##                     ntreeTry = 2000,
##                     verbose = TRUE,
##                     plot = FALSE,
##                     ...) {
##     ## if tuneMtry = TRUE and tuneNtree = TRUE,
##     ##        ntree is used for tuneRF, mtryFactor is ignored
##     ## if tuneMtry = TRUE and tuneNtree = FALSE
##     ##        ntree used for tuneRF and final ntree, and mtry factor is ignored
##     ## if tuneMtry = FALSE and tuneRF = TRUE,
##     ##        ntree is used only as starting value for search of ntree and
##     ##        mtryFactor is used for mtry
##     ## if tuneMtry = FALSE and tuneRF = FALSE,
##     ##        mtryFactor and ntree are the ntree and mtry used.

##     if(plot) {
##         op <- par(mfrow = c(1,2))
##         on.exit(par(op))
##     }
##     if(tuneMtry) {
##         tunedMtry <- tuneRF(x, y, stepFactor = stepFactorMtry,
##                             ntreeTry = ntreeTry,
##                             mtryStart = floor(sqrt(ncol(x)) * startMtryFactor),
##                             doBest = FALSE,
##                             plot = plot,
##                             trace = verbose)
##         tunedMtry <-
##             tunedMtry[which.min(tunedMtry[, 2]), 1]
##     } else {
##         tunedMtry <- floor(sqrt(ncol(x)) * startMtryFactor)
##     }

##     tunedNtree <- startNtree
##     if(tuneNtree) {
##         if(verbose)
##             cat("\n Tunning ntree: initial forest construction \n")
##         f1 <- randomForest(x, y, mtry = tunedMtry,
##                            importance = TRUE,
##                            keep.forest = FALSE,
##                            ntree = tunedNtree)
##         repeat {
##             f2 <- randomForest(x, y, mtry = tunedMtry,
##                                importance = TRUE,
##                                keep.forest = FALSE,
##                                ntree = round(stepFactorNtree * tunedNtree))
##             m1 <- cbind(
##                         importance(f1, type = 1, scale = FALSE),
##                         importance(f2, type = 1, scale = FALSE))
##             m1[m1 < quantile(m1, 1-quantNtree)] <- NA
##             m1 <- na.omit(m1)
##             m2 <- m1
##             m2[m2 <= 0] <- NA
##             m2 <- na.omit(m2)
##             mc1 <- cor(m1)
##             mc.rob <- cov.rob(m2, method = "mcd", cor = TRUE)$cor
##             if(verbose) {
##                 cat("\n ... tunning ntree; cor. importances successive ntrees\n")
##                 colnames(mc1) <- c(tunedNtree, round(stepFactorNtree * tunedNtree))
##                 rownames(mc1) <- c(tunedNtree, round(stepFactorNtree * tunedNtree))
##                 colnames(mc.rob) <- c(tunedNtree, round(stepFactorNtree * tunedNtree))
##                 rownames(mc.rob) <- c(tunedNtree, round(stepFactorNtree * tunedNtree))
##                 cat("\n         correlation matrix\n")
##                 print(round(mc1, 4))
##                 cat("\n         robust correlation matrix\n")
##                 print(round(mc.rob, 4))
##                 cat(" \n using ", dim(m2)[1], "observations\n")
##                 cat("\n")
##             }
##             if((min(mc1) > minCorNtree) &
##                (min(mc.rob) > minCorNtree)) {
##                 plot(m1[, 1], m1[, 2], xlab = paste("Ntree", tunedNtree),
##                      ylab = paste("Ntree", round(stepFactorNtree * tunedNtree)),
##                      main =
##                      paste("Upper ", quantNtree,
##                            "th quantile of importances", sep = ""))
                     
##                 break
##             }
##             tunedNtree <- round(stepFactorNtree * tunedNtree)
##             f1 <- f2
##         }
##     }
##     if(returnForest) {
##         if(tuneNtree)
##             return(f1)
##         else
##             return(randomForest(x, y, mtry = tunedMtry,
##                                 importance = TRUE,
##                                 keep.forest = FALSE,
##                                 ntree = ntree)
##                    )
##     }
##     else {
##         return(c(tunedMtry = tunedMtry, tunedNtree = tunedNtree))
##     }
## }




boot.imp <- function(data, class, ntree = 20000, B = 200) {
## just checking the output from the pg.plots. is correct.
    N <- length(class)
    mat.out <- matrix(NA, nrow = dim(data)[2], ncol = B)
    for(i in 1:B) {
        sample.again <- TRUE
        while(sample.again) {
            bootsample <- unlist(tapply(1:N,  class,
                                        function(x)
                                        sample(x, size = length(x),
                                               replace = TRUE)))
            ## sure, this isn't the fastest, but will do for now.
            nobootsample <- setdiff(1:N, bootsample)
            if(!length(nobootsample))
                sample.again <- TRUE
            else sample.again <- FALSE
        }
        ## this is an ugly hack to prevetn nobootsamples
        ## of size 0.
        train.data <- data[bootsample, , drop = FALSE]
        train.class <- class[bootsample]
        rftmp <- randomForest(x = train.data, y = train.class,
                              ntree = ntree,
                              keep.forest = FALSE,
                              importance = TRUE)
        mat.out[, i] <- importance(rftmp, type = 1, scale = FALSE)
    }
    return(mat.out)
}


####gl.b.i <- boot.imp(gl.data, gl.class, B = 20)
        
####pairs(cbind(gl.b.i, importance(gl.20000.rf, type = 1, scale = FALSE)),
####      pch = ".")

####pairs(cbind(gl.b.i[, 1:10], importance(gl.20000.rf, type = 1, scale = FALSE)),
####      pch = ".")

####summary(cbind(gl.b.i, importance(gl.20000.rf, type = 1, scale = FALSE)))
####boxplot(data.frame(cbind(gl.b.i, importance(gl.20000.rf, type = 1, scale = FALSE))))




#####f.mtr <- function(x, mf = 3) {
#####  x2 <- floor(sqrt(x) * mf)
#####  x2[x2 > x] <- x[x2 > x]
#####  return(x2)
#####}
#####plot(x = sqrt(2:200), f.mtr(2:200, mf = 5), type = "l")
