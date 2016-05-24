svmlight<-function (x, ...) 
    UseMethod("svmlight")


svmlight.default <- function(x, grouping, temp.dir=NULL, pathsvm=NULL, del=TRUE, 
    type="C", class.type = "oaa" ,svm.options=NULL,  prior=NULL, out=FALSE, ...)
{
y <- grouping
 pick<-NULL     
 if (type!="R")
   {
    ### Construct Dummymatrix for 1-a-a classification: 
    ### ncol= no. of classes, (i,j)=+1 iff obj i comes from class j, -1 else 
    if(is.null(prior)) prior <- table(y) / length(y)
    ys <- as.factor(y)
    tys <- table(ys)
    lev <- levels(ys)
    if (class.type !="oao")
      {
       class.type<-"oaa"
       ymat <- matrix(-1, nrow = nrow(x), ncol = length(tys))
       ymat[cbind(seq(along = ys), sapply(ys, function(x) which(x == lev)))] <- 1
      }
     else
       { ## Classification: one against one
        nclass <- length(table(levels(y)))
        m <- (nclass - 1)
        minus <- nclass + 1 - sequence(m:1)
        plus <- rep(1:m, m:1)
        pick <- rbind(plus, minus)
        xsplit <- split(data.frame(x), ys)
        ymat <- list()
        xlist <- list()
        for(k in 1:ncol(pick)){
            ymat[[k]] <- c(rep(1, nrow(xsplit[[ pick[1, k] ]])), rep(-1, nrow(xsplit[[ pick[2, k] ]])))
            xlist[[k]] <- rbind(xsplit[[ pick[1, k] ]], xsplit[[ pick[2, k] ]])
            }
        } 
    counts <- as.vector(tys)     
  }
 else
   {
   ### Regression: So transform vector to 1-col matrix and set options, J ...
     ymat <- matrix(y,ncol=1)
     lev <- NULL
     counts <- 1
     svm.options <- paste("-z r ",svm.options)
     J <- 1 
     }
    svm.model <- list()

### call svm_learn
    cmd <- if (is.null(pathsvm)) 
        "svm_learn"
    else file.path(pathsvm, "svm_learn")
    
    if(is.matrix(ymat)) J <- 1:ncol(ymat)
    if(is.list(ymat)) J <- 1:length(ymat)
    
### "paste" file names for training data and model
    train.filename <- paste(temp.dir, "_train_", J, ".dat", sep = "")
    model.filename <- paste(temp.dir, "_model_", J, ".txt", sep = "")
    PWin <- .Platform$OS.type == "windows"
 for (j in J){
      ### construct matrix to learn svm in a format, so that svmlight can read it
      if (class.type !="oao")
         {
          train <- svmlight.file(cbind(ymat[,j], x), train = TRUE)    
          }
        else
          {
           train <- svmlight.file(cbind(ymat[[j]], xlist[[j]]), train = TRUE)
           }  
      ### save to disk
          write.table(train, file = train.filename[j], row.names = FALSE, 
          col.names = FALSE, quote = FALSE)
          if (PWin) 
              system(paste(cmd, svm.options, train.filename[j], model.filename[j]), 
                  show.output.on.console = out)
          else 
              system(paste(cmd,svm.options, train.filename[j], model.filename[j]))
      ### store learned model
          svm.model[[j]] <- readLines(model.filename[j])
      }
    if (del) 
       file.remove(c(train.filename, model.filename))

    cl <- match.call()
    cl[[1]] <- as.name("svmlight")
    structure(list(prior = prior, counts = counts, lev = lev, 
        temp.dir = temp.dir, pathsvm = pathsvm, del = del, 
        type=type, class.type=class.type, pick=pick, J=J,
        svm.model = svm.model, svm.options = svm.options, call = cl), class = "svmlight")
}


#### svmlight interface for different calls. Copied and adapted from lda.xxx
svmlight.formula <- function(formula, data = NULL, ..., subset, na.action = na.fail) 
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
    res <- svmlight.default(x, grouping, ...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1]] <- as.name("svmlight")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        res$na.action <- attr(m, "na.action")
    res
}

svmlight.matrix <- function(x, grouping, ..., subset, na.action = na.fail) 
{
    if (!missing(subset)) {
        x <- x[subset, , drop = FALSE]
        grouping <- grouping[subset]
    }
    if (!missing(na.action)) {
        dfr <- na.action(structure(list(g = grouping, x = x), 
            class = "data.frame"))
        grouping <- dfr$g
        x <- dfr$x
    }
    res <- svmlight.default(x, grouping, ...)
    cl <- match.call()
    cl[[1]] <- as.name("svmlight")
    res$call <- cl
    res
}

svmlight.data.frame <- function (x, ...) 
{
   res <- svmlight.matrix(structure(data.matrix(x), class = "matrix"), 
        ...)
    cl <- match.call()
    cl[[1]] <- as.name("svmlight")
    res$call <- cl
    res
}


### predict method for svmlight
predict.svmlight <- function(object, newdata, scal=TRUE, ...)
{
### utility function for oao classification
sf<-function(x,pick)
{
erg<-numeric(max(pick))
names(erg)<-1:max(pick)
dummy<-table(diag(pick[(x<0)+1,]))
erg[names(dummy)]<-dummy
return(erg)
}




### copied and adapted form predict.lda
    if (!inherits(object, "svmlight")) 
        stop("object not of class", "'svmlight'")
    if (!is.null(Terms <- object$terms)) {
        if (missing(newdata)) 
            newdata <- model.frame(object)
        else {
            newdata <- model.frame(as.formula(delete.response(Terms)), 
                newdata, na.action = function(x) x, xlev = object$xlevels)
        }
        x <- model.matrix(delete.response(Terms), newdata, contrasts = object$contrasts)
        xint <- match("(Intercept)", colnames(x), nomatch = 0)
        if (xint > 0) 
            x <- x[, -xint, drop = FALSE]
    }
    else {
        if (missing(newdata)) {
            if (!is.null(sub <- object$call$subset)) 
                newdata <- eval.parent(parse(text = paste(deparse(object$call$x, 
                  backtick = TRUE), "[", deparse(sub, backtick = TRUE), 
                  ",]")))
            else newdata <- eval.parent(object$call$x)
            if (!is.null(nas <- object$call$na.action)) 
                newdata <- eval(call(nas, newdata))
        }
        if (is.null(dim(newdata))) 
            dim(newdata) <- c(1, length(newdata))
        x <- as.matrix(newdata)
    }
#######################################

### save test data on disk
    x <- svmlight.file(cbind(rep(0,nrow(x)),x), train = TRUE)
    #x <- svmlight.file(x, train = FALSE)
    test.filename <- paste(object$temp.dir, "_test_.dat", sep = "")
    write.table(x, file = test.filename, row.names = FALSE, 
        col.names = FALSE, quote = FALSE)
    werte <- NULL
    
    cmd <- if (is.null(object$pathsvm)) 
        "svm_classify"
    else file.path(object$pathsvm, "svm_classify")
    J <- object$J
    model.filename <- paste(object$temp.dir, "_model_", J, ".txt", sep = "")
    pred.filename <- paste(object$temp.dir, "_pred_", J, ".txt", sep = "")
### for all learned models predict call svm_classify
    for (j in J){
        writeLines(object$svm.model[[j]], model.filename[j])
        system(paste(cmd, test.filename, model.filename[j], pred.filename[j]))
    ### read predicted values
        prognose <- read.table(pred.filename[j], header = FALSE)[ , 1]
        werte <- cbind(werte, prognose)
    }
    if (object$del) 
        file.remove(c(test.filename, pred.filename, model.filename))
    if (object$type=="C")
    {
    ### Classification: choose class with highest decision value f(x)
    if  (object$class.type=="oao")   
        {
        werte2<-apply(werte,1,sf,pick=object$pick) 
        werte<-t(werte2)
        }
    
    
    
    classes <- factor(max.col(werte), levels = seq(along = object$lev), 
        labels = object$lev)
    colnames(werte)<-object$lev 
    
    if (scal) werte<-e.scal(werte)$sv
    return(list(class = classes, posterior = werte))
  }
   else {return(as.vector(werte))} ### Regression
}


svmlight.file <- function(x, train = FALSE, ...)
{
    if(is.vector(x)) x <- t(x)
    erg <- x
    sn <- 1:nrow(x) 
    if(!train) erg[sn, 1] <- paste("1:", x[sn, 1], sep = "")
    if(ncol(x) > 1){
        j <- 2:ncol(x)
        erg[ , -1] <- matrix(paste(j - train, t(x[,j]), sep = ":"), ncol = ncol(x)-1, byrow = TRUE)
    }
    return(erg)
} 
