meclight <- function(x, ...) 
    UseMethod("meclight")
##################################################################################
##################################################################################


# meclight: INPUT

# x: matrix containing the data, values of one observation = one row
# grouping: vector containing group memberships, matching to X
# d: wanted dimension after projetion
# number: number of cross-validations for estimating the misclassification rate, default = 10


# meclight: OUTPUT

# G.optimal: estimated projection matrix G
# error.rate: misclassification rate (per cent) using G.optimal
# improvement: improvement (per cent) using G.optimal compared to using LDA


meclight.default <- function(x, grouping, r = 1, fold = 10, ...){
    
    
# error: INPUT

# G : starting matrix, e.g. from LDA, G has to be put in as an vector (column by column!)
# x: matrix containing the data, values of one observation = one row
# grouping: vector containing group memberships, matching to x
# d: wanted dimension after projetion
# fold: fold of cross-validations for estimating the misclassification rate


# error: OUTPUT

# wrong.percent: estimated misclassification rate (per cent) using LDA


    error <- function(G, r, d.sample, x, grouping, fold, ...){         
        
        G <- matrix(G, ncol = r)
        # splitting up the data
        d <- length(d.sample)
        length.iter <- floor(d/fold)
        
        # projecting the data into the lower dimensional space
        x.low <- t(G) %*% t(x)
        wrong <- numeric(fold)
        
        if(is.null(dim(G))) G <- as.matrix(G)
        G.dimension <- dim(G)[2]
    
            for (i in 1:fold){
                
                # test set
                d.aside <- d.sample[1:length.iter]         
                if (G.dimension == 1){
                    data.training <- data.frame(cbind(x.low[ , -d.aside], 
                        grouping = as.factor(grouping[-d.aside])))
                    data.test <- data.frame(cbind(V1 = x.low[,d.aside], 
                        grouping = as.factor(grouping[d.aside])))
                }
                else{
                    # the data has to be transposed if dimension > 2
                    data.training <- data.frame(cbind(t(x.low[ , -d.aside]), 
                        grouping = as.factor(grouping[-d.aside])))
                    data.test <- data.frame(cbind(V1 = t(x.low[ , d.aside]), 
                        grouping = as.factor(grouping[d.aside])))
                }
            
                # renew LDA using training set
                z.new <- lda(grouping ~ ., data = data.training, ...)
                
                # estimate group membership of test set
                pr <- predict(z.new, data.test)$class
            
                # fold of wrong results
                right <- sum(pr == data.test$grouping)
                wrong[i] <- length.iter - right
                
                d.sample <- d.sample[-(1:length.iter)]
            }
        
        # error rate:
        wrong.percent <- sum(wrong) / (fold*length.iter)
        
    return(wrong.percent)
    }
    
    
    # creating starting matrix by using LDA
    data.lda <- data.frame(x, grouping = grouping)
    n <- nrow(x)
    if(n < fold){
        warning("not enough observations for ", fold, "-fold, instead setting fold = ", n)
        fold <- n
    }
    lda_obj <- lda(grouping ~ ., data.lda, ...)$scaling
    n_lda <- ncol(lda_obj)
    if(r > n_lda){
        warning("dimension reduced to r = min(n.levels - 1, n.variables) = ", n_lda)
        r <- n_lda
    }

    G.start <- as.vector(lda_obj[ , 1:r])
    daten <- data.matrix(cbind(x, grouping))
    
    d.sample <- sample(dim(x)[1])
    error.start <- error(G.start, r = r, d.sample = d.sample, x = x, 
        grouping = grouping, fold = fold, ...)
    
    # optimizing the misclassification rate using Nelder-Mead
    result <- optim(G.start, error, r = r, d.sample = d.sample, x = x, 
        grouping = grouping, fold = fold, ...)
    
    G.result <- matrix(result$par, ncol = r)
    G.error.rate <- result$value
    improvement <- error.start - G.error.rate
    
    rownames(G.result) <- colnames(x)
    colnames(G.result) <- paste("Proj.dim", 1:ncol(G.result))
    projd <- daten[ , 1:(ncol(daten)-1)] %*% G.result
    colnames(projd) <- paste("Proj.dim", 1:ncol(projd))
    daten.proj <- data.frame(grouping, projd)
    model <- lda(grouping ~ ., data = daten.proj, ...)

    cl <- match.call()
    cl[[1]] <- as.name("meclight")
    structure(list(model = model, B.error = G.error.rate, 
        B.impro = improvement, Proj.matrix = G.result, Proj.dim = r,
        call = cl, newx = daten.proj), class = "meclight")
}

meclight.formula <- function (formula, data = NULL, ..., subset, na.action = na.fail) 
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
    res <- meclight.default(x, grouping, ...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1]] <- as.name("meclight")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        res$na.action <- attr(m, "na.action")
    res
}

meclight.matrix <- function (x, grouping, ..., subset, na.action = na.fail) 
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
    res <- NextMethod("meclight")
    cl <- match.call()
    cl[[1]] <- as.name("meclight")
    res$call <- cl
    res
}

meclight.data.frame <- function(x, ...) 
{
    res <- meclight(structure(data.matrix(x), class = "matrix"), ...)
    cl <- match.call()
    cl[[1]] <- as.name("meclight")
    res$call <- cl
    res
}  

###########################################################################################################



predict.meclight <- function(object, newdata, ...){
    if(missing(newdata)) 
        return(predict(object$model, ...))
    else 
        if(length(grep("~", as.character(object$call))) == 0) 
            projnewdata <- data.matrix(newdata) %*% object$Proj.matrix
        else projnewdata <-
            data.matrix(newdata[ , attr(terms(object), "term.label")]) %*% 
                object$Proj.matrix
    colnames(projnewdata) <- paste("Proj.dim", 1:ncol(projnewdata))
    prediction <- predict(object$model, data.frame(projnewdata), ...)
    return(prediction)   
}

print.meclight <- function(x,  ...){
  cat("Dimension of projection:  ", x$Proj.dim, "\n")
  cat("est. bootstrap error rate:", x$B.error,  "\n")
  cat("est. improvement to LDA:  ", x$B.impro,  "\n\n")
  cat("Projection matrix:\n")
  print(x$Proj.matrix, ...)
  invisible(x)
}
