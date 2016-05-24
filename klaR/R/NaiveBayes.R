NaiveBayes <- function (x, ...) 
    UseMethod("NaiveBayes")

NaiveBayes.formula <- function (formula, data, ..., subset, na.action = na.pass) 
{
    call <- match.call()
    Yname <- as.character(formula[[2]])
    if (is.data.frame(data)) {
        m <- match.call(expand.dots = FALSE)
        m$... <- NULL
        m$na.action <- na.action
        m[[1]] <- as.name("model.frame")
        m <- eval(m, parent.frame())
        Terms <- attr(m, "terms")
        if (any(attr(Terms, "order") > 1)) 
            stop("NaiveBayes cannot handle interaction terms")
        Y <- model.extract(m, "response")
        X <- m[ , -attr(Terms, "response"), drop=FALSE]
        return(NaiveBayes(X, Y, ...))
    }
    else stop("NaiveBayes formula interface handles data frames only")
}


NaiveBayes.default <- function (x, grouping, prior = NULL, usekernel = FALSE, fL = 0, ...) 
{
    x <- data.frame(x)
    if(!is.factor(grouping))
        stop("grouping/classes object must be a factor")
    if (is.null(prior)) 
        apriori <- table(grouping) / length(grouping)
    else 
        apriori <- as.table(prior / sum(prior))
    call <- match.call()
    Yname <- "grouping"
    LaplaceEst <- function(x, f = 0)  
        t(apply(x, 1, function(u) (u + f)/(sum(u) + (length(u) * f))))
        
    est <- function(var){
        if(is.numeric(var)) {
            temp <- if (usekernel)
                lapply(split(var, grouping), FUN = function(xx) density(xx, ...))
            else
                cbind(tapply(var, grouping, mean), tapply(var, grouping, sd))
        } 
        else LaplaceEst(table(grouping, var), f = fL)
    }
    
    tables <- lapply(x, est)
    
    if(!usekernel){
        num <- sapply(x, is.numeric)
        temp <- as.matrix(sapply(tables, function(x) x[,2]))
        temp[,!num] <- 1
        temp <- apply(temp, 2, function(x) any(!x))
        if(any(temp))
            stop("Zero variances for at least one class in variables: ", 
                paste(names(tables)[temp], collapse=", "))
    }
    names(dimnames(apriori)) <- Yname
    structure(list(apriori = apriori, tables = tables, levels = levels(grouping), 
        call = call, x = x, usekernel = usekernel, varnames = colnames(x)), 
        class = "NaiveBayes")
}
