#################################################################
###
### WOE
###
### Functions to perform WOE transforms factor ~> numeric 
### (for categorical variables and a binary target)
###
###############################################################

### (sub-)functions:
# compute woe - coefficients for a single factor variable x
computewoes <- function(x, y, adj){
  if(! is.factor(x)) stop("Input variable x must be a factor!")
  if(! is.factor(y)) stop("Target variable y must be a factor!")
  if(sum(table(x,y)==0) > 0) warning("At least one empty cell (class x level) does exist")
  # class wise event rates 
  xtab <- table(y, x)
  # correction for empty class levels
  xtab[which(xtab == 0)] <- xtab[which(xtab == 0)] + adj
  ncl      <- table(x)
  # compute woes for alle classes
  fxgegy <- xtab
  fy  <- rowSums(fxgegy)
  for(i in 1:2) fxgegy[i,] <- fxgegy[i,] / fy[i] 
  woes <- log(fxgegy[1,]/fxgegy[2,])
  if(any(fxgegy[2,] == 0)) warning("Empty cells result in infinite woes, zeroadj should be specified > 0!")
  
  # add IV 
  # difference of class wise relative frequencies
  bdiff <- fxgegy[1,] - fxgegy[2,]
  # calculate information value
  IV <- sum(bdiff * woes)
  woes <- c(woes, IV)
  }

# apply woe transformation to one single variable using prespecified woe coefficients 
applywoes <- function(woe.obj, x.vec){
    # woe.obj: woe coefficients as returned from compute woes => a single alement of train.woes
    if ( sum(sapply(unique(x.vec), function(x) return(sum(x == unique(names(woe.obj))) == 0))) > 0 ) stop("Factor Levels do not match!")
    xwoe <- sapply(x.vec, function(z) return(woe.obj[which(names(woe.obj) == z)]))  
    return(as.numeric(xwoe))
    }


### main WOE functions
# 'train' woe - coefficients for all factor variables of a data frame
woe.default <- function(x, grouping, zeroadj = 0, ids = NULL, appont = TRUE, ...){
  if (is.factor(x)){ 
    warning("Only one single input variable. Variable name in resulting object$woe is only conserved in formula call.")
    x <- as.data.frame(x)
    }
  if (!is.data.frame(x)) stop("x should be of type data frame.")
  if(is.numeric(ids)){
    if(max(ids) > ncol(x)) stop("Uncorrect coloumn ids specified.")    
    fact.ids <- ids
    } 

  if(is.character(ids)){
    fact.ids <- which(colnames(x) %in% ids)
    if (length(fact.ids) < 1) stop("Uncorrect variable names specified!")
  }
  
  if(is.null(ids)) fact.ids <- which(sapply(x,is.factor))

  if (length(fact.ids) < 1) stop("No factor variables to be transformed!")
  if(!all(sapply(x[,ids], is.factor))) stop("ids should specify only factor variables")
  if (length(table(grouping)) != 2) stop("WOE transformation is only for binary targets!")

  x.fact <- x[, fact.ids]
  # case of more than factor variable to be transformed
  if(length(fact.ids) > 1){
    x.woes <- lapply(x.fact, computewoes, y = grouping, adj = zeroadj)
    # separate woes and IVs
    IVs <-  sapply(x.woes, function(x) return(x[length(x)]))
    x.woes <- lapply(x.woes, function(x) return(x[-length(x)]))
  }
  # case of only one factor variable to be transformed
  if(length(fact.ids)==1){
    x.woes <- computewoes(x=x.fact, y=grouping, adj = zeroadj)
    # separate woes and IVs
    IVs <-  x.woes[length(x.woes)]
    x.woes <- x.woes[-length(x.woes)]
    x.woes <- list(x.woes)
    names(x.woes) <- names(x)[fact.ids]
    }
  
  res <- list("woe" = x.woes, "IV" = IVs)
  class(res) <- "woe"
  
  if(appont){
    xnew <- predict(res, x, ...)
    res$xnew <- xnew
    }
  
  return(res)
  }

# definition of method woe
woe <- function (x, ...) 
    UseMethod("woe")

# formula interface
woe.formula <- function(formula, data = NULL, ...)
{
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval.parent(m$data))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.frame(Terms, m)
    #cat(as.character(attr(Terms, "variables"))[1:2])    
    x <- x[,-which(names(x) == as.character(attr(Terms, "variables"))[2])] 
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
    #return(list(x, grouping, xvars))
    res <- woe(x, grouping, ...)
    res$terms <- Terms
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    # if x is a single factor variable its name has to be stored within model
    if(is.factor(x)) names(res$woe)[1] <- names(xlev)
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        res$na.action <- attr(m, "na.action")
    res
}


# predict woes to a data set
predict.woe <- function(object, newdata, replace = TRUE, ...){
  
  if (!is.data.frame(newdata)) stop("newdata should be of type data frame.")
  # object:
  object <- object$woe
  if (sum(sapply(newdata, is.factor)) == 0) stop("No factor variables to be transformed!")
  #fact.ids <- unlist(sapply(names(newdata), function(x) return(which(names(object) == x))))
  fact.ids <- unlist(sapply(names(object), function(x) return(which(names(newdata) == x))))

  if (length(fact.ids) < 1) stop("No factor variables to be transformed!")  

  x.fact <- newdata[, fact.ids]
  
  if(length(fact.ids) > 1){
    x.woes <- as.data.frame(sapply(seq_along(fact.ids), 
            function(i) return(applywoes(object[[i]], x.fact[,which(names(x.fact) == names(object)[i])])))) 
  
    names(x.woes) <- paste("woe", names(fact.ids),sep=".")
    x <- data.frame(newdata, x.woes)
    }
  # special case of only one variable (vector) to be transformed
  if(length(fact.ids) == 1){
    x.woes <- applywoes(object[[1]], x.fact) 
    x <- data.frame(newdata, x.woes)
    names(x)[ncol(x)] <- paste("woe", names(fact.ids),sep=".")
  }
  
  if (replace) x <- x[,-fact.ids]  
  return(x)
  }


plot.woe <- function(x, type = c("IV", "woes"), ...){
  type <- match.arg(type)
  if(type=="IV"){
    x <- x$IV
    barplot(height = x, ylab = "Information value", ...)
    xmax <- length(x) * 1.2
    lines(c(0, xmax), rep(0.02, 2), col = "red") 
    lines(c(0, xmax), rep(0.10, 2), col = "yellow") 
    lines(c(0, xmax), rep(0.30, 2), col = "green") 
  }
  if(type=="woes"){
    if(!any(names(x) == "xnew")) stop("Plot of type 'woe' only possible for appont == TRUE.")
    vids <- which(substr(names(x$xnew),1,4) == "woe.")
    if(length(vids) == 0) stop("Data contains no woe variables to be displayed.")
    
    par(ask="TRUE")
    n <- nrow(x$xnew)
    for(i in vids){
      tab <- table(x$xnew[,i])/n
      plot(as.numeric(names(tab)), tab, type = "h", xlab = names(x$xnew)[i], 
           ylim=c(0,1), yaxt = "n", ylab ="Relative Frequency", ...)
      axis(2, at = seq(0, 1, 0.2))  
      }
    par(ask="FALSE")
    invisible()
  }
}


print.woe <- function(x, ...){
    cat("Information values of transformed variables: \n\n")
    IV   <- sort(x$IV, decreasing = TRUE)
    names <- names(IV)
    print(cbind(IV))
}
