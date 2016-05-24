workhorse.inbagg <- function(object, y, X, W, 
  cFUN, w.training.set, y.training.set, bcontrol, control, ...)
{
  formula.list <- object
  data <- data.frame(y, X, W)
  mtrees <- vector(mode="list", length=bcontrol$nbagg)
  if(w.training.set[1] == "all") fit.vals <- 1:length(y)

  for (i in 1:bcontrol$nbagg) {
    bindx <- sample(1:length(y), bcontrol$ns, replace=bcontrol$replace)
    if(w.training.set[1] == "oob") fit.vals <- (-bindx)
    if(w.training.set[1] == "bag") fit.vals <- bindx

    objs <- vector(mode="list", length=length(formula.list))	#prediction models for intermediate variables
    names(objs) <- names(formula.list)

    addclass <- function() {					##START addclass <- function()
      for (j in 1:length(formula.list)) {			##Fitting prediction models for intermediates
        oX <- data[fit.vals, c(paste(formula.list[[j]]$formula[[2]]), 
                   attr(terms(formula.list[[j]]$formula, dataa = data), "term.labels"))]
        foo <- try(formula.list[[j]]$model(formula.list[[j]]$formula, data = oX))
        objs[[j]] <- foo
      }
     
      fct <- function(newdata) {				##START fct <- function(newdata)
        if (!is.data.frame(newdata))
          newdata <- as.data.frame(newdata)
        add.predictors <- rep(0, nrow(newdata))

        for (j in 1:length(formula.list)){			## predict additional intermediates using fitted models
          oXnewdata <- newdata[,attr(terms(formula.list[[j]]$formula, data = data), "term.labels")]
          if(is.null(formula.list[[j]]$predict)) {
            res <- try(predict(objs[[j]], newdata  = oXnewdata))
          } else {
            res <- try(formula.list[[j]]$predict(objs[[j]], newdata  = oXnewdata))
            }
###FIX: action for class(res) == "try-error"
          add.predictors <- data.frame(add.predictors, res)
        }
        add.predictors <- add.predictors[,-1]
        if(is.null(dim(add.predictors))) add.predictors <- matrix(add.predictors, ncol = 1)
        colnames(add.predictors) <- names(formula.list)
        add.predictors
      }
        					##END fct <- function(newdata)      
      return(fct)
    }						##END addclass <- function()


    bfct <- addclass()				###bfct is a function (addclass)

    if (!is.null(bfct)) {
      expl.cFUN <- attr(terms(as.formula(cFUN$formula), data = data), "term.labels")

      if(!is.null(cFUN$fixed.function)) {
         btree <- cFUN
      } else {
        W.new <- bfct(X)
        W.new.names <- sub(".[0-9]$", "", colnames(W.new))

        if(y.training.set[1] == "fitted.bag") {	###contstruct on bag
          oX <- data.frame(y, X, W.new)[bindx,]
          right.side <- paste(c(expl.cFUN[!(expl.cFUN %in% W.new.names)], colnames(W.new)[W.new.names %in% expl.cFUN]), collapse = "+")
          cFUN$formula <- as.formula(paste(cFUN$formula[[2]], "~", right.side))
        }
		
        if(y.training.set[1] == "original") {	###construct on original variables
          if(length(W.new.names)> length(colnames(W))) stop("If classifying function is trained on original intermediate, only one predictive model per intermediate can be constructed.")
          oX <- data.frame(y, X, W[,W.new.names])
          names(oX)[(ncol(oX)-ncol(W)+1):ncol(oX)] <- colnames(W.new)
        }

        if(y.training.set[1] == "fitted.subset") {		###construct on subset
          oX <- data.frame(y, X, W.new)[!subset,]		
          right.side <- paste(c(expl.cFUN[!(expl.cFUN %in% W.new.names)], colnames(W.new)[W.new.names %in% expl.cFUN]), collapse = "+")
          cFUN$formula <- as.formula(paste(cFUN$formula[[2]], "~", right.side))
        }
        names(oX)[names(oX) == "y"] <- paste(cFUN$formula[[2]])
        btree <- cFUN$model(cFUN$formula, data = oX, ...)
        btree <- list(model = btree, predict = cFUN$predict)
      }

      this <- list(bindx = bindx, btree = btree, bfct=bfct)
    } else {
      stop("Predictive function for intermediates not executable: Classifying function can not be applied.")
    }
    class(this) <- "thisclass"
    mtrees[[i]] <- this
  }
  mtrees
}




inbagg <- function(formula, data, ...) UseMethod("inbagg", data)

inbagg.default <- function(formula, data,...)
{
  stop(paste("Do not know how to handle objects of class", class(data)))
}  


inbagg.data.frame <- function(formula, data, pFUN=NULL, 
 cFUN=list(model = NULL, predict = NULL, training.set = NULL), 
 nbagg = 25, ns = 0.5, replace = FALSE, ...)
{
  if(!is.function(cFUN)) {
    if(is.null(cFUN$model)) {
    cFUN$model <-  function(formula, data) 
            rpart(formula, data, control = rpart.control(minsplit=2, cp=0, xval=0))
    if(is.null(cFUN$predict)) cFUN$predict <- function(object, newdata) predict(object, newdata, type = "class")
    if(is.null(cFUN$training.set))  cFUN$trainig.set <- "fitted.bag"
    }
 }

##check formula
  if(missing(formula)
    || (length(formula) != 3)
    || (length(attr(terms(formula[-2], data = data), "term.labels")) < 1))
    stop("formula missing or incorrect")

  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)

##editing formula
  if(length(formula[[2]])==3) {
    if(!is.function(cFUN)) {
        if (is.null(cFUN$formula)) y.formula <- as.formula(formula[[2]]) else y.formula <- cFUN$formula
    }
  
    w.formula <- XX~YY
    w.formula[[2]] <- formula[[2]][[3]]
    w.formula[[3]] <- formula[[3]]

    response <-  paste(formula[[2]][[2]])
    w.names <- attr(terms(as.formula(formula[[2]]), data = data), "term.labels")
    x.names <- attr(terms(as.formula(formula), data = data), "term.labels")

    if(length(x.names == 1) && x.names == ".") x.names <- colnames(data)[!(colnames(data) %in% c(response, w.names))]
    y <- data[, response]
    X <- data[, x.names]
    W <- data[, w.names]
    if(is.null(dim(X))) X <- matrix(X, ncol = 1, dimnames = list(rownames(W), x.names))
    if(is.null(dim(W))) W <- matrix(W, ncol = 1, dimnames = list(rownames(X), w.names))
                                    
    if(is.function(cFUN)) {
      y.formula <- as.formula(paste(formula[[2]][[2]], "~", paste(c(x.names, w.names), collapse = "+")))
      fixed.function <- cFUN
      cFUN <- list()
      cFUN$fixed.function <- fixed.function
    }
   cFUN$formula <- y.formula

  } else {
    stop(paste("Specified formula has to be of type y~x~w"))   
  }
##remove settings of training.set
  if(is.null(pFUN$training.set)) w.training.set <- "oob" else w.training.set <- pFUN$training.set[1]
  pFUN$training.set <- NULL

  if(is.null(cFUN$training.set)) y.training.set <- "fitted.bag" else y.training.set <- cFUN$training.set[1]
  cFUN$training.set <- NULL

  bcontrol <- list(nbagg = nbagg, ns = length(y)*ns, replace = replace)

  if(is.null(w.formula)) stop("no formula for prediction model specified")

  ##formula.list : list of lists which specify an abitrary number of models for intermediate variables:
  ##w1.1, w2.1, w3.1, ...., w2.1, w2.2, w3.1, .... where 'w*' is the variable and '.*' describes the model

  P <- length(pFUN)
  number.models <- c() 
  for(i in 1:P) {
    if(is.null(pFUN[[i]]$formula)) pFUN[[i]]$formula <- w.formula
    number.models <- c(number.models, 
                       paste(attr(terms(pFUN[[i]]$formula[-3], data = data), "term.labels"), 
                             ".", i, sep = ""))
  }

  formula.list <- vector(mode = "list", length= length(number.models))
  names(formula.list) <- paste(number.models)

  for(i in 1:P) {  
    res <- list()  
    Qi <- length(attr(terms(pFUN[[i]]$formula[-3], data = data), "term.labels"))
    for(j in 1:Qi) {
      res$formula <- w.formula
      res$formula[[2]] <- as.name(attr(terms(res$formula[-3], data = data), "term.labels")[j])
      res$formula[[3]] <- pFUN[[i]]$formula[[3]]

      if(res$formula[[3]] == ".") res$formula <- as.formula(paste(res$formula[[2]], "~", paste(x.names, collapse= "+")))
      res$model <- pFUN[[i]]$model
      res$predict <- pFUN[[i]]$predict
      formula.list[[paste(res$formula[[2]], ".", i, sep = "")]] <- res
    }       
  }    

##apply
  res <- workhorse.inbagg(object = formula.list, y = y, X = X, W = W, 
    cFUN = cFUN, w.training.set = w.training.set, y.training.set = y.training.set, 
    bcontrol = bcontrol, ...)
  RET <- list(mtrees = res, y = y, W = W, X = X)
  class(RET) <- "inbagg"
  RET
}
 

print.inbagg <- function(x, ...)
{
  q <- length(x$mtrees)
  intermediates <- attr(x$W, "names")
  text.intermediates <- paste("Indirect bagging, with", q, 
    "bootstrap samples and intermediate variables: \n", 
    paste(intermediates, collapse = " "))
  cat("\n", text.intermediates, "\n")
}


summary.inbagg <- function(object, ...)
{
  class(object) <- "summary.inbagg"
  object
}


print.summary.inbagg <- function(x, ...)
{
  q <- length(x$mtrees)
  intermediates <- attr(x$W, "names")
 
  text.intermediates <- paste("Indirect bagging, with", q,
"bootstrap samples and intermediate variables:", paste(intermediates, collapse = " "))

  cat("\n", text.intermediates, "\n")
  for(i in 1:length(x)) {
    print(x$mtrees[[i]])
  }
}

