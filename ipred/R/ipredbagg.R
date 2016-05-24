#$Id: ipredbagg.R,v 1.13 2003/06/11 10:40:17 peters Exp $

workhorse <- function(y, X, control, comb, bcontrol, thisclass, ...) {
  # This is double-bagging (comb is lda) or bundling (any arbritrary
  # model in comb)
  if (!is.data.frame(X)) X <- as.data.frame(X)

  # check user supplied functions
  if (!is.list(comb)) stop("comb not a list")

  N <- nrow(X)

  mydata <- cbind(data.frame(y), X)
  mtrees <- vector(mode="list", length=bcontrol$nbagg)

  for (i in 1:bcontrol$nbagg) {
    # double-bagging or bundling
    # comb is a list of lists, each of them having two elements:
    # model and predict

    bindx <- sample(1:N, bcontrol$ns, replace=bcontrol$replace)

    objs <- vector(mode="list", length=length(comb))
    addclass <- function() {
      myindx <- 1:length(comb)
      for (k in 1:length(comb)) {
        # put the user supplied models into a try statement
        # if this fails, simply ignore it.
        # options(show.error.messages = FALSE)
        oX <- mydata[-bindx,]
        foo <- try(comb[[k]]$model(y ~ ., data=oX))
        if (inherits(foo, "try-error")) {
          warning("could not build model:")
          print(foo[1])
          foo <- NA
          myindx <- myindx[-k]
        } 
        objs[[k]] <- foo
        # options(show.error.messages = TRUE)
      }
      fct <- function(newdata) {
        # use lexical scoping: return this function for the computation of 
        # the additional predictors
        if (!is.data.frame(newdata))
          newdata <- as.data.frame(newdata)
        addpred <- c()
        # the user supplied model failed, ignore it here.
        if (length(myindx) < 1) {
          RET <- NULL 
        } else {
          # compute additional predictors for user supplied models
          for (k in myindx)
            addpred <- cbind(addpred, comb[[k]]$predict(objs[[k]], newdata))
          # <FIXME>: more informative names???
          colnames(addpred) <- paste("addpred", 1:ncol(addpred), sep="")
          # </FIXME>
          RET <- addpred
        }
        RET
      }
      if (length(myindx) < 1) return(NULL) else return(fct)
    }
    bfct <- addclass()
    # may have failed
    if (!is.null(bfct)) {
      # grow a tree using the original predictors
      # from the bootstrap sample and the additional predictors computed on
      # the bootstrap sample.
      oX <- cbind(mydata, bfct(X))[bindx,]
      btree <- rpart(y ~., data=oX, control = control,...)
      # return this object
      this <- list(bindx = bindx, btree = btree, bfct=bfct)
    } else {
      # return a simple tree if the user supplied model failed.
      oX <- mydata[bindx,]
      btree <- rpart(y ~., data=oX, control = control,...)
      this <- list(bindx = bindx, btree = btree)
    }
    class(this) <- thisclass
    mtrees[[i]] <- this
  }
  mtrees
}


ipredbagg <- function(y, ...) {
  if(is.null(class(y))) 
    class(y) <- data.class(y)
#  UseMethod("ipredbagg", y, ...)
  UseMethod("ipredbagg", y)
}

ipredbagg.default <- function(y, ...) {
  stop(paste("Do not know how to handle objects of class", class(y)))
}

ipredbagg.integer <- function(y, ...) {
  ipredbagg.numeric(y,...)
}


ipredbagg.factor <- function(y, X=NULL, nbagg=25, control=
                             rpart.control(minsplit=2, cp=0, xval=0), 
                             comb=NULL, coob=FALSE, ns=length(y), keepX = 
                             TRUE, ...) {
  # bagging classification trees

  if (!is.null(comb) && coob) 
    stop("cannot compute out-of-bag estimate for combined models")

  if (nbagg == 1 && coob) 
    stop("cannot compute out-of-bag estimate for single tree")

  # check nbagg
  if (nbagg < 1) stop("nbagg is not a positive integer")
  # bagging only if nbagg greater 1, else use the whole sample, i.e. one
  # simple tree
  if (nbagg == 1) { 
    REPLACE <- FALSE 
  } else {
    if (ns < length(y)) {
      # this is "subagging", i.e. sampling ns out of length(y) WITHOUT
      # replacement
      REPLACE <- FALSE 
    } else {
      # the usual bootstrap: n out of n with replacement
      REPLACE <- TRUE
    }
  }

  if (!is.null(comb)) {
    # this is rather slow but we need to be as general as possible
    # with respect to classifiers as well as outcome of prediction (classes,
    # linear discriminant functions, conditional class probabilities, random
    # noise, if you like)
    mtrees <- workhorse(y, X, control, comb,
                        bcontrol=list(nbagg=nbagg, ns=ns, replace=REPLACE),
                        thisclass="sclass")
  } else {
    # use an optimized version
    mydata <- cbind(data.frame(y), X)
    mtrees <- irpart(y ~ ., data=mydata, control=control,
                     bcontrol=list(nbagg=nbagg, ns=ns, replace=REPLACE))
  }
  # always keep response and predictors as well as a list of nbagg objects
  # of class "sclass" 
  if (keepX) 
    RET <- list(y=y, X=X, mtrees=mtrees, OOB=coob, comb=!is.null(comb))
  else 
    RET <- list(y=y, X=NULL, mtrees=mtrees, OOB=coob, comb=!is.null(comb))
  class(RET) <- "classbagg"

  if (coob) {
    pred <- predict(RET)
    ae <- all.equal(levels(pred), levels(RET$y))
    if (is.logical(ae) && ae)
       RET$err <- mean(pred != RET$y, na.rm=TRUE)
    else
       RET$err <- mean(as.character(pred) != as.character(RET$y), 
                       na.rm=TRUE)
   }
   RET
}

ipredbagg.numeric <- function(y, X=NULL, nbagg=25, control=
                             rpart.control(xval=0), 
                             comb=NULL, coob=FALSE, ns=length(y), keepX =
                             TRUE, ...) {
  # <FIXME> is control meaningful here??? </FIXME>

  # bagging regression trees

  if (!is.null(comb) && coob) 
    stop("cannot compute out-of-bag estimate for combined models")

  if (nbagg == 1 && coob) 
    stop("cannot compute out-of-bag estimate for single tree") 

  # check nbagg
  if (nbagg < 1) stop("nbagg is not a positive integer")
  # only bagg if nbagg greater 1, else use the whole sample 
  if (nbagg == 1) {
    REPLACE <- FALSE
  } else {
    if (ns < length(y)) {
      # this is "subagging", i.e. sampling ns out of length(y) WITHOUT
      # replacement
      REPLACE <- FALSE
    } else {
      # the usual bootstrap: n out of n with replacement
      REPLACE <- TRUE
    }
  }

  if (!is.null(comb)) {
    mtrees <- workhorse(y, X, control, comb,
                        bcontrol=list(nbagg=nbagg, ns=ns, replace=REPLACE),
                        thisclass="sreg")
  } else {
    mydata <- cbind(data.frame(y), X)
    mtrees <- irpart(y ~ ., data=mydata, control=control,
                     bcontrol=list(nbagg=nbagg, ns=ns, replace=REPLACE))
  }

  if (keepX) 
    RET <- list(y=y, X=X, mtrees=mtrees, OOB=coob, comb=!is.null(comb))
  else 
    RET <- list(y=y, X=NULL, mtrees=mtrees, OOB=coob, comb=!is.null(comb))
  class(RET) <- "regbagg"

  if (coob)
    RET$err <- sqrt(mean((predict(RET) - RET$y)^2, na.rm=TRUE))
  RET
}


ipredbagg.Surv <- function(y, X=NULL, nbagg=25, control=
                             rpart.control(xval=0), 
                             comb=NULL, coob=FALSE, ns=dim(y)[1], keepX =
                             TRUE, ...) {
  # <FIXME> is control meaningful here??? </FIXME>

  # bagging survival trees

  if (!is.null(comb) && coob) 
    stop("cannot compute out-of-bag estimate for combined models")

  if (nbagg == 1 && coob) 
    stop("cannot compute out-of-bag estimate for single tree") 

  # check nbagg
  if (nbagg < 1) stop("nbagg is not a positive integer")
  # only bagg if nbagg greater 1, else use the whole sample 
  if (nbagg == 1) {
    REPLACE <- FALSE
  } else {
    if (ns < dim(y)[1]) {
      # this is "subagging", i.e. sampling ns out of length(y) WITHOUT
      # replacement
      REPLACE <- FALSE
    } else {
      # the usual bootstrap: n out of n with replacement
      REPLACE <- TRUE
    }
  }

  if (!is.null(comb)) {
    mtrees <- workhorse(y, X, control, comb,
                        bcontrol=list(nbagg=nbagg, ns=ns, replace=REPLACE),
                        thisclass="ssurv")
  } else {
    mydata <- cbind(data.frame(y), X)
    mtrees <- irpart(y ~ ., data=mydata, control=control,
                     bcontrol=list(nbagg=nbagg, ns=ns, replace=REPLACE))
  }
  if (keepX) 
    RET <- list(y=y, X=X, mtrees=mtrees, OOB=coob, comb=!is.null(comb))
  else 
    RET <- list(y=y, X=NULL, mtrees=mtrees, OOB=coob, comb=!is.null(comb))
  class(RET) <- "survbagg"
  
  if (coob) 
    RET$err <- sbrier(RET$y, predict(RET))
  RET
}

