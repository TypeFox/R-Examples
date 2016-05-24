# $Id: predict.bagging.R,v 1.17 2009/03/27 16:18:38 hothorn Exp $

uwhich.max <- function(x) {
  # need to determine all maxima in order to sample from them
  wm <- (1:length(x))[x == max(x)]
  if (length(wm) > 1)
    wm <- wm[sample(length(wm), 1)]
  wm
}

predict.classbagg <- function(object, newdata=NULL, type=c("class", "prob"),
                              aggregation=c("majority", "average", "weighted"), ...) {
  type <- match.arg(type)
  agg <- match.arg(aggregation)
  if (missing(newdata)) {
    if (length(object$mtrees) < 10) 
      stop("cannot compute out-of-bag predictions for small number of trees")
    OOB <- TRUE
    if (!is.null(object$X))
      newdata <- object$X
    else
      stop("cannot compute out-of-bag predictions without object$X!")
  } else {
    OOB <- FALSE
  }
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  N <- nrow(newdata)
  if (!object$comb) {
    tree <- object$mtrees[[1]]$btree
    Terms <- delete.response(tree$terms)
    act <- (tree$call)$na.action
    if (is.null(act)) act<- na.rpart
    newdata <- model.frame(Terms, newdata, na.action = act,
                           xlev=attr(tree, "xlevels"))
    newdata <- getFromNamespace("rpart.matrix", ns = "rpart")(newdata)
  }
  classes <- levels(object$y)
  switch(agg, "majority" = {
    vote <- matrix(0, nrow=N, ncol=length(classes))
    for (i in 1:length(object$mtrees)) {
      if (OOB) {
        bindx <- object$mtrees[[i]]$bindx
        if (!is.null(object$mtrees[[i]]$bfct)) 
          stop("cannot compute out-of-bag estimate for combined models!")
        pred <- predict(object$mtrees[[i]], newdata, type="class")
        tindx <- cbind((1:N), pred)[-bindx,]
      } else {
        tindx <- cbind(1:N, predict(object$mtrees[[i]], newdata,
                                    type="class"))
      }
      vote[tindx] <- vote[tindx] + 1
    }
    if (type=="class") {
      RET <- factor(classes[apply(vote, 1, uwhich.max)])
    } else {
      RET <- vote/apply(vote, 1, sum)
      colnames(RET) <- classes
    }
  }, 
  "average" = {
    cprob <- matrix(0, nrow=N, ncol=length(classes))
    if (OOB) ncount <- rep(0,N) else ncount <- length(object$mtrees)
    for (i in 1:length(object$mtrees)) {
      if (OOB) {
        bindx <- object$mtrees[[i]]$bindx
        pred <- predict(object$mtrees[[i]], newdata, type="prob")[-bindx,]
        tindx <- (1:N)[-bindx]
        ncount[tindx] <- ncount[tindx] + 1
      } else {
        pred <- predict(object$mtrees[[i]], newdata, type="prob")
        tindx <- 1:N
      }
      cprob[tindx,] <- cprob[tindx,] + pred
    }
    switch(type, "class" = {
      RET <- as.factor(apply(cprob, 1, uwhich.max))
      levels(RET) <- classes
    }, 
    "prob" = {
      ncount[ncount < 1] <- NA
      RET <- cprob / ncount
      colnames(RET) <- classes
    })
  },
  "weighted" = {
    agglsample <- matrix(0, ncol=length(classes), nrow=N)
    for (i in 1:length(object$mtrees)) {
      bdata <- object$y[object$mtrees[[i]]$bindx]
      newpart <- getpartition(object$mtrees[[i]], newdata)
      oldpart <- object$mtrees[[i]]$btree$where
      if (OOB)
        tindx <- (1:N)[-object$mtrees[[i]]$bindx]
      else
        tindx <- 1:N
      for (j in tindx) {
        aggobs <- table(bdata[oldpart == newpart[j]])
        agglsample[j,] <- agglsample[j,] + aggobs
      }
    }
    switch(type, "class" = {
      RET <- c()
      for (j in 1:N)
        RET <- as.factor(c(RET, uwhich.max(agglsample[j,])))
      levels(RET) <- classes
    },
    "prob" = {
      RET <- agglsample / apply(agglsample, 1, sum)
      colnames(RET) <- classes
    })
  })
  RET
}

predict.sclass <- function(object, newdata=NULL, type=c("class", "prob"),
...) {
  if (!is.null(object$bfct))
    newdata <- cbind(newdata, object$bfct(newdata))
  pred <- predict.irpart(object$btree, newdata, type=type)
  RET <- pred
  if (type == "class") RET <- as.integer(pred)
  if (type == "prob" && is.vector(pred)) RET <- cbind(pred, 1 - pred)
  RET
}


predict.regbagg <- function(object, newdata=NULL, aggregation=c("average",
"weighted"), ...) {
  agg <- match.arg(aggregation)
  if (missing(newdata)) {
    if (length(object$mtrees) < 10) 
      stop("cannot compute out-of-bag predictions for small number of trees")
    OOB <- TRUE
    if (!is.null(object$X))
      newdata <- object$X
    else 
      stop("cannot compute out-of-bag predictions without object$X!")
  } else {
    OOB <- FALSE
  }
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  N <- nrow(newdata)
  if (!object$comb) {
    tree <- object$mtrees[[1]]$btree
    Terms <- delete.response(tree$terms)
    act <- (tree$call)$na.action
    if (is.null(act)) act<- na.rpart
    newdata <- model.frame(Terms, newdata, na.action = act,
                           xlev=attr(tree, "xlevels"))
    newdata <- getFromNamespace("rpart.matrix", ns = "rpart")(newdata)
  }
  switch(agg, "average" = {
    cprob <- rep(0, N)
    if (OOB) ncount <- rep(0,N) else ncount <- length(object$mtrees)
    for (i in 1:length(object$mtrees)) {
      if (OOB) {
        bindx <- object$mtrees[[i]]$bindx
        if (!is.null(object$mtrees[[i]]$bfct))
          stop("cannot compute out-of-bag estimate for combined models!")
        pred <- predict(object$mtrees[[i]], newdata)[-bindx]
        tindx <- (1:N)[-bindx]
        ncount[tindx] <- ncount[tindx] + 1
      } else {
        pred <- predict(object$mtrees[[i]], newdata)
        tindx <- 1:N
      }
      cprob[tindx] <- cprob[tindx] + pred
    }
    ncount[ncount < 1] <- NA
    RET <- cprob / ncount
  },
  "weighted" = {
    agglsample <- rep(0, N)
    ncount <- rep(0, N)
    for (i in 1:length(object$mtrees)) {
      bdata <- object$y[object$mtrees[[i]]$bindx]
      newpart <- getpartition(object$mtrees[[i]], newdata)
      oldpart <- object$mtrees[[i]]$btree$where
      if (OOB)
        tindx <- (1:N)[-object$mtrees[[i]]$bindx]
      else
        tindx <- 1:N
      for (j in tindx) {
        aggobs <- bdata[oldpart == newpart[j]]
        agglsample[j] <-  agglsample[j] + sum(aggobs)
        ncount[j] <- ncount[j] + length(aggobs)
      }
    }
    ncount[ncount < 1] <- NA
    RET <- agglsample / ncount
  })
  RET
}


predict.sreg <- function(object, newdata=NULL, ...) {
  if (!is.null(object$bfct))
    newdata <- cbind(newdata, object$bfct(newdata))
  predict.irpart(object$btree, newdata)
}


predict.survbagg <- function(object, newdata=NULL, ...) {
  if (missing(newdata)) {
    if (length(object$mtrees) < 10) 
      stop("cannot compute out-of-bag predictions for small number of trees")
    OOB <- TRUE
    if (!is.null(object$X))
      newdata <- object$X
    else 
      stop("cannot compute out-of-bag predictions without object$X!")
  } else {
    OOB <- FALSE
  }
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  N <- nrow(newdata)
  if (!object$comb) {
    tree <- object$mtrees[[1]]$btree
    Terms <- delete.response(tree$terms)
    act <- (tree$call)$na.action
    if (is.null(act)) act<- na.rpart
    newdata <- model.frame(Terms, newdata, na.action = act,
                           xlev=attr(tree, "xlevels"))
    newdata <- getFromNamespace("rpart.matrix", ns = "rpart")(newdata)
  }
  agglsample <- list()
  aggcens <- list()
  for (j in 1:N) { 
    agglsample <- c(agglsample, list(c()))
    aggcens <- c(aggcens, list(c()))
  }
  for (i in 1:length(object$mtrees)) {
    bdata <- object$y[object$mtrees[[i]]$bindx]
    newpart <- getpartition(object$mtrees[[i]], newdata)
    oldpart <- object$mtrees[[i]]$btree$where
    if (OOB) {
      if (!is.null(object$mtrees[[i]]$bfct))
        stop("cannot compute out-of-bag estimate for combined models!")
      tindx <- (1:N)[-object$mtrees[[i]]$bindx]
    } else {
      tindx <- 1:N
    }
    for (j in tindx) {
        aggobs <- bdata[oldpart == newpart[j],1]
        agglsample[[j]] <- c(agglsample[[j]], aggobs)
        aggobs <- bdata[oldpart == newpart[j],2]
        aggcens[[j]] <- c(aggcens[[j]], aggobs)
    }
  }
  RET <- list()
  for (j in 1:N)
    RET <- c(RET, list(survfit(Surv(agglsample[[j]], aggcens[[j]]) ~ 1)))
  RET
}

getpartition <- function(object, newdata=NULL) {
  if (!is.null(object$bfct)) {
    newdata <- cbind(newdata, object$bfct(newdata))
    Terms <- delete.response(object$btree$terms)
    act <- (object$btree$call)$na.action
    if (is.null(act)) act<- na.rpart
    newdata <- model.frame(Terms, newdata, na.action = act,
                             xlev=attr(object$btree, "xlevels"))
    newdata <- getFromNamespace("rpart.matrix", ns = "rpart")(newdata)
  }
  getFromNamespace("pred.rpart", ns = "rpart")(object$btree, newdata)
}

