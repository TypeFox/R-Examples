library(e1071)
library(catnet)
library(lars)

classification.cv10fold <- function(cdata, cls, pathway, ddata = NULL, npars=2, nodePars = NULL, ncats = 3, nodeCats = NULL, nodeOrder=NULL, ker="radial", gam=0.01, folds=10, echo = FALSE) {

  if(ncol(cdata)*nrow(cdata) <= 0)
    return(list(nb=0, catnet0=0, catnet.bic=0, catnet.aic=0, svm=0, lasso.cp=0, lar.cp=0, time=0))
  
  data <- ddata
  if(is.null(data))
     data <- cnDiscretize(cdata, ncats, mode="uniform")
  
  if(is.null(nodeOrder))
    nodeOrder <- 1:nrow(cdata)

  clslevs <- range(cls)

  ind <- sample(1:ncol(cdata), ncol(cdata))

  numnodes <- nrow(cdata)
  
  t.start <- proc.time()
  bres0 <- NULL
  bres1 <- NULL
  bres2 <- NULL
  bres4 <- NULL
  bres8 <- NULL
  bres9 <- NULL
  bres10 <- NULL
  j <- 1
  while(j < ncol(cdata)) {
    jj <- floor(j + ncol(cdata)/folds)
    if(jj > ncol(cdata))
      jj <- ncol(cdata)
    ##cat(ind[j:jj],"\n")

    train.data <- data[,-ind[j:jj]]
    dat1 <- train.data[, cls[-ind[j:jj]]==clslevs[1]]   
    dat2 <- train.data[, cls[-ind[j:jj]]==clslevs[2]]
    
    ##########################################################
    ## naive bayes

    bnodes <- c("class", pathway)
    bpars <- vector("list", (numnodes+1))
    bcats <- vector("list", (numnodes+1))
    bcats[[1]] <- clslevs
    for(ni in 1:numnodes) {
      bpars[[ni+1]] <- c(1)
      bcats[[ni+1]] <- nodeCats[[ni]]
    }
    bnet <- cnNew(bnodes, bcats, bpars)
    bnet <- cnSetProb(bnet, rbind(class=cls[-ind[j:jj]], train.data), nodeCats=bcats) 

    for(k in ind[j:jj]) {
      rr <- matrix(c(class=clslevs[1],data[,k]), ncol=1)
      if(nrow(rr) == length(bnet@nodes))
        rownames(rr) <- bnet@nodes
      l1 <- cnNodeLoglik(bnet, sapply(1:bnet@numnodes, function(i) return(bnet@nodes[i])), rr)
      rr <- matrix(c(class=clslevs[2],data[,k]), ncol=1)
      if(nrow(rr) == length(bnet@nodes))
        rownames(rr) <- bnet@nodes
      l2 <- cnNodeLoglik(bnet, sapply(1:bnet@numnodes, function(i) return(bnet@nodes[i])), rr)
      dl <- l1 - l2
      sel <- sum(dl[!is.nan(dl)&abs(dl)<Inf])
      if(sel >= 0)
        sel <- clslevs[1]
      else
        sel <- clslevs[2]
      bres10 <- c(bres10, sel == cls[k])
    }
    
    ##########################################################
    ## catnet
    
    catnet.res <- cnSearchOrder(train.data, NULL, npars, nodePars, 0, 
                         nodeOrder=nodeOrder, nodeCats=nodeCats, 
                         parentsPool=NULL, NULL, echo=FALSE)
    bnet0 <- catnet.res@nets[[1]]
    bnet1 <- cnFindBIC(catnet.res)
    bnet2 <- cnFindAIC(catnet.res)

    net1 <- cnSetProb(bnet0, dat1, nodeCats=nodeCats) 
    net2 <- cnSetProb(bnet0, dat2, nodeCats=nodeCats)

    for(k in ind[j:jj]) {
      dd <- matrix(data[,k], ncol=1)
      if(nrow(dd) == length(bnet0@nodes))
        rownames(dd) <- bnet0@nodes
      l1 <- cnNodeLoglik(net1, sapply(1:net1@numnodes, function(i) return(net1@nodes[i])), dd)
      l2 <- cnNodeLoglik(net2, sapply(1:net2@numnodes, function(i) return(net2@nodes[i])), dd)
      dl <- l1 - l2
      sel <- sum(dl[!is.nan(dl)&abs(dl)<Inf])
      if(sel >= 0)
        sel <- clslevs[1]
      else
        sel <- clslevs[2]
      bres0 <- c(bres0, sel == cls[k])
    }
    
    net1 <- cnSetProb(bnet1, dat1, nodeCats=nodeCats) 
    net2 <- cnSetProb(bnet1, dat2, nodeCats=nodeCats)

    for(k in ind[j:jj]) {
      dd <- matrix(data[,k], ncol=1)
      if(nrow(dd) == length(bnet1@nodes))
        rownames(dd) <- bnet1@nodes
      l1 <- cnNodeLoglik(net1, sapply(1:net1@numnodes, function(i) return(net1@nodes[i])), dd)
      l2 <- cnNodeLoglik(net2, sapply(1:net2@numnodes, function(i) return(net2@nodes[i])), dd)
      dl <- l1 - l2
      sel <- sum(dl[!is.nan(dl)&abs(dl)<Inf])
      if(sel >= 0)
        sel <- clslevs[1]
      else
        sel <- clslevs[2]
      bres1 <- c(bres1, sel == cls[k])
    }
    
    net1 <- cnSetProb(bnet2, dat1, nodeCats=nodeCats) 
    net2 <- cnSetProb(bnet2, dat2, nodeCats=nodeCats)

    for(k in ind[j:jj]) {
      dd <- matrix(data[,k], ncol=1)
      if(nrow(dd) == length(bnet2@nodes))
        rownames(dd) <- bnet2@nodes
      l1 <- cnNodeLoglik(net1, sapply(1:net1@numnodes, function(i) return(net1@nodes[i])), dd)
      l2 <- cnNodeLoglik(net2, sapply(1:net2@numnodes, function(i) return(net2@nodes[i])), dd)
      dl <- l1 - l2
      sel <- sum(dl[!is.nan(dl)&abs(dl)<Inf])
      if(sel >= 0)
        sel <- clslevs[1]
      else
        sel <- clslevs[2]
      bres2 <- c(bres2, sel == cls[k])
    }
   
    ##########################################################
    ## svm

    jdata <- cdata[,-ind[j:jj]]
    jdata <- t(jdata)
    jcls <- cls[-ind[j:jj]]
    
    svm.fit <- svm(jcls ~ ., data = data.frame(cbind(jdata,jcls)),
               method = "C-classification", kernel = ker,
               cost = 1, gamma = gam, probability=TRUE)

    for(k in ind[j:jj]) {
      jpred <- predict(svm.fit, t(c(cdata[,k])) , probabilities = TRUE)
      sel <- abs(jpred-clslevs[1])-abs(jpred-clslevs[2])
      if(sel > 0)
        sel <- clslevs[2]
      else
        sel <- clslevs[1]
      bres4 <- c(bres4, sel == cls[k])
    }

    ##########################################################
    ## lars

    x <- t(as.matrix(cdata[,-ind[j:jj]]))
    y <- as.numeric(cls[-ind[j:jj]])
    lars.fit <- lars(x, y, type="lasso", normalize=TRUE, intercept=TRUE)
    ##sc <- -log(lars.fit$RSS) - apply(lars.fit$beta != 0, 1, sum)
    ## Cp RSS/\hat{sigma} - n + 2 *df
    sc <- lars.fit$RSS/lars.fit$RSS[length(lars.fit$RSS)-1] - nrow(x)+2*lars.fit$df
    if(!is.nan(lars.fit$Cp[1]))
      sc <- lars.fit$Cp
    id <- which(sc==min(sc))[1]
    ##cat(id,", ", dim(lars.fit$beta)[1], ",", dim(lars.fit$beta)[2], "\n")
    pred <- predict.lars(lars.fit, t(cdata[, ind[j:jj]]), type="fit",s=id, mode="step")
    predcls <- pred$fit
    predcls[pred$fit<mean(clslevs)] <- clslevs[1]
    predcls[pred$fit>=mean(clslevs)] <- clslevs[2]
    bres8 <- c(bres8, predcls==cls[ind[j:jj]])

    lars.fit <- lars(x, y, type="lar", normalize=TRUE, intercept=TRUE)
    sc <- lars.fit$RSS/lars.fit$RSS[length(lars.fit$RSS)-1] - nrow(x)+2*lars.fit$df
    if(!is.nan(lars.fit$Cp[1]))
      sc <- lars.fit$Cp
    id <- which(sc==min(sc))[1]
    pred <- predict.lars(lars.fit, t(cdata[, ind[j:jj]]), type="fit",s=id, mode="step")
    predcls <- pred$fit
    predcls[pred$fit<mean(clslevs)] <- clslevs[1]
    predcls[pred$fit>=mean(clslevs)] <- clslevs[2]
    bres9 <- c(bres9, predcls==cls[ind[j:jj]])
    
    if(echo)
      cat("[", jj,"\\",ncol(cdata),"] ", "\n")

    j <- jj + 1
  }
  if(echo)
    cat("perf = ", sum(bres1)/length(bres1), "\n")
  t.end <- proc.time()
  eval.time <- as.numeric(t.end[3] - t.start[3])

  return(list(nb=100*sum(bres10)/length(bres10), catnet0=100*sum(bres0)/length(bres0), catnet.bic=100*sum(bres1)/length(bres1), catnet.aic=100*sum(bres2)/length(bres2), svm=100*sum(bres4)/length(bres4), lasso.cp=100*sum(bres8)/length(bres8), lar.cp=100*sum(bres9)/length(bres9), time=eval.time))
 }

