
sdnConsolidate <- function(train, test) {
  
  train.pathway <- rownames(train$cdata)
  test.pathway <- rownames(test$cdata)

  rgen <- NULL
  for(gen in 1:length(train.pathway))
    if(length(which(test.pathway == train.pathway[gen]))==1 && length(which(train.pathway == train.pathway[gen]))==1)
      rgen <- c(rgen, gen)
  if(length(rgen) < 1) {
    warning("no genes found")
    return(NULL)
  }
  train <- as.matrix(train[rgen,], nrow=length(rgen))
  
  train.pathway <- rownames(train$cdata)
  rgen <- NULL
  for(gen in 1:length(train.pathway)) {
    id <- which(test.pathway==train.pathway[gen])
    if(length(id)>=1)
      rgen <- c(rgen, id[1])
  }
  if(length(rgen) < 1) {
    warning("no genes found")
    return(NULL)
  }
  test <- as.matrix(test[rgen,], nrow=length(rgen))

  return(list(train, test))
}


## train and test have the following components 
## geneset: gene subset
## cdata: expression data as gene-row matrix
## cls:: class labels 

sdnLearn <- function(data, cls, clslevs = NULL, ncats = 3, nodeCats = NULL, quant="uniform", std=TRUE) {

  if(!is.numeric(data)||!is.matrix(data))
    stop("Numeric matrix expected")
  if(nrow(data)*ncol(data) <= 0)
    stop("Insufficient data")
  if(length(cls) != ncol(data))
    stop("Labels should match data columns")
  geneset <- rownames(data)
  if(is.null(geneset))
    geneset <- 1:nrow(data)
  
  if(is.numeric(cls))	
    cls <- as.factor(as.integer(cls))
  if(is.null(clslevs))
    clslevs <- levels(as.factor(cls))
  for(cl in cls)
    if(!cl%in%clslevs)
      stop("Wrong class labels")
  if(sum(is.na(clslevs))>0)
    stop("Wrong train$cls")
  ncls <- length(clslevs)
                 
  N <- nrow(data)
  if(is.null(nodeCats))
    nodeCats <- lapply(1:N, function(nn) c(1:ncats))
  names(nodeCats) <- geneset

  n1 <- sum(cls==clslevs[1])
  n2 <- sum(cls==clslevs[2])
  dprior <- log(n1) - log(n2)
  
  if(std) {
    ## genes z-transformation
    data <- data - apply(data, 1, "mean")
    vv <- apply(data, 1, "var")
    vv[is.na(vv)] <- 1
    data <- data/sqrt(vv/ncol(data))
  }
  
  ## soft discretization
  qdata <- sdnet::cnDiscretize(data, ncats, mode="soft", marginal=quant, learnset=1:ncol(data), cover=0.95)
  pdata <- qdata$pdata
  rm(qdata)

  bnet <- sdnet::cnNew(geneset, cats=nodeCats, pars=vector("list",N), probs=NULL, dagonly=TRUE)

  nets <- lapply(clslevs, function(cl) {
    dat <- pdata[, cls==cl]
    ## sets P(X=k) \propto \sum_s q_k(y^s) = P(y^s|X=k)/(\sum_m P(y^s|X=m))
    return(sdnet::cnSetProb(bnet, dat, nodeCats=nodeCats, softmode=TRUE))
  })

  res <- list(geneset=geneset, clslevs=clslevs, nets=nets, nodeCats=nodeCats, quant=quant, dprior=dprior)
  return(res)
}

sdnPredict <- function(model, data, std=TRUE) {
  
  if(nrow(data)*ncol(data) <= 0)
    stop("Insufficient data")
  if(is.null(model$geneset) || is.null(model$nodeCats))
    stop("Wrong model")
  geneset <- model$geneset
  numnodes <- length(geneset)
  nets <- model$nets
  for(net in nets) {
    if(net@numnodes != numnodes)
      stop("Wrong model's nets")
  }
      
  rgen <- NULL
  testgenes <- rownames(data)
  if(is.null(testgenes))
    testgenes <- 1:nrow(data)
  for(gen in 1:length(geneset)) {
    id <- which(testgenes==geneset[gen])
    if(length(id)>=1)
      rgen <- c(rgen, id[1])
  }
  if(length(rgen) != length(geneset)) {
    warning("Incompatible model and test data")
    return(NULL)
  }
  testgeneset <- geneset[rgen]
  data <- data[rgen,]
    
  clslevs <- model$clslevs
  if(sum(is.na(clslevs))>0)
    stop("Wrong clslevs")
  
  ncats <- sapply(model$nodeCats, function(cc) length(cc))

  if(std) {
    ## genes z-transformation
    data <- data - apply(data, 1, "mean")
    vv <- apply(data, 1, "var")
    vv[is.na(vv)|vv==0] <- 1
    data <- data/sqrt(vv/ncol(data))
  }
  
  bres1 <- NULL
  bres1cls <- NULL
  
  ## soft discretization
  qdata <- sdnet::cnDiscretize(data, ncats, mode="soft", marginal=model$quant, learnset=1:ncol(data), cover=0.95)
  pmdata <- qdata$pdata
  rm(qdata)
  
  pl <- NULL
  ind <- NULL
  for(k in 1:ncol(pmdata)) {
    dd <- t(matrix(pmdata[,k], nrow=ncats))
    ##r <- sapply(nets, function(net) sapply(1:net@numnodes, function(j) log(sum(dd[j,]*net@probs[[j]]))))
    r <- sapply(nets, function(net) sapply(1:net@numnodes, function(j) sum(net@probs[[j]]*log(dd[j,]/net@probs[[j]]))))
    ind <- apply(r, 1, function(rr) sum(abs(rr)>=Inf | is.nan(rr) | is.na(rr)))
    ind <- ind <= 0
    r <- r[ind,]
    r <- apply(r, 2, sum)
    pl <- cbind(pl, r)
  }
  rm(pmdata)
  rm(nets)

  colnames(pl) <- colnames(data)
  rownames(pl) <- clslevs
  if(nrow(pl) == 2) {
    pred <- pl[2,]-pl[1,]
    return(pred)
  }
  return(pl)
}


auroc <- function(pred, tres){
  if(is.na(sum(tres)) || length(pred) != length(tres))
    return(NA)
  n <- length(pred)
  ind <- order(pred, decreasing=TRUE)
  
  tp <- NULL
  tn <- NULL
  fp <- NULL
  fn <- NULL
  for(i in ind) {
    x <- pred[i]
    tp <- c(tp, sum(pred >= x & tres))
    fp <- c(fp, sum(pred >= x & !tres))
    tn <- c(tn, sum(pred < x & !tres))
    fn <- c(fn, sum(pred < x & tres))
  }
  fp <- c(0,fp)
  tp <- c(0,tp)
  tn <- c(sum(!tres),tn)
  fn <- c(sum(tres),fn)
  spec <- tn/(tn+fp)
  rec <- tp/(tp+fn)
  ff <- sum((rec[-1]-rec[1:(length(rec)-1)])*(spec[-1]+spec[1:(length(spec)-1)])/2)
  return(ff)
}

sdnEvaluate <- function(train, test, ncats = 3, nodeCats = NULL, std=FALSE) {

  bseparatedisc <- FALSE
  quant <- "quantile"

  if(nrow(train$cdata)*ncol(train$cdata) <= 0 || nrow(test$cdata)*ncol(test$cdata) <= 0)
        stop("Insufficient data")

  if(is.numeric(train$cls))	
    train$cls <- as.factor(as.integer(train$cls))
  if(!is.factor(train$cls))
    stop("train$cls should be factor")
  
  clslevs <- levels(train$cls)
  if(sum(is.na(clslevs))>0)
    stop("Wrong train$cls")
  if(length(clslevs) != 2)
    stop("sdnEvaluate handles 2-class problems only")
  
  if(is.null(train$geneset))
     train$geneset <- rownames(train$cdata)
  if(is.null(train$geneset))
    train$geneset <- 1:nrow(train$cdata)
  if(is.null(test$geneset))
     test$geneset <- rownames(test$cdata)
  if(is.null(test$geneset))
    test$geneset <- 1:nrow(test$cdata)
  
  rgen <- NULL
  for(gen in 1:length(train$geneset))
    if(length(which(test$geneset == train$geneset[gen]))==1)
      rgen <- c(rgen, gen)
  if(length(rgen) < 1) {
    warning("No common genes found")
    return(NULL)
  }
  train$cdata <- train$cdata[rgen,]
  nodePars <- vector("list",length(rgen))
  for(i in 1:length(rgen))
    if(!is.null(train$nodePars[[rgen[i]]]))
      nodePars[[i]] <- train$nodePars[[rgen[i]]]
  train$nodePars <- nodePars
  train$geneset <- train$geneset[rgen]
  
  rgen <- NULL
  for(gen in 1:length(train$geneset)) {
    id <- which(test$geneset==train$geneset[gen])
    if(length(id)>=1)
      rgen <- c(rgen, id[1])
  }
  if(length(rgen) < 1) {
    warning("No common genes found")
    return(NULL)
  }
  nodePars <- vector("list",length(rgen))
  for(i in 1:length(rgen))
    if(!is.null(test$nodePars[[rgen[i]]]))
      nodePars[[i]] <- test$nodePars[[rgen[i]]]
  test$nodePars <- nodePars
  test$geneset <- test$geneset[rgen]
  test$cdata <- test$cdata[rgen,]

  geneset <- train$geneset
  numnodes <- length(geneset)

  if(is.null(nodeCats))
    nodeCats <- lapply(1:numnodes, function(nn) c(1:ncats))
  names(nodeCats) <- geneset
  
  if(std) {
    ## genes z-transformation
    train$cdata <- train$cdata - apply(train$cdata, 1, "mean")
    vv <- apply(train$cdata, 1, "var")
    vv[is.na(vv)] <- 1
    train$cdata <- train$cdata/sqrt(vv/ncol(train$cdata))
    test$cdata <- test$cdata - apply(test$cdata, 1, "mean")
    vv <- apply(test$cdata, 1, "var")
    vv[is.na(vv)|vv==0] <- 1
    test$cdata <- test$cdata/sqrt(vv/ncol(test$cdata))
  }

  N <- numnodes
  n1 <- sum(train$cls==clslevs[1])
  n2 <- sum(train$cls==clslevs[2])
  dprior <- log(n1) - log(n2)
  
  ## soft discretization
  if(bseparatedisc) {
    qdata <- sdnet::cnDiscretize(train$cdata, ncats, mode="soft", marginal=quant, learnset=1:ncol(train$cdata), cover=0.95)
    pdata <- qdata$pdata
    rm(qdata)
    qdata <- sdnet::cnDiscretize(test$cdata, ncats, mode="soft", marginal=quant, learnset=1:ncol(test$cdata), cover=0.95)
    pmdata <- qdata$pdata
    rm(qdata)
  }
  else { 
    qdata <- sdnet::cnDiscretize(cbind(train$cdata,test$cdata), ncats, mode="soft", marginal=quant, learnset=1:ncol(train$cdata), cover=0.95)
    pmdata <- qdata$pdata[,(ncol(train$cdata)+1):(ncol(train$cdata)+ncol(test$cdata))]
    pdata <- qdata$pdata[,1:ncol(train$cdata)]
    rm(qdata)
  }

  dat1 <- pdata[, train$cls==clslevs[1]]
  dat2 <- pdata[, train$cls==clslevs[2]]
  
  bnet <- sdnet::cnNew(geneset, cats=nodeCats, pars=vector("list",numnodes), probs=NULL, dagonly=TRUE)
 
  ## sets P(X=k) \propto \sum_s q_k(y^s) = P(y^s|X=k)/(\sum_m P(y^s|X=m))
  bnet <- sdnet::cnSetProb(bnet, pdata, nodeCats=nodeCats[bnet@nodes], softmode=TRUE)
  net1 <- sdnet::cnSetProb(bnet, dat1, nodeCats=nodeCats, softmode=TRUE)
  net2 <- sdnet::cnSetProb(bnet, dat2, nodeCats=nodeCats, softmode=TRUE)
  
  pl1 <- -sapply(1:net1@numnodes, function(j) sum(bnet@probs[[j]]*log(net1@probs[[j]]/bnet@probs[[j]])))
  pl1 <- 2*(n1+n2)*(n1/n2)*pl1
  ind <- order(pl1, decreasing = TRUE)
  hc <- sapply(1:length(pl1), function(k) sum(pl1[ind[1:k]]) - 0.5*log(n1+n2)*2*k)
  kmax <- which(hc==max(hc[!is.nan(hc)]))[1]
  hc <- hc[1:kmax]
  if(is.na(kmax) || kmax<2) kmax <- 2
  ind <- ind[1:kmax]
  rm(bnet)
  
  bres1 <- NULL
  bres1cls <- NULL
  tres <- test$cls
  
  for(k in 1:ncol(pmdata)) {
    dd <- t(matrix(pmdata[,k], nrow=ncats))
    pl1 <- sapply(ind, function(j) sum(net1@probs[[j]]*log(dd[j,]/net1@probs[[j]])))
    pl2 <- sapply(ind, function(j) sum(net2@probs[[j]]*log(dd[j,]/net2@probs[[j]])))
    pdl <- pl2-pl1
    pdl <- sum(pdl[!is.nan(pdl)&abs(pdl)<Inf])
    bres1 <- c(bres1, pdl)
    sel <- clslevs[1]
    if(pdl > 0)
      sel <- clslevs[2]
    bres1cls <- c(bres1cls, sel)
  }
  rm(pmdata)
  rm(net1)
  rm(net2)

  ##dp <- quantile(bres1, sum(tres==1)/length(tres))
  ##bres1cls <- bres1
  ##bres1cls[bres1 < dp] <- clslevs[1]
  ##bres1cls[bres1 >= dp] <- clslevs[2]
  names(bres1cls) <- colnames(test$cdata)
  auc1 <- auroc(bres1, tres)

  acc1 <- 0.5*(sum(bres1cls==clslevs[1]&tres==clslevs[1])/sum(tres==clslevs[1]) + sum(bres1cls==clslevs[2]&tres==clslevs[2])/sum(tres==clslevs[2]))
    
  res <- c(acc1, auc1)
  names(res) <- c("acc", "auc")
  return(list(res, bres1, bres1cls, geneset=geneset[ind]))
}
