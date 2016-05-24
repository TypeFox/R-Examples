ubRacing <-
function(formula, data, algo, positive=1, ncore = 1, nFold=10, maxFold=10, 
                   maxExp=100, stat.test="friedman", metric="f1", ubConf, verbose=FALSE, ...){
  
  stopifnot(class(formula)=="formula", NROW(data)>1, NCOL(data)>1, is.logical(verbose), ncore>0, nFold>1)
  metric <- match.arg(metric, c("f1","gmean", "auc"))
  stat.test <- match.arg(stat.test, c("friedman","t.bonferroni", "t.holm", "t.none", "no"))
  stopifnot(class(ubConf)=="list", names(ubConf) %in% c("type", "percOver", "percUnder", "k", "perc", "method", "w"))
  
  target <- as.character(formula[[2]])
  tgt <- which(names(data) == target)
  if(length(tgt) == 0)
    stop("target variable not defined")
  
  
  predCandidate <- function(Xtr, Ytr, Xts, Yts, algo, balType, positive, ubConf, metric, verbose, ...){
    
    if (balType != "unbal"){
      #re-balance the dataset
      data <- ubBalance(Xtr, factor(Ytr), type=balType, positive=positive, 
                        ubConf$percOver, ubConf$percUnder, ubConf$k, 
                        ubConf$perc, ubConf$method, ubConf$w, verbose)
      TR <- data.frame(data$X, Y=as.factor(data$Y)) 
    } else 
      TR <- data.frame(Xtr, Y=as.factor(Ytr)) 
    
    #single prediction without parameter tuning
    #library(mlr)
    lrnTask <- makeClassifTask(id=paste(algo, balType, sep="_"), data=TR, target="Y", positive=positive)
    lrnTask <- removeConstantFeatures(lrnTask, show.info=verbose)
    L <- lrnTask$task.desc$class.levels
    negative <- setdiff(L, positive)
    if(length(L) > 2)
      stop("only binary classification supported yet")
    
    race.lrn <- paste(lrnTask$task.desc$type, algo, sep=".")
    lrn <- makeLearner(race.lrn, predict.type = "prob", ...)
    mod <- train(lrn, lrnTask)  
    newdata <- data.frame(Xts, Y=as.factor(Yts)) 
    pred <- predict(mod, newdata=newdata)
    
    perf <- mlr::performance(pred, measures = list(gmean, f1, mlr::auc))
    res.metric <- as.numeric(perf[metric])
    #revert metric since racing is a minimizing algorithm
    res <- 1 - res.metric
    
    res
  }
  
  
  #test each methods on the same observations
  #return a vector of errors (1 or 0 for each method)
  testCandidates <- function(Xtr, Xts, Ytr, Yts, algo, positive, balanceTypes, ubConf, metric, ncore, verbose, ...){
    nBalanceTypes <- length(balanceTypes)
    
    #library(doParallel)
    if (ncore > 1) {  
      doParal <- `%dopar%`
    } else {
      doParal <- `%do%`
    }
    
    j <- NULL 
    
    #We assume that a cluster is already registered!
    #library(foreach)
    #make predictions
    error <- doParal( foreach(j=1:nBalanceTypes, .combine=c, .packages=c('mlr'), 
    						.export=c('predCandidate', 'ubBalance', 'ubUnder', 'ubCNN', 'ubENN', 
  									'ubNCL', 'ubOSS', 'ubOver', 'ubSMOTE', 'ubSmoteExs', 'ubTomek')), 
                      predCandidate(Xtr, Ytr, Xts, Yts, algo, balanceTypes[j], positive, ubConf, metric, verbose, ...))
    
    #   error <- NULL
    #   for(j in 1:nBalanceTypes) {  
    #     error <- c(error, predCandidate(Xtr, Ytr, algo, balanceTypes[j], positive, ubConf, metric, verbose, ...))
    #   }
    

    names(error) <- balanceTypes
    return(error)
  }
  
  
  
 
  
  
  #from race package but edited
  aux2.friedman.edit<- function(y, I = 1:ncol(y), n=nrow(y), conf.level = 0.95,interactive=F) {
    k <- length(I)
    r <- t(apply(y[1:n, I, drop=FALSE], 1, rank))
    A <- sum(as.vector(r)^2)
    R <- apply(r, 2, sum)
    J <- I[order(R)]
    alpha <- 1 - conf.level
    TIES <- tapply(r, row(r), table)
    STATISTIC <- ((12 * sum((R - n * (k + 1)/2)^2))/(n * 
                                                       k * (k + 1) - (sum(unlist(lapply(TIES, function(u) {
                                                         u^3 - u
                                                       })))/(k - 1))))
    PARAMETER <- k - 1
    PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
    if (!is.nan(PVAL) && (PVAL < alpha)) {
      if (interactive) 
        cat("|-|")
      t <- qt(1 - alpha/2, (n - 1) * (k - 1)) * (2 * (n * 
                                                        A - sum(R^2))/((n - 1) * (k - 1)))^(1/2)
      o <- order(R)
      J <- I[o[1]]
      for (j in 2:k) if (abs(R[o[j]] - R[o[1]]) > t) 
        break
      else J <- c(J, I[o[j]])
    }
    else {
      if (interactive) 
        cat("|=|")
    }
    return(J)
  }
  
  #from race package but edited    
  aux.friedman.edit<- function(Results,no.subtasks.sofar,alive,conf.level=0.95,interactive=F) {
    which.alive <- which(alive)
    no.alive<-length(which.alive)
    
    if (no.alive == 2) {
      V1 <- Results[1:(no.subtasks.sofar), which.alive[1]]
      V2 <- Results[1:(no.subtasks.sofar), which.alive[2]]
      PVAL <- wilcox.test(V1, V2, paired = TRUE, exact = FALSE)$p.value
      D <- V1 - V2
      w.stat <- sum(sign(D) * rank(abs(D)))
      if (!is.nan(PVAL) && !is.na(PVAL) && (PVAL < 1 - conf.level)) {
        if (interactive) 
          cat("|-|")
        if (w.stat < 0) {
          best <- which.alive[1]
          alive[which.alive[2]] <- FALSE
        }
        else {
          best <- which.alive[2]
          alive[which.alive[1]] <- FALSE
        }
      }
      else {
        if (interactive) 
          cat("|=|")
        if (w.stat < 0) {
          best <- which.alive[1]
        }
        else {
          best <- which.alive[2]
        }
      }
    }
    else {
      J <- aux2.friedman.edit(Results[1:(no.subtasks.sofar), ,drop=FALSE], I=which.alive, conf.level=conf.level, interactive=interactive)
      alive[-J] <- FALSE
      best <- J[1]
    }
    
    return(list(alive=alive,best=best))
  }
  
  #from race package but edited
  aux.ttest.edit <- function(Results,no.subtasks.sofar,alive,adjust=c("none", "bonferroni", "holm"),conf.level=0.95,interactive=F) {
    which.alive <- which(alive)
    no.alive<-length(which.alive)
    
    adjust <- match.arg(adjust)
    mean.all <- array(0, c(ncol(Results)))
    for (j in 1:ncol(Results)) mean.all[j] <- sum(Results[1:no.subtasks.sofar, 
                                                          j]/no.subtasks.sofar)
    best <- match(min(mean.all[alive]), mean.all)
    PJ <- array(0, dim = c(2, 0))
    for (j in which.alive) {
      Vb <- Results[1:no.subtasks.sofar, best]
      Vj <- Results[1:no.subtasks.sofar, j]
      p <- t.test(Vb, Vj, paired = TRUE)$p.value
      if (!is.nan(p) & !is.na(p)) 
        PJ <- array(c(PJ, j, p), dim = dim(PJ) + c(0, 
                                                   1))
    }
    PJ[2, ] <- p.adjust(PJ[2, ], method = adjust)
    dropped.any <- FALSE
    for (j in 1:ncol(PJ)) if (PJ[2, j] < (1 - conf.level)) {
      alive[PJ[1, j]]<- FALSE
      dropped.any <- TRUE
    }
    if (interactive) {
      if (dropped.any) 
        cat("|-|")
      else cat("|=|")
    }
    return(list(alive=alive,best=best))
  }
  
  #from race package but edited
  format.precis <- function(title, value,title.width,value.width) {
    dots <- function(n) return(paste(rep(".", n), collapse = ""))
    spaces <- function(n) return(paste(rep(" ", n), collapse = ""))
    string <- paste(title, dots(title.width - nchar(title)), sep = "")
    if (nchar(value) <= value.width) {
      string <- paste(string, dots(value.width - nchar(value)),value, sep = "")
    }
    else {
      f.vec <- strwrap(value, width = value.width)
      first <- f.vec[1]
      first <- paste(dots(title.width - nchar(first)), 
                     first, sep = "")
      rest <- format(f.vec[-1])
      rest <- paste(spaces(title.width + value.width - max(nchar(rest))), rest, sep = "", collapse = "\n")
      string <- paste(string, paste(first, rest, sep = "\n"), sep = "")
    }
    return(paste(string, "\n"))
  }
  
  
  #from function createFolds in caret package but edited
  idFolds <- function(y, k = 10){
    
    stopifnot(k > 0)
    
    if (is.numeric(y)) {
      cuts <- floor(length(y)/k)
      if (cuts < 2) 
        cuts <- 2
      if (cuts > 5) 
        cuts <- 5
      y <- cut(y, unique(quantile(y, probs = seq(0, 1, length = cuts))), include.lowest = TRUE)
    }
    if (k < length(y)) {
      y <- factor(as.character(y))
      numInClass <- table(y)
      foldVector <- vector(mode = "integer", length(y))
      for (i in 1:length(numInClass)) {
        seqVector <- rep(1:k, numInClass[i]%/%k)
        if (numInClass[i]%%k > 0) 
          seqVector <- c(seqVector, sample(1:k, numInClass[i]%%k))
        foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
      }
    } else 
      foldVector <- seq(along = y)
    
    
    foldVector
  }
  
  foldVector <- idFolds(data[ ,tgt], nFold)
  input <- data[ ,-tgt]
  output <-  data[ ,tgt]
  
  balanceTypes <- c("unbal", "ubOver", "ubUnder", "ubSMOTE", "ubOSS","ubCNN", "ubENN", "ubNCL", "ubTomek")
  nBalanceTypes <- length(balanceTypes)
  title <- paste("Racing for unbalanced methods selection in", nFold, "fold CV")
  title.width<-32
  value.width<-title.width
  precis <- paste("\n", title,"\n",
                  format.precis("Number of candidates", nBalanceTypes, title.width, value.width), 
                  format.precis("Max number of folds in the CV", maxFold,title.width,value.width), 
                  format.precis("Max number of experiments", ifelse(maxExp, maxExp, "unlimited"),title.width,value.width), 
                  #format.precis("Metric", metric, title.width, value.width), 
                  format.precis("Statistical test", 
                                switch(stat.test, 
                                       friedman = "Friedman test", 
                                       t.bonferroni = "t-test with Bonferroni's correction for multiple comparisons", 
                                       t.holm = "t-test with Holm's correction for multiple comparisons", 
                                       t.none = "t-test with no correction for multiple comparisons",
                                       no = "No test"),title.width,value.width)
  )
  cat(precis) 
  if (nBalanceTypes > maxExp & maxExp!=0) 
    stop("Max number of experiments is smaller than number of candidates")
  if (maxFold > nFold)
    maxFold <- nFold 
  
  Results <- matrix(NA, nrow=nFold, ncol=nBalanceTypes)
  colnames(Results) <- balanceTypes
  alive <- rep(TRUE,nBalanceTypes)
  n.Exp <- 0
  n.alive <- nBalanceTypes
  types.alive <- balanceTypes
  
  #open a cluster
  if (ncore > 1) {  
    cl <- makeCluster(ncore)
    cat(" Parallel execution with", ncore, "cores \n")
    registerDoParallel(cl)  
  } else
	makeCluster <- NULL
  
  cat("                                                               \n",
      "                            Markers:                           \n",
      "                               x No test is performed.         \n",
      "                               - The test is performed and     \n",
      "                                 some candidates are discarded.\n", 
      "                               = The test is performed but     \n",
      "                                 no candidate is discarded.    \n",
      "                                                               \n",
      "+-+-----------+-----------+-----------+-----------+-----------+\n", 
      "| |       Fold|      Alive|       Best|  Mean best| Exp so far|\n", 
      "+-+-----------+-----------+-----------+-----------+-----------+\n")
  
  for(i in 1:nFold){
    n.Exp <- n.Exp + n.alive
    
    if ((n.Exp>maxExp | maxFold<i) & maxExp!=0 | n.alive==1) 
      break()
    
    max.tasks <- i
    id.fold <- which(foldVector == i)
    
    Xtr <- input[foldVector == i, ]
    Xts <- input[foldVector != i, ]
    Ytr <- output[foldVector == i]
    Yts <- output[foldVector != i]
    
    if(any(is.na(types.alive)))
      stop("NAs in balanceTypes")
    
    err <- testCandidates(Xtr, Xts, Ytr, Yts, algo, positive, types.alive, ubConf, metric, ncore, verbose, ...)
    Results[i, match(types.alive, balanceTypes)] <- err
    
    #the best candidate is the one that obtains the smallest results on the various tasks considered
    if (i==1)
      best <- as.numeric(which(err==min(err)))[1]
    #make a test on the error
    if (i>1){	
      if (stat.test=="friedman")
        test <- aux.friedman.edit(Results, no.subtasks.sofar=i, alive, conf.level=0.95, interactive=TRUE) 
      if (stat.test=="t.none")
        test<-aux.ttest.edit(Results,no.subtasks.sofar=i,alive,"none", conf.level=0.95)
      if (stat.test=="t.holm")
        test<-aux.ttest.edit(Results,no.subtasks.sofar=i,alive,"holm", conf.level=0.95)
      if (stat.test=="t.bonferroni")
        test<-aux.ttest.edit(Results,no.subtasks.sofar=i,alive,"bonferroni", conf.level=0.95)
      if (stat.test=="no"){
        test <- list(alive=rep(TRUE, nBalanceTypes), best=as.numeric(which(err==min(err)))[1])
        cat("|x|")
      }
      
      alive <- test$alive
      best <- test$best
    }
    #remove significantly worse mehtods
    which.alive <- which(alive)
    n.alive <- length(which.alive)
    types.alive <- balanceTypes[which.alive]	
    if(i==1){
      avg.best <- min(err)
      cat(" |x|")
    }
    else 
      avg.best <- as.numeric(apply(Results, 2, mean, na.rm=TRUE)[best])
    
    #revert metric since racing is a minimizing algorithm
    avg.best <- 1 - avg.best
    
    cat(paste(formatC(i, width = 11), "|", 
              formatC(n.alive, width = 11), "|", 
              formatC(best, width = 11), "|", 
              formatC(avg.best, width = 11), "|",
              formatC(n.Exp, width = 11), "|\n"," ",
              sep = ""))
    
  }
  
  cat("+-+-----------+-----------+-----------+-----------+-----------+\n")
  colnames(Results) <- balanceTypes
  Results <- Results[1:max.tasks, ]
  #revert metric since racing is a minimizing algorithm
  Results <- 1 - Results
  #write.csv(Results,file="Race.csv")
  
  #close the cluster
  if (ncore > 1) {
    stopCluster(cl)
    rm(cl)
  } else
	stopCluster <- NULL
  
  if(max.tasks>1){
    avg.err <- apply(Results, 2, mean, na.rm=TRUE)
    sd.err <- apply(Results, 2, sd, na.rm=TRUE)
  }
  else {	#max.tasks==1
    avg.err <- Results
    sd.err <- Results
  }
  avg.best <- as.numeric(avg.err[best])
  sd.best <- as.numeric(sd.err[best])
  best.desc <- balanceTypes[best]
  #percentage of computational gain
  percGain <- round((1 - n.Exp/(nFold*nBalanceTypes))*100, 0)
  
  cat("Selected candidate:",best.desc,"\t metric:", metric, "\t mean value:", round(avg.best, 4),"\n")
  return(list(best=best.desc, avg=avg.best, sd=sd.best, N.test=n.Exp, Gain=percGain, Race=Results))
  
}
