
fct_interselect <- function(res, indmain, m=100, mat, cutoff = 0){
      vim <-  sapply(res, function(arg) return(arg)) 
      
      if (class(cutoff) == "function") {
        cutoff <- do.call(cutoff, list(vim))
      }
      
      Y <- apply(vim, 2, function(arg){
                        arg > cutoff})
      Z <- matrix(0, ncol = nrow(Y), nrow = nrow(Y))
      for (w in 1:ncol(vim)){
        for (i in 2:ncol(mat)){
          for (j in 1:(i-1)){
              if (Y[i,w] == 1 & Y[j,w] == 1) Z[i,j] <- Z[i,j] +1
            }
          }
      }
      lpaare <- m 
      hfgktpaare <- cbind(which(Z >= 1, arr.ind = TRUE),
                          Z[which(Z >= 1, arr.ind = TRUE)])[order(Z[which(Z >= 1, arr.ind = TRUE)],decreasing=TRUE),]
      paare <- which(Z >= max(min(sort(Z, decreasing = T)[1:min(length(Z),lpaare)]),1), arr.ind = TRUE)
      xcombmat <- apply(paare,1, function(arg){mat[,arg[1]]*mat[,arg[2]]})
      if(length(paare)!=0) {colnames(xcombmat) <-  apply(paare,1, function(arg){paste(colnames(mat)[arg[2]],':',colnames(mat)[arg[1]], sep = '')})}
      out <- list()
      out$xcombmat0 <- xcombmat
      out$hfgktpaare <- hfgktpaare
      return(out)
 }



fct_indexlist <- function(set = 123, r = 200, x){ 
  # Indizes Subsamples 1  
  indexinb <- indexoob <- list()                                  
  set.seed(set)                                       
  seed1 <- as.list(sample(1:1000000, r, replace = FALSE)) 
  for(i in 1:length(unlist(seed1))){                               
    set.seed(seed1[[i]])                                         
    ind <-  sample(1:nrow(x), length(1:nrow(x))*0.632)  
    indexinb[[i]] <- ind                                
    indexoob[[i]] <- c(1:nrow(x))[-ind]                               
  }           
  
  out <- list()                                                   
  out$indexinb <- indexinb                                        
  out$indexoob <- indexoob       
  return(out)                                                     
}          


resample.sprinter <- function(x, 
                              time, 
                              status,
                              fold = 5,
                              oob.rel = 0.632,
                              mandatory,
                              repetitions = 25, 
                              n.inter.candidates = 1000, 
                              screen.main, 
                              screen.inter = fit.rf, 
                              fit.final,
                              args.screen.main = list(), 
                              args.screen.inter = list(),
                              args.fit.final = args.screen.main, 
                              orthogonalize = TRUE,
                              parallel = FALSE, 
                              mc.cores = detectCores(), ...){
  
  
  seed <- sample(1:1000000, fold, replace = FALSE)  
  set.seed(seed+100)
  cv <- list()
  for (i in 1:fold){
    cv[[i]] <- sample(1:length(time), length(time)*oob.rel)
  }
  
   
  
  out <- lapply(as.list(1:fold), 
                    function(arg){
                      sprinter(x = x[cv[[arg]],],
                               time = time[cv[[arg]]],
                               status = status[cv[[arg]]],
                               mandatory = mandatory,
                               repetitions = repetitions,
                               n.inter.candidates = n.inter.candidates,
                               screen.main = screen.main, 
                               screen.inter = screen.inter, 
                               fit.final = fit.final,
                               args.screen.main = args.screen.main, 
                               args.screen.inter = args.screen.inter,
                               args.fit.final = args.fit.final, 
                               orthogonalize = orthogonalize,
                               parallel = parallel,
                               mc.cores = mc.cores,...
                               )
                      }
                    )  
  class(out) <- "resample.sprinter"
  return(out)  
}
  
   


 sprinter <- function(x, 
                      time, 
                      status= rep(0, nrow(x)),
                      mandatory= NULL,
                      repetitions = 25, 
                      n.inter.candidates =1000, 
                      screen.main, 
                      screen.inter=fit.rf,
                      fit.final = screen.main, 
                      args.screen.main = list(), 
                      args.screen.inter = list(),
                      args.fit.final = args.screen.main, 
                      orthogonalize = TRUE,
                      cutoff = 0,
                      parallel = FALSE, mc.cores = detectCores(), ...){   
    # Einstellungen:
      xmat <- as.matrix(x)
      xdf <- as.data.frame(x)
              
    # Berechnen der Haupteffekte mit Hilfe von Cb1:
      args.screen.main$time <- time
      args.screen.main$status <- status
      args.screen.main$x <- xmat
      if (!is.null(mandatory)){ 
        unpen.index <- which(colnames(x) %in% mandatory)
        args.screen.main$unpen.index <- unpen.index
      }
      
      cat('Screening for main effects \n')
      res1 <- do.call(screen.main,args.screen.main)
  
    # Berechnen der Residuen:
      inds1 <- c(1:ncol(xmat))[-res1$indmain]
      haupt1 <- res1$indmain
      
      if (length(haupt1 > 0) & orthogonalize == TRUE) {
        xmats <- apply(
          xmat[,inds1], 
          2, 
          function(arg){
            predict(glm(arg ~ xmat[,haupt1]))
            })
        colnames(xmats) <- colnames(xmat)[inds1]
        zmat <- xmat[,inds1] - xmats
        zdf  <- as.data.frame(
                  cbind(
                    xmat[,haupt1] ,zmat, status, time
                    )
                  )
        colnames(zdf)[1:length(haupt1)] <- colnames(xmat)[haupt1]
      } else { 
        zdf <- as.data.frame(cbind(xdf, status, time))
        colnames(zdf)[(ncol(zdf)-1):ncol(zdf)] <- c('status', 'time')
      }

    # Erstellen der Indizes fuer innere Subsamples (repetitions):
      isub <- list()
      length(isub) <- repetitions
      for (i in 1:repetitions){
        isub[[i]]   <- sample(1:length(time), round(length(time)*0.632) , replace = FALSE)
      }
    
    # Durchfuehren der inneren Subsamples (50):
      seed.inter <- sample(1:1000000,repetitions)
      cat('Screening for interactions: \n')
      cat('Resample run: ')
      if (parallel){
        if(require(parallel))
          
        res <- mclapply(1:repetitions, 
                          function(arg){
                            args.screen.inter <- c(args.screen.inter,list(nr = arg, data = zdf, indices = isub[[arg]], seed.interselect = seed.inter[arg]))
                            do.call(screen.inter,args.screen.inter)
                          }, mc.cores = mc.cores
                       )
      } else {
        res <- lapply(1:repetitions, 
                        function(arg){
                          args.screen.inter <- c(args.screen.inter,list(nr = arg, data = zdf, indices = isub[[arg]], seed.interselect = seed.inter[arg]))
                          do.call(screen.inter,args.screen.inter)
                        }
                      )
      }

   
    # Erstellen einer Matrix mit Interaktionen:
      interselectout <- fct_interselect(res, indmain=haupt1, mat = xmat, m=n.inter.candidates, cutoff = cutoff)
      xcombmat0 <- interselectout$xcombmat

      if (length(haupt1) >  1){
      hauptpaare <- combn(haupt1,2)
      xcombmathe <- matrix(NA, nrow(xmat), ncol(hauptpaare))
              for(j in 1:ncol(hauptpaare)){
                  xcombmathe[,j] <-  xmat[,hauptpaare[1,j]]*xmat[,hauptpaare[2,j]]
                  }
      colnames(xcombmathe) <-  apply(hauptpaare, 2, function(arg){paste(colnames(xmat)[arg[1]],':',colnames(xmat)[arg[2]], sep = '')})
      }

    # Merge Daten
      if (length(haupt1) > 1){
        xcombmat <- cbind(xmat[,haupt1], xcombmathe, xcombmat0)
      } else {
        xcombmat <- cbind(xmat[,haupt1], xcombmat0)
      }
      if (ncol(xcombmat)>1) {
        xcombmat <- xcombmat[,!duplicated(colnames(xcombmat))]
      }

    
    # New Unpen.index:
      if(!is.null(mandatory)){
        unpen.index2 <- rep(NA, length(unpen.index))
        for (i in 1:length(unpen.index)){ unpen.index2[i] <- which(colnames(xcombmat)== mandatory[i])}
      } else {unpen.index2 <- NULL}
    

    # Berechne ein CoxBoost CB2 mit Haupteffekten aus CB1 und den Interaktionen
                      # ueber die Variablen mit den hoechsten IFs:
      args.fit.final$time <- time
      args.fit.final$status <- status
      args.fit.final$x <- xcombmat
      if(!is.null(mandatory)){ args.fit.final$unpen.index <- unpen.index2}
    
    if (sum(xcombmat) == 0){
      cat('No Main effects and Interactions are preselected')
      
    }
    
      cat('\n Fitting joint model \n')
      res2 <- do.call(fit.final,args.fit.final)
      

    # Output:
      out <- list()
      out$call <- match.call()
      out$xnames <- colnames(x)
      out$mandatory <- mandatory
      out$repetitions <- repetitions
      out$n.inter.candidates <- n.inter.candidates
      if(is.matrix(interselectout$hfgktpaare) & sum(interselectout$hfgktpaare)!=0) {
        out$inter.candidates <- interselectout$hfgktpaare[1:min(100,nrow(interselectout$hfgktpaare)),]
      } else {
        out$inter.candidates <- interselectout$hfgktpaare
      }
      out$main.model <- res1
      out$final.model <- res2
      class(out) <- c("sprinter",class(out))
      return(out)
  }


  
 fit.CoxBoost <- function(time, 
                    status, 
                    x, 
                    unpen.index=NULL, 
                    seed = 123, 
                    stepno = NULL,
                    K = 10, 
                    criterion = 'pscore', 
                    nu = 0.05,
                    maxstepno=200,
                    standardize = T, 
                    trace = T,... 
                    ){
    penalty <- sum(status)*(1/nu-1) 
    if (is.null(stepno)){
      set.seed(seed)
      cv.subsample <- cv.CoxBoost(time=time,
                                   status=status,
                                   x=x,
                                   unpen.index = unpen.index,
                                    maxstepno=200,K=K,penalty=penalty,
                                   criterion = criterion,...)
      stepno <- cv.subsample$optimal.step
      } 
      set.seed(seed)
      cb <- CoxBoost(time=time,
                      status=status,
                      x=x,
                      unpen.index = unpen.index, 
                      stepno=stepno,penalty=penalty,
                      criterion = criterion,...)
      indmain <- which(cb$coefficients[nrow(cb$coefficients),] != 0)
      res <- list()
      res$model <- cb
      res$indmain <- indmain
      res$xnames <- colnames(x)[indmain]
      res$beta <- cb$coefficients[nrow(cb$coefficients),][which(cb$coefficients[nrow(cb$coefficients),] != 0) ]
      return(res)
    }
     
fit.GAMBoost <- function(time, 
                         status, 
                         x, 
                         unpen.index=NULL, 
                         seed = 123, 
                         stepno = NULL,
                         penalty = 100,
                         maxstepno=200,
                         standardize = T, 
                         criterion = 'deviance',
                         family = gaussian(),
                         trace = T,... 
){
  penalty <- rep(penalty, ncol(x))
  if (!is.null(unpen.index)){penalty[unpen.index] <- 0}
  
  if (is.null(stepno)){
    set.seed(seed)
    gamb <- cv.GAMBoost(y=time,
                            x.linear = x,
                            maxstepno = maxstepno, 
                            penalty.linear = penalty,
                            criterion = criterion,
                            family = family, 
                            ...)
  } else {
    set.seed(seed)
    gamb <- GAMBoost(y=time,
                 x.linear = x,
                 penalty.linear=penalty,
                 criterion = criterion, family = family,stepno = stepno,...)
  }
  selected <- getGAMBoostSelected(gamb)
  indmain <- selected$parametric
  res <- list()
  res$model <- gamb
  res$indmain <- indmain
  res$xnames <- colnames(x)[indmain]
  res$beta <- gamb$beta.linear[gamb$stepno + 1, selected$parametric]
  return(res)
}


fit.uniCox <- function(time, status, x, 
    unpen.index = NULL, 
    method = 'bonferroni', sig = 0.05,...){
  
  # Identifying main effects:
    pvalue <- beta <- rep(NA, ncol(x))
    for (i in 1:ncol(x)){
       res <- coxph(Surv(time,status)~ x[,i])
       pvalue[i] <- summary(res)$coefficients[,5]
       beta[i] <- summary(res)$coefficients[,1]
    }
    pvalueadjust <- p.adjust(p=pvalue, 
                             method = method)
    indmain <- unique(c(which(pvalueadjust < sig), 
                         unpen.index))
    if (is.null(indmain)) 
      stop("Error: no main effects are determined. \n Increase the significance level.")
  # Fitting Cox Model:
    datamulti <- as.data.frame(x[,indmain])
    cox <- coxph(Surv(time,status)~., 
                 data= datamulti)
  # Results:
    res <- list()
    res$model <- cox
    res$xnames <- colnames(x)[indmain]
    res$indmain <- indmain
    res$beta <- coefficients(cox)
    
    return(res)
    
}



fit.uniGlm <- function(time, status, x, 
                       unpen.index = NULL, 
                       method = 'bonferroni', family = gaussian(), sig = 0.05,...){
  
  # Identifying main effects:
  pvalue <- beta <- rep(NA, ncol(x))
  for (i in 1:ncol(x)){
    if (sum(x[,i]) != 0){
      res <- glm(time~ x[,i], family = family)
      pvalue[i] <- summary(res)$coefficients[2,4]
      beta[i] <- summary(res)$coefficients[2,1]
    } else {
      pvalue[i] <- 1
      beta[i] <- NA
    }  
  }
  pvalueadjust <- p.adjust(p=pvalue, 
                           method = method)
  indmain <- unique(c(which(pvalueadjust < sig), 
                      unpen.index))
  # Fitting Cox Model:
  datamulti <- as.data.frame(x[,indmain])
  glm <- glm(time~., 
               data= datamulti)
  # Results:
  res <- list()
  res$model <- glm
  res$xnames <- colnames(x)[indmain]
  res$indmain <- indmain
  res$beta <- coefficients(glm)
  
  return(res)
  
}


fit.hierarch <- function(nr, data, indices, seed.interselect,
                        method = 'bonferroni', family = gaussian(), sig = 0.05,...){
  
  # Identifying main effects:
  time <- data$time
  status <- data$status
  pvalue <- beta <- rep(NA, ncol(data)-2)
  
  if (sum(data$status) != 0){
    for (i in 1:(ncol(data)-2)){
      res <- coxph(Surv(time,status)~ data[,i])
      pvalue[i] <- summary(res)$coefficients[,5]
      beta[i] <- summary(res)$coefficients[,1]
    }
  } else {
    for (i in 1:(ncol(data)-2)){
      res <- glm(time~ data[,i], family = family)
      pvalue[i] <- summary(res)$coefficients[2,4]
      beta[i] <- summary(res)$coefficients[2,1]  
    }
  }
  
  pvalueadjust <- p.adjust(p=pvalue, 
                           method = method)
  indmain <- unique(which(pvalueadjust < sig))
  vim <- rep(0, ncol(data))
  vim[indmain] <- 1
  return(vim)
}




fit.rf <- function(nr, data, indices, 
                   seed.interselect,...){
  cat(nr)
  if(sum(data$status) == 0){ 
    rf  <- rfsrc(
      time ~ ., 
      data = data[indices,-(ncol(data)-1)],
      big.data = T, 
      seed = seed.interselect
    )
    } else {
    rf  <- rfsrc(
      Surv(time, status) ~ ., 
      data = data[indices,],
      big.data = T, 
      seed = seed.interselect
      ) 
  }
  return(rf$importance)
}


#!! noch nicht für stetige Outcomes modifiziert:
fit.rf.select <- function(nr, data, indices, 
                   seed.interselect,n.select, ...){
  pvalue <- rep(NA, ncol(data))
  for (i in 1:ncol(data)){
    res <- coxph(Surv(data$time,data$status)~ data[,i])
    pvalue[i] <- summary(res)$coefficients[,5]
  }
  
  sig <- sort(pvalue)[min(n.select,ncol(data)-2)]
  indmain <- which(pvalue <=  sig )
  
  data.select <- data[,c(indmain,(ncol(data)-1):ncol(data))]
  cat(nr)
  rsf  <- rfsrc(
    Surv(time, status) ~ ., 
    data = data.select[indices,],
    big.data = T, 
    seed = seed.interselect
  )
  return(rsf$importance)
}



fit.logicReg <- function(nr, 
                         data, 
                         indices, 
                         seed.interselect,
                         type=4,
                         nleaves=16,
                         ntrees=4,...){
  cat(nr)
  # Daten umformen in binäre Variablen:
  d <- sapply(data[,1:(ncol(data)-2)],
            function(arg){
                bin <- rep(NA, length(arg))
                bin[which(arg < median(arg))] <- 0
                bin[which(arg >= median(arg))] <- 1
                return(bin)
              }
    )
  d <- as.data.frame(d)
  obs.status <- data[,ncol(data)-1]
  obs.time <- data[,ncol(data)]
  
  # Durchführung logicRegression:
  set.seed(seed.interselect)
  myanneal <- logreg.anneal.control(start = -1, end = -4, iter = 2500, update = 100)
  if (sum(obs.status) == 0){
    logicR <- logreg(resp = obs.time[indices], bin=d[indices,],
                     anneal.control = myanneal,select = 2,
                     type=type,
                     nleaves=nleaves,
                     ntrees=ntrees) 
  } else {
    logicR <- logreg(resp = obs.time[indices], bin=d[indices,], cens = obs.status[indices], 
                   anneal.control = myanneal,select = 2,
                   type=type,
                   nleaves=nleaves,
                   ntrees=ntrees) 
  }
  
  # Bewerten der Variablen, die in logicRegression ausgewählt werden mit vimp = 1:
  # Code läuft nur bei select = 2:
  vim <- rep(0,ncol(d))
  for (i in 1:ntrees){
    knot <- logicR$alltrees[[1]]$trees[[i]]$trees$knot 
    vim[knot[which(knot > 0)]] <- 1
  }
  
  # Ergebnisse zusammenstellen:
  return(vim)
  
}


#!! noch nicht für stetige Outcomes modifiziert:
fit.logicReg.select <- function(nr, 
                         data, 
                         indices, 
                         seed.interselect,
                         type=4,
                         nleaves=16,
                         ntrees=4,
                         n.select,...){
  cat(nr)
  # restrict Data
  pvalue <- rep(NA, ncol(data))
  for (i in 1:ncol(data)){
    res <- coxph(Surv(data$time,data$status)~ data[,i])
    pvalue[i] <- summary(res)$coefficients[,5]
  }
  
  sig <- sort(pvalue)[min(n.select,ncol(data)-2)]
  indmain <- which(pvalue <=  sig )
  
  data.select <- data[,c(indmain,(ncol(data)-1):ncol(data))]
  
  
  # Daten umformen in binäre Variablen:
  d <- sapply(data.select[,1:(ncol(data.select)-2)],
              function(arg){
                bin <- rep(NA, length(arg))
                bin[which(arg < median(arg))] <- 0
                bin[which(arg >= median(arg))] <- 1
                return(bin)
              }
  )
  d <- as.data.frame(d)
  obs.status <- data[,ncol(data)-1]
  obs.time <- data[,ncol(data)]
  
  # Durchführung logicRegression:
  set.seed(seed.interselect)
  myanneal <- logreg.anneal.control(start = -1, end = -4, iter = 2500, update = 100)
  logicR <- logreg(resp = obs.time[indices], bin=d[indices,], cens = obs.status[indices], 
                   anneal.control = myanneal,select = 2,
                   type=type,
                   nleaves=nleaves,
                   ntrees=ntrees) 
  
  # Bewerten der Variablen, die in logicRegression ausgewählt werden mit vimp = 1:
  # Code läuft nur bei select = 2:
  vim <- rep(0,ncol(d))
  for (i in 1:ntrees){
    knot <- logicR$alltrees[[1]]$trees[[i]]$trees$knot 
    vim[knot[which(knot > 0)]] <- 1
  }
  
  # Ergebnisse zusammenstellen:
  return(vim)
  
}



simul.int <- function(seed, n = 100, p = 1000,
                              n.main = 2,
                              n.int = 2,
                              beta.main=2, 
                              beta.int = 4, 
                              censparam = 1/5, 
                              lambda = 1/20){
  # Create seeds:
  set.seed(seed+70)
  seeds <-  sample(1:100000, 4)
    
  # Simulate X:
  set.seed(seeds[1])
  data <- matrix(rnorm(n*p),nrow=n,ncol=p)
  
  # Create interactions of variables without effect:
  i <- 1:n.int *2-1 +n.main 
  j <- i+1
  
  data.int <- data[,i] * data[,j]
  data.pred <- cbind(data,data.int)
  
    
  # Create Survival times (parametrisches Cox-Modell mit Exponentialverteilung):
  beta.main.vec <- c(rep(beta.main,floor(n.main/2)),rep(-beta.main,ceiling(n.main/2)))
  beta.int.vec <- c(rep(beta.int,floor(n.int/2)),rep(-beta.int,ceiling(n.int/2)))
  beta    <- c( beta.main.vec,          # main effects
                rep(0,p-n.main),        # noise
                beta.int.vec)           # interactions
  linpred  <- exp(data.pred %*% beta)
  set.seed(seeds[3])
  runifdata <- runif(n,0,1)
  obs.time <- array()
  obs.time <- (-log(runif(n,0,1))/(lambda*linpred))
  
  # Censering:
  set.seed(seeds[1])
  cens.time <- (-log(runif(n,0,1))/(censparam))
  obs.status <- ifelse(obs.time <= cens.time,1,0)
  obs.time <- ifelse(obs.time <= cens.time,obs.time,cens.time)
  
  obs.status[obs.time > 1/censparam*2] <- 0
  obs.time[obs.time > 1/censparam*2] <- 1/censparam*2
  obs.time[which(obs.time == 9)] <- 10 ^-5
  
  
  # Summarize:
  simdata <- as.data.frame(cbind(data, obs.time, obs.status))
  colnames(simdata) <- c(paste('ID',1:p,sep=''),'obs.time', 'obs.status')
  
  # Data information:
  info <- as.data.frame(
      rbind(cbind(colnames(simdata)[1:p], beta[1:p])[which(beta[1:p]!=0) ,],
        cbind(paste(colnames(simdata)[i],':',colnames(simdata)[j], sep = ""), beta.int.vec))
  )
  colnames(info) <- c('ID','Effect size')
  
  out <- list()
  out$data <- simdata
  out$info  <- info
  return(out)
}


get.info <-function(object){
  return(object$infos)
}

get.data <- function(object){
  return(object$simdata)
}

predict.sprinter <- 
function(object, newdata=NULL, ...){
  if (is.null(object)) {
    stop("object is empty!")
  }
  if (!inherits(object, c("sprinter"))) {
    stop("This function only works for objects of class `sprinter'.")
  }
  
  if (is.null(newdata)) {
    linear.predictor <- predict(object$final.model$model)
  } else {
    
    #Create Dataset with interactions:
    if (!is.null(object$inter.candidates)){
      if (length(coef(object$final.model$model)) != length(object$final.model$indmain)){
          list.final <- object$final.model$model$xnames
        } else {
          list.final <- object$final.model$xnames
        }
      list.final.split <- strsplit(list.final, ':')
      
      newdata.inter <- matrix(NA, nrow = nrow(newdata), ncol = length(list.final))
      
      for (i in 1:length(object$main.model$indmain)){
        newdata.inter[,i] <- newdata[,object$main.model$indmain[i]]
      }
      
      if (length(object$main.model$indmain)!=length(list.final.split)){
        for( i in (length(object$main.model$indmain)+1):length(list.final.split)) {
          newdata.inter[,i] <- newdata[,which(colnames(newdata) ==list.final.split[[i]][[1]])] * 
            newdata[,which(colnames(newdata) ==list.final.split[[i]][[2]])] 
        }
      }
      colnames(newdata.inter) <- list.final
      newdata <- as.data.frame(newdata.inter)
    } 
    linear.predictor <- predict(object$final.model$model, newdata = newdata)
  }
  return(linear.predictor)
}

print.sprinter <- function(x,...){
  if (is.null(x)) {
    stop("object x is empty!")
  }
  if (!inherits(x, c("sprinter"))) {
    stop("This function only works for objects of class `sprinter'.")
  }
  
  cat('Top 20 Interaction Candidates with Inclusion Frequencies: \n ')
  if (is.matrix(x$inter.candidates)){#& nrow(x$inter.candidates)>0){
    inter.candidates <- x$inter.candidates[1:min(20, nrow(x$inter.candidates)),3]
    names(inter.candidates) <- paste(x$xnames[x$inter.candidates[1:min(20, nrow(x$inter.candidates)),1]],':',
                                     x$xnames[x$inter.candidates[1:min(20, nrow(x$inter.candidates)),2]], sep='')
  } else {
    inter.candidates <- x$inter.candidates[3]
    names(inter.candidates) <- paste(x$xnames[x$inter.candidates[1]],':',
                                     x$xnames[x$inter.candidates[2]], sep='')
  }
  print(inter.candidates)
  
  cat('\n Final model: \n \n')
  print(x$final.model$model)
  
}


summary.sprinter <- function(object, ...){
  if (is.null(object)) {
    stop("object is empty!")
  }
  if (!inherits(object, c("sprinter"))) {
    stop("This function only works for objects of class `sprinter'.")
  }
  
  main.candidates <- coef(object$main.model$model)[which(coef(object$main.model$model) != 0)]
  
  inter.candidates <- object$inter.candidates[,3]
  names(inter.candidates) <- paste(object$xnames[object$inter.candidates[,1]],':',
                                   object$xnames[object$inter.candidates[,2]], sep='')


  cat('Main candidates with coefficients: \n')
  print(main.candidates)
  cat('\n Top 20 Interaction Candidates with Inclusion Frequencies: \n')
  print(inter.candidates[1:20])
  cat('\n Final Model: \n')
  summary(object$final.model$model)
}




summary.resample.sprinter <- function(object, print = TRUE, plot = TRUE, optional.only = FALSE, threshold.vif = 0,...){
  
  if (is.null(object)) {
    stop("object is empty!")
  }
  if (!inherits(object, c("resample.sprinter"))) {
    stop("This function only works for objects of class `resample.sprinter'.")
  }
  
  
  # ausgew?hlte Variablen mit VIFs
  fold <- length(object)
  tab.var <-                                          
    sort(
      table(
        unlist(
          lapply(object, 
                 function(arg){
                     arg$final.model$xnames
                  }
                 )
        )
      )
    , decreasing = T
    )/fold
  
  
  # Interaktionen mit VIFs
  splitnames.var <- strsplit(names(tab.var), ':')   
  tab.var.inter.index <- unlist(lapply(splitnames.var, length)) == 2
  tab.inter <- tab.var[tab.var.inter.index ]       
  
  
  # Betas aller Variablen
  betas.var <- 
    unlist(
      lapply(object, 
             function(arg){
              if (length(coef(arg$final.model$model)) != length(arg$final.model$indmain)){
                main <- coef(arg$final.model$model)[arg$final.model$indmain]
              } else {
                main <- coef(arg$final.model$model)
                names(main) <- arg$final.model$xnames
              }
              return(main)
            } 
      )
    )
  
  
  # Mean betas aller Variablen
  betas.split <- split(betas.var, names(betas.var))
  d <- sapply(betas.split, function(arg){
    vif <- length(arg)/fold
    m <- mean(arg)
    return(c(vif,m))
  })
  
  # korrespondierende VIFs und Mean betas zu betas.var   
  dmerge <- as.data.frame(list(betas = betas.var, 
                             vifs = d[1,match(as.character(names(betas.var)), as.character(colnames(d)))],
                             betamean = d[2,match(as.character(names(betas.var)), as.character(colnames(d)))],
                             names = names(betas.var)))
  
    
  # Betas und vifs der Interaktionen
  splitnames.vif.var <- strsplit(as.character(dmerge$names), ':')   
  vif.var.inter.index <- unlist(lapply(splitnames.vif.var, length)) == 2
  betas.inter <- betas.var[vif.var.inter.index]
  vif.inter <- dmerge[vif.var.inter.index,2]
  names(vif.inter) <- names(betas.var[vif.var.inter.index])
  
  
  
  
  if (sum(vif.var.inter.index)>0){
    if (plot){
      if (optional.only | threshold.vif != 0){
        if (threshold.vif>0){
          index <- -which(dmerge$vifs < threshold.vif)
          }
        if (optional.only){
          index <- -which(dmerge$names %in% object[[1]]$mandatory)
        }
        if ((threshold.vif >0) & optional.only){
          index <- c(-which(dmerge$vifs < threshold.vif), index)
        }
        par(mar = c(5,4,10,2)+0.1)
        plot(x=dmerge[index,c(3,1)], axes = T, xlab = 'Mean coefficient', ylab = 'coefficient', col = rgb(0,0,0, 0.6))  
        axis(3, labels = dmerge$names[index], at = dmerge$betamean[index], las = 3)
        axis(2)
      } else {
        par(mar = c(5,4,10,2)+0.1)
        plot(x=dmerge[,c(3,1)], axes = T, xlab = 'Mean coefficient', ylab = 'coefficient', col = rgb(0,0,0, 0.6))  
        axis(3, labels = dmerge$names, at = dmerge$betamean, las = 3)
        axis(2)
      }
      }
  
    if (print){
      dmerge2 <- dmerge[vif.var.inter.index,c(4,2,3)]
      inter <- dmerge2[!duplicated(dmerge2[,1]),]
      rownames(inter) <- inter[,1]
      colnames(inter) <- c('Names','Variable Inclusion Frequency','Mean Coefficient')
      
      inter <- inter[order(inter[,2], decreasing = TRUE),c(2,3)]
      interprint <- inter[which(inter[,1] > 1/length(object)),]
      if(nrow(interprint)!=0){
        print(interprint)
      } else {
        cat('No interactions with variable inclusion frequency larger than 1/fold are selected.')
      }
      
    }
  } else {
    message("No interactions have been selected.")
  }
}


print.resample.sprinter<- function(x,...){
  summary.resample.sprinter(x, print = TRUE, plot = FALSE)
}

plot.resample.sprinter <- function(x, optional.only = FALSE, threshold.vif = 0, ...){
  summary.resample.sprinter(x, print = FALSE, plot = TRUE, optional.only, threshold.vif)
}

