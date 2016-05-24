EstimateLPCFDR <- function(x,y, type, nreps=100,soft.thresh=NULL,censoring.status=NULL){
  call <- match.call()
  dat <- list(x=x,y=y,censoring.status=censoring.status)
  CheckEstimateLPCFDRFormat(dat,type,nreps,soft.thresh)
  t.est <- EstimateTFDR(x,y, type,censoring.status=censoring.status)
  if(is.null(soft.thresh)){
    lpc.obj <- LPC(x,y,type=type,censoring.status=censoring.status)
    soft.thresh <- lpc.obj$soft.thresh
  }  
  pi0 <- t.est$pi0
  # Equation in LPC paper: from left to right and top to bottom, call the probabilities p1, p2, p3, and p4
  # so, my fdr estimate will be (1-pi0)(p1-p2)/(p3-p4)
  p1s <- NULL
  p2s <- NULL
  ttstars <- NULL
  if(type=="regression"){
    for(i in 1:100) ttstars <- c(ttstars, quantitative.func(dat$x, sample(dat$y), .05)$tt)
  } else if (type=="two class") {
    for(i in 1:100) ttstars <- c(ttstars, ttest.func(dat$x, sample(dat$y), .05)$tt)
  } else if (type=="multiclass"){
    for(i in 1:100) ttstars <- c(ttstars, multiclass.func(dat$x, sample(dat$y), .05)$tt)
  } else if (type=="survival") {
    for(i in 1:100){
      oo <- sample(ncol(dat$x))
      ttstars <- c(ttstars, cox.func(dat$x, dat$y[oo], dat$censoring.status[oo], .05)$tt)
    }  
  }
  p4 <- mean(abs(ttstars))
  if(type=="regression") tt.tots <- quantitative.func(dat$x, dat$y, .05)$tt
  if(type=="two class") tt.tots <- ttest.func(dat$x, dat$y, .05)$tt
  if(type=="survival") tt.tots <- cox.func(dat$x, dat$y, dat$censoring.status, .05)$tt
  if(type=="multiclass") tt.tots <- multiclass.func(dat$x, dat$y, .05)$tt
  p3 <- mean(abs(tt.tots))
  for(i in 1:nreps){
    cat(i,fill=F)
    split <- CreateTrainTestSet(dat, type=type)
    train <- split$train
    test <- split$test
    if(type=="regression"){
      tt.test <- quantitative.func(test$x, test$y, .05)$tt
      tt.train <- quantitative.func(train$x, train$y, .05)$tt
    } else if (type=="two class"){
      tt.test <- ttest.func(test$x, test$y, .05)$tt
      tt.train <- ttest.func(train$x, train$y, .05)$tt
    } else if (type=="survival"){
      tt.test <- cox.func(test$x, test$y, test$censoring.status, .05)$tt
      tt.train <- cox.func(train$x, train$y, train$censoring.status, .05)$tt
    }  else if (type=="multiclass"){
      tt.test <- multiclass.func(test$x, test$y, .05)$tt
      tt.train <- multiclass.func(train$x, train$y, .05)$tt
    }
    lpc.train <- LPC(train$x, train$y, type=type, soft.thresh=soft.thresh,censoring.status=train$censoring.status)$lpcscores
    p1s.tmp <- NULL
    p2s.tmp <- NULL
    for(q in 1:nrow(dat$x)){
      qt <- abs(tt.train[q])
      qlpc <- abs(lpc.train[q]) 
      bigt <- abs(tt.train)>=qt
      biglpc <- abs(lpc.train)>=qlpc
      if(length(bigt)==0 || length(biglpc)==0) print(q)
      p1s.tmp <- c(p1s.tmp, mean(abs(tt.test[biglpc])))
      p2s.tmp <- c(p2s.tmp, mean(abs(tt.test[bigt])))
      if(is.na(p1s.tmp[q])) p1s.tmp[q] <- max(abs(tt.test))
      if(is.na(p2s.tmp[q])) p2s.tmp[q] <- max(abs(tt.test))
    }
    p1s <- cbind(p1s, p1s.tmp)
    p2s <- cbind(p2s, p2s.tmp)
  }
  cat("",fill=T)
  p1s <- apply(p1s, 1, mean)
  p2s <- apply(p2s, 1, mean)
  tminuslpcfdr <- (1-pi0)*(p1s-p2s)/(p3-p4)
  lpcfdrobj <- list(fdrlpc=pmax(0,pmin(1, t.est$fdrt - tminuslpcfdr)), fdrdiff=tminuslpcfdr, fdrt=t.est$fdrt, pi0=pi0,call=call,soft.thresh=soft.thresh)
  class(lpcfdrobj) <- "lpcfdrobj"
  return(lpcfdrobj)
}

