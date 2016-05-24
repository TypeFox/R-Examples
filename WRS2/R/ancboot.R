ancboot <- function(formula, data, tr = 0.2, nboot = 599, fr1 = 1, fr2 = 1, pts = NA){
 
  # Confidence intervals are computed using a percentile t bootstrap
  # method. Comparisons are made at five empirically chosen design points.
  #
  #  Assume data are in x1 y1 x2 and y2
  #
  alpha <- .05
  xout <- FALSE
  LP <- TRUE
  pr <- TRUE
  sm <- FALSE
  
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  mcl <- match.call()
  
  if(is.factor(mf[,2])) {
    datfac <- 2
    datcov <- 3
  } else {
    datfac <- 3
    datcov <- 2
  }
  
  grnames <- levels(mf[,datfac])
  if (is.null(grnames)) stop("Group variable needs to be provided as factor!")
  if (length(grnames) > 2) stop("Robust ANCOVA implemented for 2 groups only!")
  yy <- split(mf[,1], mf[, datfac])
  y1 <- yy[[1]]
  y2 <- yy[[2]]
  xx <- split(mf[,datcov], mf[, datfac])
  x1 <- xx[[1]]
  x2 <- xx[[2]]
  
  ## bugfix
  if(length(x1) < length(x2)) {              
    dummy <- x1
    x1 <- x2
    x2 <- dummy
    dummy <- y1
    y1 <- y2
    y2 <- dummy
    change = TRUE
  } else {
    change <- FALSE
  }
  
  
  if(is.na(pts[1])){
    isub<-c(1:5)  # Initialize isub
    test<-c(1:5)
    m1=elimna(cbind(x1,y1))
    x1=m1[,1]
    y1=m1[,2]
    m1=elimna(cbind(x2,y2))
    x2=m1[,1]
    y2=m1[,2]
    xorder<-order(x1)
    y1<-y1[xorder]
    x1<-x1[xorder]
    xorder<-order(x2)
    y2<-y2[xorder]
    x2<-x2[xorder]
    n1<-1
    n2<-1
    vecn<-1
    for(i in 1:length(x1))n1[i]<-length(y1[near(x1,x1[i],fr1)])
    for(i in 1:length(x1))n2[i]<-length(y2[near(x2,x1[i],fr2)])
    for(i in 1:length(x1))vecn[i]<-min(n1[i],n2[i])
    sub<-c(1:length(x1))
    isub[1]<-min(sub[vecn>=12])
    isub[5]<-max(sub[vecn>=12])
    isub[3]<-floor((isub[1]+isub[5])/2)
    isub[2]<-floor((isub[1]+isub[3])/2)
    isub[4]<-floor((isub[3]+isub[5])/2)
    mat<-matrix(NA,5,8)
    dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","TEST","ci.low","ci.hi",
                               "p.value"))
    gv1<-vector("list")
    for (i in 1:5){
      j<-i+5
      temp1<-y1[near(x1,x1[isub[i]],fr1)]
      temp2<-y2[near(x2,x1[isub[i]],fr2)]
      temp1<-temp1[!is.na(temp1)]
      temp2<-temp2[!is.na(temp2)]
      mat[i,2]<-length(temp1)
      mat[i,3]<-length(temp2)
      gv1[[i]]<-temp1
      gv1[[j]]<-temp2
    }
    I1<-diag(5)
    I2<-0-I1
    con<-rbind(I1,I2)
    test<-linconb(gv1,con=con,tr=tr,nboot=nboot)
    for(i in 1:5){
      mat[i,1]<-x1[isub[i]]
    }
    mat[,4]<-test$psihat[,2]
    mat[,5]<-test$test[,2]
    mat[,6]<-test$psihat[,3]
    mat[,7]<-test$psihat[,4]
    mat[,8]<-test$test[,4]
  }
  if(!is.na(pts[1])){
    n1<-1
    n2<-1
    vecn<-1
    for(i in 1:length(pts)){
      n1[i]<-length(y1[near(x1,pts[i],fr1)])
      n2[i]<-length(y2[near(x2,pts[i],fr2)])
      if(n1[i]<=5)paste("Warning, there are",n1[i]," points corresponding to the design point X=",pts[i])
      if(n2[i]<=5)paste("Warning, there are",n2[i]," points corresponding to the design point X=",pts[i])
    }
    mat<-matrix(NA,length(pts),9)
    dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","TEST","se","ci.low","ci.hi",
                               "p.value"))
    gv<-vector("list",2*length(pts))
    for (i in 1:length(pts)){
      g1<-y1[near(x1,pts[i],fr1)]
      g2<-y2[near(x2,pts[i],fr2)]
      g1<-g1[!is.na(g1)]
      g2<-g2[!is.na(g2)]
      j<-i+length(pts)
      gv[[i]]<-g1
      gv[[j]]<-g2
    }
    I1<-diag(length(pts))
    I2<-0-I1
    con<-rbind(I1,I2)
    test<-linconb(gv,con=con,tr=tr,nboot=nboot)
    mat[,1]<-pts
    mat[,2]<-n1
    mat[,3]<-n2
    mat[,4]<-test$psihat[,2]
    mat[,5]<-test$test[,2]
    mat[,6]<-test$test[,3]
    mat[,7]<-test$psihat[,3]
    mat[,8]<-test$psihat[,4]
    mat[,9]<-test$test[,4]
  }
  
  mat <- as.data.frame(mat)
  if(change) grnames <- grnames[2:1]
  result <- list(evalpts = mat$X, n1 = mat$n1, n2 = mat$n2, trDiff = mat$DIF, ci.low = mat$ci.low, ci.hi = mat$ci.hi, test = mat$TEST,  crit.vals = mat$crit.val, p.vals = mat$p.value, cnames = colnames(mf[,c(1, datfac, datcov)]), faclevels = grnames, call = mcl)
  class(result) <- c("ancova")
  result
}
