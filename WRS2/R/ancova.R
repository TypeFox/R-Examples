ancova <- function(formula, data, tr = 0.2, fr1 = 1, fr2 = 1, pts = NA){

  alpha <- .05
  xout <- FALSE
  LP <- TRUE
  pr <- TRUE
  
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
    change <- TRUE
  } else {
    change <- FALSE
  }
 
  xy=elimna(cbind(x1,y1))
  x1=xy[,1]
  y1=xy[,2]
  xy=elimna(cbind(x2,y2))
  x2=xy[,1]
  y2=xy[,2]
 
  ## --- smoothing for each group
  fitted1 <- runmean(x1, y1, fr = fr1, tr = tr)
  fitted2 <- runmean(x2, y2, fr = fr2, tr = tr)
  fitted.values <- list(fitted1, fitted2)
  names(fitted.values) <- grnames
  
  if(is.na(pts[1])){
    npt<-5
    isub<-c(1:5)  # Initialize isub
    test<-c(1:5)
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
    mat<-matrix(NA,5,10)
    dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","TEST","se","ci.low","ci.hi","p.value","crit.val"))
    
    for (i in 1:5){
      g1<-y1[near(x1,x1[isub[i]],fr1)]
      g2<-y2[near(x2,x1[isub[i]],fr2)]
      g1<-g1[!is.na(g1)]
      g2<-g2[!is.na(g2)]
      test<-yuen1(g1,g2,tr=tr)
      mat[i,1]<-x1[isub[i]]
      mat[i,2]<-length(g1)
      mat[i,3]<-length(g2)
      mat[i,4]<-test$dif
      mat[i,5]<-test$teststat
      mat[i,6]<-test$se
      critv<-NA
      if(alpha==.05)critv<-smmcrit(test$df,5)
      if(alpha==.01)critv<-smmcrit01(test$df,5)
      if(is.na(critv))critv<-smmval(test$df,5,alpha=alpha)
      mat[i,10]<-critv
      cilow<-test$dif-critv*test$se
      cihi<-test$dif+critv*test$se
      mat[i,7]<-cilow
      mat[i,8]<-cihi
      mat[i,9]<-test$p.value
    }}
  
  if(!is.na(pts[1])){
    if(length(pts)>=29)stop("At most 28 points can be compared")
    n1<-1
    n2<-1
    vecn<-1
    for(i in 1:length(pts)){
      n1[i]<-length(y1[near(x1,pts[i],fr1)])
      n2[i]<-length(y2[near(x2,pts[i],fr2)])
    }
    mat<-matrix(NA,length(pts),10)
    dimnames(mat)<-list(NULL,c("X","n1","n2","DIF","TEST","se","ci.low","ci.hi",
                               "p.value","crit.val"))
    for (i in 1:length(pts)){
      g1<-y1[near(x1,pts[i],fr1)]
      g2<-y2[near(x2,pts[i],fr2)]
      g1<-g1[!is.na(g1)]
      g2<-g2[!is.na(g2)]
      test<-yuen1(g1,g2,tr=tr)
      mat[i,1]<-pts[i]
      mat[i,2]<-length(g1)
      mat[i,3]<-length(g2)
      if(length(g1)<=5)print(paste("Warning, there are",length(g1)," points corresponding to the design point X=",pts[i]))
      if(length(g2)<=5)print(paste("Warning, there are",length(g2)," points corresponding to the design point X=",pts[i]))
      mat[i,4]<-test$dif
      mat[i,5]<-test$teststat
      mat[i,6]<-test$se
      if(length(pts)>=2)critv<-smmcrit(test$df,length(pts))
      if(length(pts)==1)critv<-qt(.975,test$df)
      cilow<-test$dif-critv*test$se
      cihi<-test$dif+critv*test$se
      mat[i,7]<-cilow
      mat[i,8]<-cihi
      mat[i,9]<-test$p.value
      mat[i,10]<-critv
    }}
  
  
  mat <- as.data.frame(mat)
  if(change) grnames <- grnames[2:1]
  result <- list(evalpts = mat$X, n1 = mat$n1, n2 = mat$n2, trDiff = mat$DIF, se = mat$se, ci.low = mat$ci.low, ci.hi = mat$ci.hi, 
                 test = mat$TEST,  crit.vals = mat$crit.val, p.vals = mat$p.value, fitted.values = fitted.values,
                 cnames = colnames(mf[,c(1, datfac, datcov)]), 
                 faclevels = grnames, call = mcl)
  class(result) <- c("ancova")
  result
}
