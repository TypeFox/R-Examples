intervalICC.default <-
function(r1, r2, predefined.classes=FALSE, classes, c.limits, optim.method=1){
  
  # Check if all arguments are correctly specified  
  if(missing(r1) | missing(r2))
    stop("Both r1 and r2 need to be specified")
  
  if(optim.method!=1 & optim.method!=2)
    stop("Misspecified optimization method: must be 1 or 2")
  
  if(predefined.classes){
    
    if(missing(classes))
      stop("Unspecified classes")
    if(missing(c.limits))
      stop("Unspecified classes limits")
    if(!is.matrix(c.limits) & !is.data.frame(c.limits))
      stop("Classes limits should be a matrix or a data frame")
    if(!is.vector(classes))
      stop("Classes should be a vector")
    if(!is.vector(r1) | !is.vector(r2))
      stop("r1 and r2 should be vectors")
    if(length(r1)!=length(r2))
        stop("r1 and r2 should have equal lenghts")
    
    ratings <- as.data.frame(cbind(r1,r2))     
    c.limits <- as.data.frame(c.limits)
      
    if(ncol(c.limits)!=2)
      stop("Classes limits should be a matrix or a data frame with 2 columns")
    if(length(classes)!=nrow(c.limits))
      stop("Number of classes differs from number of intervals given in c.limits")
      
    names(ratings) <- c("t1","t2")
    names(c.limits) <- c("lower","upper")
    
    if(!is.numeric(c.limits$lower) | !is.numeric(c.limits$upper))
      stop("Classes limits must be numeric")
    if(any(is.na(classes)))
      stop("Missing values in classes")
    if(any(is.na(c.limits)))
      stop("Missing values in classes limits")
    if(any(is.na(ratings))){
      warning("Missing values detected: data rows omitted from calculation")
    ratings <- na.omit(ratings)
      }
    
    ratings$t1 <- factor(ratings$t1, levels=classes)
    ratings$t2 <- factor(ratings$t2, levels=classes)
    
    if(any(is.na(ratings)))
      stop("Unrecognized class in r1 or r2")
    if(any(c.limits$lower >= c.limits$upper))
      stop("Misspecified classes limits: lower bound equal to upper bound or greater detected")
    
    c.means <- rowMeans(c.limits)
    n.classes <- nrow(c.limits)
    n.resp <- nrow(ratings)

    t.means <- mat.or.vec(n.classes,n.classes)
    t.sd <- mat.or.vec(n.classes,n.classes)
    for(i in 1:n.classes){
      for(j in i:n.classes){
        t.means[i,j] <- 0.5*(c.means[i] + c.means[j])
        t.means[j,i] <- t.means[i,j]
        t.sd[i,j] <- sqrt((c.means[i] - t.means[i,j])^2 + (c.means[j] - t.means[i,j])^2)
        t.sd[j,i] <- t.sd[i,j]
      }
    }

    t.r <- table(ratings$t1,ratings$t2)
    theta0 <- c(0,0, sum(t.means*t.r)/n.resp)
    theta0[1:2] <- c(max(sqrt(sum((t.means - theta0[3])^2*t.r)/(n.resp-1)),.Machine$double.eps+1e-10), max(sum(t.sd*t.r)/n.resp,.Machine$double.eps+1e-10))
    
    if(optim.method==1){
      est <- .intervalICC.est1(ratings,classes,c.limits,theta0)
    }else{
      est <- .intervalICC.est2(ratings,classes,c.limits,theta0)
    }
  
  
        
  }else{
      if( (!is.matrix(r1) & !is.data.frame(r1)) | (!is.matrix(r2) & !is.data.frame(r2)) )
        stop("r1 and r2 should be matrices or data frames")
      if(ncol(r1)!=2 | ncol(r2)!=2)
        stop("r1 and r2 should be matrices or data frames with 2 columns")
      if(nrow(r1)!=nrow(r2))
        stop("r1 and r2 should have equal number of rows")
      if(!is.numeric(r1[,1]) | !is.numeric(r1[,2]) | !is.numeric(r2[,1]) | !is.numeric(r2[,2]))
        stop("r1 and r2 should be numeric")
      if(any(r1[,1] >= r1[,2]) | any(r2[,1] >= r2[,2]))
        stop("Misspecified limits: lower bound equal to upper bound or greater detected")
      
      ratings <- cbind(r1,r2)
      if(any(is.na(ratings))){
        warning("Missing values detected: data rows omitted from calculation")
        ratings <- na.omit(ratings)
      }
      
      r1 <- ratings[,1:2]  
      r2 <- ratings[,3:4]  
        
      ratings.num <- cbind(rowMeans(r1), rowMeans(r2))
      theta0=c(max(sd(rowMeans(ratings.num)),.Machine$double.eps+1e-10), max(mean(apply(ratings.num,1,sd)),.Machine$double.eps+1e-10), mean(ratings.num))      

      if(optim.method==1){
        est <- .intervalICC.est3(r1,r2,theta0)
      }else{
        est <- .intervalICC.est4(r1,r2,theta0)
      }    
  }  
  
  class(est) <- "ICCfit"
  est
}

