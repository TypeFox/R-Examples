`print.ppar` <-
function(x,...)
# print method for person.parameter
# x...object of class ppar
{
  cat("\n")
  cat("Person Parameters:")
  cat("\n")
  
  if (length(x$pers.ex) > 0) {    
      X <- x$X[-x$pers.ex,]                                        #list with raw scores
      sumlist <- by(x$X[-x$pers.ex,],x$gmemb,rowSums,na.rm=TRUE)
    } else {
      X <- x$X
      sumlist <- by(x$X,x$gmemb,rowSums,na.rm=TRUE)
    }
  
  if (is.null(x$pred.list)) {                                       #no spline Interpolation
    coef.list <-  mapply(function(sm,th,se) {
                           th.u <- tapply(th,sm, function(tm) {tm[1]})     #due to rounding errors, pck out first one 
                           se.u <- tapply(se,sm, function(ss) {ss[1]})
                           sm.u <- unique(sort(sm))
                           
                           smth <- cbind(sm.u,th.u,se.u)
                           return(smth)
                         },sumlist,x$thetapar,x$se,SIMPLIFY=FALSE)
  } else {                                                          #if spline Interpolation
    #TFvec <- sapply(x$pred.list,is.null)                            #for these NA groups no spline interpolation was computed
    #predind <- (1:length(x$pred.list))[!TFvec]
    #x$pred.list <- x$pred.list[predind]
    
    coef.list <- mapply(function(sm,pl,se) {
                            se.u <- tapply(se,sm, function(ss) {ss[1]})
                            sm.u <- unique(sort(sm))
                            
                            TFvec <- pl$x %in% sm.u
                            se.ind <- 1:length(TFvec)
                            se.all <- rep(NA,length(se.ind))
                            se.all[se.ind[TFvec]] <- se.u              
                            
                            cbind(pl$x,pl$y,se.all)
                            },sumlist,x$pred.list,x$se,SIMPLIFY=FALSE)
  }
  
  if (dim(coef.list[[1]])[2] == 2) {                            #if no standard errors were computed
    coef.list <- lapply(coef.list,function(cl) {cbind(cl,NA)})
  }
  
 # if (any(is.na(x$X))) {                                       #recompute gmemb without persons excluded
 #   dichX <- ifelse(is.na(x$X),1,0)
 #   strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
 #   gmemb <- as.vector(data.matrix(data.frame(strdata)))
 # } else {
 #   gmemb <- rep(1,dim(x$X)[1])
 # }
   
  for (i in 1:length(x$thetapar)) {
    cat("\n")
    if (length(x$thetapar) > 1) {
      cat("Person NA Group:",i,"\n")
      xvec <- rep(NA, (dim(x$X)[2]))
      notNApos <- which(!is.na(as.vector(rbind(X[x$gmemb == i,])[1,])))
      xvec[notNApos] <- "x"
      cat("NA pattern:",xvec,"\n")
    }
    colnames(coef.list[[i]]) <- c("Raw Score","Estimate","Std.Error")
    rownames(coef.list[[i]]) <- rep("",dim(coef.list[[i]])[1])
    print(coef.list[[i]])  
  }
  invisible(coef.list)
}

