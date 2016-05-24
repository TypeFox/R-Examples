prm_dcv   <- function(  X,                            # matrix for X
                        Y,                            # vector for y
                        a=10,                     # max number of PLS components
                        repl=10,                       # repetitions in rdcv
                        segments0=4,                  # number of outer segments
                        segments=7,                   # number of inner segments
                        segment0.type="random",       # type of outer segmentation 
                        segment.type="random",        # type of inner segmentation
                        sdfact=2,                     # parsimony factor for a_opt
                        fairct=4,                     # factor for robust pls weights
                        trim=0.2,                     # trim factor for robust pls
                        opt="median",                 # calculation method for multivariate median
                        plot.opt=FALSE, ...)
{
#require(chemometrics)
n = nrow(X)                                             # number of samples

ncomp <- a
optcomp <- matrix(NA, nrow=segments0, ncol=repl)      # a_opt cv [1:segments0, 1:repl] 
b       <- matrix(NA, nrow=dim(X)[2], ncol=ncomp)
bAll       <- array(NA, dim=c(segments0,dim(X)[2],ncomp,repl)) # all regression coefficients
b0All       <- array(NA, dim=c(segments0,ncomp,repl))  # all intercept terms
pred    <- array(NA, dim=c(n, ncomp, repl))            # y_test [1:n, 1:ncomp, 1:repl]
predopt <- matrix(NA, nrow=n, ncol=repl)              # y_test for each a_opt [1:n, 1:repl] 


  # (1) --- start repetition loop --- #
  for(i in 1:repl)                                         
  {
    cat("\n", "repl: ", i, "\n")
    segment0 <- cvsegments(n, k=segments0, type=segment0.type)      # create outer segments
  
            
            # (2) --- start outer loop --- #
            for(n.seg0 in 1:length(segment0))                       
            {
            cat("\n", "seg-nr: ", n.seg0, "\n")                     
            seg0 <- segment0[[n.seg0]]                              # test set samples
            obsuse <- as.numeric(unlist(segment0[-n.seg0]))         # calibration set samples
            d1 <- list(X=X[obsuse,],Y=Y[obsuse])                    # data for calibration
  
                # (3) --- start inner loop --- #
                  res <-  prm_cv(d1$X, d1$Y,                   # robust PLS with (inner segments)-fold CV with calibration data
                          a=ncomp,fairct=fairct, 
                          opt=opt, segments=segments,
                          segment.type=segment.type, trim=trim,
                          sdfact=sdfact, plot.opt=plot.opt)
                  
                  optcomp[n.seg0,i] <- res$optcomp                  # a_opt cv for current calibration set
                # (3) --- end inner loop --- #
          
             
             # Extract robust PLS models' regression coefficients 
             b0 <- vector(length=ncomp)                              # robust intercept
             for(n.comp in 1:ncomp)                                  # for all numbers of PLS components 
             {
                rcal <- prm(d1$X, d1$Y, a=n.comp, fairct=fairct, opt=opt)   # robust PLS with entire current calibration set
    
                b[,n.comp]  <- rcal$coef                             # robust regr.coeff
                b0[n.comp] <- rcal$intercept                         # robust intercept
             }
             b0All[n.seg0,,i] <- b0
             #b0 <- matrix(rep(b0, nrow(X[-obsuse,])), ncol=ncomp, byrow=TRUE) 
             bAll[n.seg0,,,i] <- b 
             # test set predicted y  
             # pred[-obsuse,,i]  <-  as.matrix(scale(X[-obsuse,], center=rcal$mx, scale=FALSE ))%*% b + b0 

             # Intercept correction for test set predicted y values:
             for(n.comp in 1:ncomp){
               pred[-obsuse,n.comp,i]  <-  as.matrix(X[-obsuse,])%*% b[,n.comp] 
               b0corr <- median(Y[-obsuse] - pred[-obsuse,n.comp,i])
               pred[-obsuse,n.comp,i] <- pred[-obsuse,n.comp,i] + b0corr
             }


             predopt[-obsuse,i] <- drop(pred[-obsuse,optcomp[n.seg0,i],i])
             }
             # (2) --- end outer loop --- #
  }
  # (1) --- end repetition loop --- #


resopt <- predopt - c(Y)
residcomp <- pred - c(Y)


# a_final as most frequent a_opt cv
afinaldistr <- table(optcomp)/sum(table(optcomp))
afinal <- as.numeric(names(which.max(afinaldistr)))



###### SEP values:
SEPall <- matrix(NA,nrow=ncomp,ncol=repl) # all untrimmed SEP values
SEPtrim <- matrix(NA,nrow=ncomp,ncol=repl) # all trimmed SEP values
SEPcomp <- rep(NA,ncomp) # trimmed SEP for each number of components
for (n.comp in 1:ncomp){
  SEPall[n.comp,] <- apply(residcomp[,n.comp,],2,sd)
  SEPtrim[n.comp,] <- sd_trim(residcomp[,n.comp,],trim=trim)
  SEPcomp[n.comp] <- sd_trim(as.vector(residcomp[,n.comp,]),trim=trim)
}

bAll=apply(bAll,c(2,3),mean)
b0All=apply(b0All,c(2),mean)
list(b=bAll,intercept=b0All,resopt=resopt, predopt=predopt, optcomp=optcomp, residcomp=residcomp, 
     pred=pred, SEPall=SEPall, SEPtrim=SEPtrim, SEPcomp=SEPcomp,  afinal=afinal, SEPopt=SEPcomp[afinal])
}
 
 
