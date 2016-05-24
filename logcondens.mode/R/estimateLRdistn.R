##Charles Doss
## Likelihood ratio of location-of-mode for log-concave class.

##mode is location of constrained mode; LR goes to  Op(1) or infty if true or false,resp.
## N.MC is number of Monte Carlos of the LR.
## n.SS is the Sample Size for each Monte Carlo of the LR.
##Returned:
## LR is vector of length N.MC of likelihood ratios.
## TLLR stands for "two-log-likelihood ratio" = 2*log(LR)

## debugging turns on and off debugging.  It is NULL for no debugging or a
## character string, rather than a logical, though.  Usually it does not end
## in an extension; rather, e.g. "abc" will create files abc.txt and abc.rsav 
## 

estimateLRdistn <-
  function(rdist=rnorm,mode=0, N.MC=1e2, n.SS=1e4,
           xgrid=NULL,
           prec=10^-10, ##Note: setting prec to smaller values has created errors for me!
           seedVal=NULL,
           debugging=NULL){
    aval <- mode; ##last minute renaming of variables from 'aval' to 'mode'
    if (!is.null(seedVal)) set.seed(seedVal)
    if (is.character(debugging)) { ##NULL fails
      ##sink("getLRlog.txt")
      sink(paste(debugging, "_log.txt", sep=""))
      currTime <- prevTime <- 0
      print(paste("Parameters are: N.MC=", N.MC,
                  "// n.SS=", n.SS, "// aval=", aval,
                  " and you can load rdist from the rsav file, ",
                  debugging, ".rsav.",
                  sep=""))
    }
    LRs <- TLLRs <- vector(length=N.MC)
    ##myww <- rep(1/n.SS, n.SS)
    for (i in 1:N.MC){

      myxx <- rdist(n.SS)
      tmp <- preProcess(x=myxx, xgrid=xgrid)
      myxx <- tmp$x ##preprocess makes unique
      myww <- tmp$w
      ## n.SS <- tmp$n ##same
      
      
      ## myxx <- sort(rdist(n.SS))
      ## ##while (!identical(myxx,unique(myxx))) myxx <- sort(rdist(n.SS));
      ## myxx.uniq <- rle(myxx) ##myxx is sorted
      if (is.character(debugging)) {
        prevTime <- currTime
        currTime <- Sys.time()
        save(file=paste(debugging, ".rsav", sep=""),
             myxx,
             ##myxx.uniq,
             xgrid,
             rdist) ##save most recent, for hang/crash
        print(paste("Starting the ", i, "th iteration.",
                    "Memory usage is ", sum(gc()[1:2,2]), ##conscells+heap in MB
                    " Printing the time for the last iteration ",
                    "and the current time.",
                    sep=""))
        print(currTime-prevTime)
        print(currTime);
      }
      ## myxx <- myxx.uniq$values
      ## myww <- myxx.uniq$lengths / n.SS
      res.UC <- activeSetLogCon(x=myxx,w=myww, prec=prec,print=F) ##UnConstrained
      res.MC <- activeSetLogCon.mode(x=myxx,mode=aval,w=myww,
                                     prec =prec,
                                     print=F) ##Mode Constrained
      ##TLLRs[i] <- 2*(sum(res.UC$phi) - sum(res.MC$phi))
      TLLRs[i] <- 2 * n.SS * (res.UC$L - res.MC$L);
      ## if (res.MC$dlcMode$isx) {
      ##   print("getLR Warning: res.MC$dlcMode$isx is true. This is not standard.")
      ##   TLLRs[i] <- 2*(sum(res.UC$phi - res.MC$phi))
      ## }
      ## else {TLLRs[i] <- 2*(sum(res.UC$phi - res.MC$phi[-res.MC$dlcMode$idx]))}
      if (is.character(debugging)){
        if (TLLRs[i] <= 0 || TLLRs[i] > 10e4){
          print("getLR ERROR: 2log-LR returned value <=0 or > 10e4");
          tmpfile <- paste(debugging,"tmpxxs", i, ".rsav", sep="");
          save(myxx,myww,aval, xgrid, file=tmpfile)
          print(paste("Saved myxxs to ", getwd(),  tmpfile, sep=""))
        }
      }
      ##rm(res.UC)
      ##rm(res.MC)
      ##res.UC.phi <- activeSetLogCon(x=myxx,w=myww,print=F)$phi ##"UnConstrained
      ##res.MC.phi <- activeSetLogCon.mode(x=myxx,aval=aval,w=myww, print=F)$phi ##Mode Constrained
      ##TLLRs[i] <- 2*(sum(res.UC.phi) - sum(res.MC.phi))
      ##rm(res.UC.phi,res.MC.phi,myxx)
    }
    if (is.character(debugging)) sink()
    ##return(list(LRs=exp(TLLRs/2), TLLRs=TLLRs, lastRes.MC =res.MC, lastRes.UC =res.UC, lastSeed =tmp, xs=myxx))
    return(list(LRs=exp(TLLRs/2), TLLRs=TLLRs))
  }

## estimateLRdistn <-
##   function(rdist=rnorm,aval=0, N.MC=1e2, n.SS=1e4,
##            prec=10^-10, ##Note: setting prec to smaller values has created errors for me!
##            seedVal=NULL,
##            debugging=NULL){
##     if (!is.null(seedVal)) set.seed(seedVal)
##     if (is.character(debugging)) { ##NULL fails
##       ##sink("getLRlog.txt")
##       sink(paste(debugging, "_log.txt", sep=""))
##       currTime <- prevTime <- 0
##       print(paste("Parameters are: N.MC=", N.MC,
##                   "// n.SS=", n.SS, "// aval=", aval,
##                   " and you can load rdist from the rsav file, ",
##                   debugging, ".rsav.",
##                   sep=""))
##     }
##     LRs <- TLLRs <- vector(length=N.MC)
##     ##myww <- rep(1/n.SS, n.SS)
##     for (i in 1:N.MC){
##       myxx <- sort(rdist(n.SS))
##       ##while (!identical(myxx,unique(myxx))) myxx <- sort(rdist(n.SS));
##       myxx.uniq <- rle(myxx) ##myxx is sorted
##       if (is.character(debugging)) {
##         prevTime <- currTime
##         currTime <- Sys.time()
##         save(file=paste(debugging, ".rsav", sep=""),
##              myxx,myxx.uniq,rdist) ##save most recent, for hang/crash
##         print(paste("Starting the ", i, "th iteration.",
##                     "Memory usage is ", sum(gc()[1:2,2]), ##conscells+heap in MB
##                     " Printing the time for the last iteration ",
##                     "and the current time.",
##                     sep=""))
##         print(currTime-prevTime)
##         print(currTime);
##       }
##       myxx <- myxx.uniq$values
##       myww <- myxx.uniq$lengths / n.SS
##       res.UC <- activeSetLogCon(x=myxx,w=myww, prec=prec,print=F) ##UnConstrained
##       res.MC <- activeSetLogCon.mode(x=myxx,aval=aval,w=myww,
##                                      prec =prec,
##                                      print=F) ##Mode Constrained
##       ##TLLRs[i] <- 2*(sum(res.UC$phi) - sum(res.MC$phi))
##       TLLRs[i] <- 2 * n.SS * (res.UC$L - res.MC$L);
##       ## if (res.MC$dlcMode$isx) {
##       ##   print("getLR Warning: res.MC$dlcMode$isx is true. This is not standard.")
##       ##   TLLRs[i] <- 2*(sum(res.UC$phi - res.MC$phi))
##       ## }
##       ## else {TLLRs[i] <- 2*(sum(res.UC$phi - res.MC$phi[-res.MC$dlcMode$idx]))}
##       if (is.character(debugging)){
##         if (TLLRs[i] <= 0 || TLLRs[i] > 10e4){
##           print("getLR ERROR: 2log-LR returned value <=0 or > 10e4");
##           tmpfile <- paste(debugging,"tmpxxs", i, ".rsav", sep="");
##           save(myxx,myww,aval, file=tmpfile)
##           print(paste("Saved myxxs to ", getwd(),  tmpfile, sep=""))
##         }
##       }
##       ##rm(res.UC)
##       ##rm(res.MC)
##       ##res.UC.phi <- activeSetLogCon(x=myxx,w=myww,print=F)$phi ##"UnConstrained
##       ##res.MC.phi <- activeSetLogCon.mode(x=myxx,aval=aval,w=myww, print=F)$phi ##Mode Constrained
##       ##TLLRs[i] <- 2*(sum(res.UC.phi) - sum(res.MC.phi))
##       ##rm(res.UC.phi,res.MC.phi,myxx)
##     }
##     if (is.character(debugging)) sink()
##     ##return(list(LRs=exp(TLLRs/2), TLLRs=TLLRs, lastRes.MC =res.MC, lastRes.UC =res.UC, lastSeed =tmp, xs=myxx))
##     return(list(LRs=exp(TLLRs/2), TLLRs=TLLRs))
##   }
