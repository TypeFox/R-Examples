Youden3Grp.Variance.Bootstrap <-
function(x,y,z,method="Normal",seed.seq=NULL,randomStart.N=5,NBOOT=10,alpha=0.05,t.minus.start=NULL,t.plus.start=NULL,...)
  {
    ###############################################Do NBOOT bootstrap estimation on Youden and cutoff and obtain basic quantile (1-alpha)% CI ##############################
    ########Changed from name "Bootstrap.Youden.Cutoff.Iter" to "Youden3Grp.Variance.Bootstrap"
    
    ###NBOOT:#of bootstrap samples to be sample to obtain bootstrap variance
    ##seed.seq: a sequence of random number generation seeds, length of NBOOT , for reproducibility
    ##alpha: significance level for CI
    ###seed.seq should be unique
    ###t.minus.start, t.plus.start:for use of method "EMP/KS/KS-SJ" in the YoudenIndex.OptimalCutoff.Normal.EMP.KS() function as random starting values for optimization
    ####...:arguments used in function Bootstrap.Youden.Cutoff.1run

    ###check if seed.seq is long enough
    len0 <- length(unique(seed.seq))

    if(missing(seed.seq) | is.null(seed.seq)) seed.seq <- 1:NBOOT 
    if(length(unique(seed.seq))<NBOOT)
      {
        cat("given seed.seq has length<NBOOT, use 1:NBOOT as seed.seq\n")
        seed.seq <- 1:NBOOT
      }
    

    ##remove NAs in data
    x <- na.exclude(x)
    y <- na.exclude(y)
    z <- na.exclude(z)
    
    temp.res <- sapply(1:NBOOT,function(kk) {out00 <- Bootstrap.Youden.Cutoff.1run(x=x,y=y,z=z,method=method,seed0=seed.seq[kk],randomStart.N=randomStart.N,t.minus.start=t.minus.start,t.plus.start=t.plus.start,...);as.numeric(out00)})#matrix(out00,nrow=1,dimnames=list(NULL,colnames(out00)))})
    #cat("bootstrap iter=",kk,"\n");

     if(mode(temp.res)=="list")##in some bootstrapping samples, the means of resampled xx, yy and zz may not in the correct order, need to remove those runs
      {
        len0 <- unlist(lapply(temp.res,length))
        rm.idx <- which(len0!=17)
        temp.res <- temp.res[-rm.idx]
        temp.res <- as.data.frame(do.call("rbind",temp.res))
      }else temp.res <- as.data.frame(t(temp.res))##all results from bootstrap samples,row are bootstrap iterations and columns are output from Bootstrap.Youden.Cutoff.1run, including t.minus   t.plus        Se        Sp        Sm    youden  youden.z  lambda.est t.minus.TN t.plus.TN 
    names(temp.res) <- c("mu.minus","mu0","mu.plus","s.minus","s0","s.plus","t.minus","t.plus","Se","Sp" ,"Sm","youden","youden.z","lambda.est","t.minus.TN","t.plus.TN","sampleSize")
       
    all.mean <- apply(temp.res,2,function(vec0) mean(as.numeric(vec0),na.rm=TRUE))
    all.var <- apply(temp.res,2,function(vec0) var(as.numeric(vec0),na.rm=TRUE,use="pairwise.complete.obs"))
       
    youden <- all.mean["youden"]##bootstrap average of Youden
    
    var.youden <- all.var["youden"]##bootstrap variance of Youden
    var.t.minus <- all.var["t.minus"]##bootstrap variance of t.minus
    var.t.plus <- all.var["t.plus"]##bootstrap variance of t.plus

    ###bootstrap percentile CI
    t.minus.CI <- quantile(as.numeric(temp.res$t.minus),prob=c(alpha/2,1-alpha/2),na.rm=TRUE)
    names(t.minus.CI) <- c(paste(alpha/2*100,"%",sep=""),paste(100-alpha/2*100,"%",sep=""))
    
    t.plus.CI <- quantile(as.numeric(temp.res$t.plus),prob=c(alpha/2,1-alpha/2),na.rm=TRUE)
    names(t.plus.CI) <- c(paste(alpha/2*100,"%",sep=""),paste(100-alpha/2*100,"%",sep=""))
    
    youden.CI <- quantile(as.numeric(temp.res$youden),prob=c(alpha/2,1-alpha/2),na.rm=TRUE)
    names(youden.CI) <- c(paste(alpha/2*100,"%",sep=""),paste(100-alpha/2*100,"%",sep=""))
    
    #############Fisher's z transformation on youden
    youden.z.vec <- FisherZ(temp.res$youden)
    youden.z.CI <- quantile(as.numeric(youden.z.vec),prob=c(alpha/2,1-alpha/2),na.rm=TRUE)
    names(youden.CI) <- c(paste(alpha/2*100,"%",sep=""),paste(100-alpha/2*100,"%",sep=""))
    
    youden.z <- mean(youden.z.vec,na.rm=TRUE)##should youden.z be mean(FisherZ(temp.res$youden)) or be FisherZ(mean(temp.res$youden))
    var.youden.z <- FisherZ.Var(youden,var.youden)
    #############

    #############will not be outptuted since sample size calculation is based on normality assumption 
    sampleSize.CI <- quantile(as.numeric(temp.res$sampleSize),prob=c(alpha/2,1-alpha/2),na.rm=TRUE)


    #out <- data.frame(t.minus=all.mean["t.minus"],t.plus=all.mean["t.plus"],youden=youden,youden.z=youden.z,var.t.minus=all.var["t.minus"],var.t.plus=all.var["t.plus"],var.youden=var.youden,var.youden.z=var.youden.z,t.minus.lower=t.minus.CI[1],t.minus.upper=t.minus.CI[2],t.plus.lower=t.plus.CI[1],t.plus.upper=t.plus.CI[2],youden.lower=youden.CI[1],youden.upper=youden.CI[2],youden.z.lower=youden.z.CI[1],youden.z.upper=youden.z.CI[2],sampleSize=ceiling(all.mean["sampleSize"]),sampleSize.lower=ceiling(sampleSize.CI[1]),sampleSize.upper=ceiling(sampleSize.CI[2]),row.names = NULL)


    #return(list(bootResult=out,boot.iter.out=temp.res))

    ###bootResult:the averaged t.minus,t.plus,youden and sample size from bootstrap samples and associated variances and lower/upper 95% CI
    ###res.bootsample: NBOOT*6 matrix, record all results from bootstrap samples
    ###all.var: variance of t.minus   t.plus        Se        Sp        Sm    youden
    ##t.minus.CI, t.plus.CI,youden.CI: CI of alpha/2 to 1-alpha/2

    out <- list(var.youden=var.youden,var.t.minus=var.t.minus,var.t.plus=var.t.plus,var.youden.z=var.youden.z ,youden.CI=youden.CI,t.minus.CI=t.minus.CI,t.plus.CI=t.plus.CI,youden.z.CI=youden.z.CI)
    

    return(out)
  }

