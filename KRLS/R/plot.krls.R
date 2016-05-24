plot.krls <-
function(x,
         which=c(1:2),
         main="distributions of pointwise marginal effects",
         setx="mean",
         ask = prod(par("mfcol")) < nplots,
         nvalues = 50,
         probs = c(.25,.75),
         ...)
      {
            
        if( class(x)!= "krls" ){
        warning("x not of class 'krls'")
        UseMethod("summary")
        return(invisible(NULL))
        }
     
        d <- ncol(x$X)
        n <- nrow(x$X)
        if(length(probs)!=2){
          stop("length(probs) must be 2")
        }
        
       # check setx
        if(is.numeric(setx)){
          if(length(setx)!=d){
           stop("length(setx) must be equal to number of predictors") 
          }
         } else {
           if(length(setx)!=1){stop("setx must be one of mean or median")}
           if(sum(setx %in% c("mean","median"))<1){stop("setx must be one of mean or median")}
           setx <- apply(x$X,2,setx)
          }
        
nplots <- 0
if(1 %in% which){ nplots <- nplots + 1}
if(2 %in% which){ nplots <- nplots + d}        
   
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
 
        if(is.null(colnames(x$X))){
          colnames(x$X) <- paste("x",1:d,sep="") 
          } 
        
        # derivatives
        if(is.null(x$derivatives)){
          cat("recompute krls x with krls(...,derivative = TRUE) to plot marginal effects\n")
         } else {
          colnames(x$derivatives) <- colnames(x$X)
         
          if(1 %in% which){ # histograms of partial derivatives
          form <-  as.formula(paste("~",paste(colnames(x$derivatives),collapse="+"),sep=""))
          require(lattice)
          
          print(histogram(form,
                    data=data.frame(x$derivatives),
                    breaks=NULL,
                    main=main
                    ,...)
                )
          #if(length(which)!=1){readline("Press any key for next plot")}
          }
         }
          
          if(2 %in% which){  # conditional expectation plots
             lengthunique    <- function(x){length(unique(x))}
             # vector with positions of binary variables
             binaryindicator <- which(apply(x$X,2,lengthunique)==2)
             quantiles <-  apply(x$X,2,quantile,probs=probs)     
             
              for(i in 1:d){
                
                if(i %in% binaryindicator){ # E[Y|X] for binary Xs
                  Xi <- c(min(x$X[,i]),max(x$X[,i]))
                  Newdata <- matrix(rep(setx,2),ncol=d,byrow=T)
                  Newdata[,i] <- Xi
                  
                } else {
                # E[Y|X] plots for cont Xs
                Xi <- seq(quantiles[1,i],quantiles[2,i],length.out=nvalues)
                Newdata <- matrix(rep(setx,nvalues),ncol=d,byrow=T)
                Newdata[,i] <- Xi
                }
                pout      <- predict(x,newdata=Newdata,se=TRUE)
                Ylo <- pout$fit-1.96*pout$se
                Yhi <- pout$fit+1.96*pout$se
                plot(y=pout$fit,x=Xi,
                     xlab=colnames(x$X)[i],
                     ylab=c("E[Y|X]"),
                     ylim=c(min(Ylo) -.25*sqrt(var(pout$fit)),
                            max(Yhi))+.25*sqrt(var(pout$fit)),pch=19
                     )
                arrows(x0=Xi,y0=Ylo,y1=Yhi,length = 0)
                
              }
          }  
          
}





