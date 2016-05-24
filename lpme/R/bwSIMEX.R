"bwSIMEX" <- function(Y, W, method="HZ", sig, error="laplace", k_fold=2, B=10, 
                      h1=NULL, h2=NULL, length.h=10, Wdiff=NULL,
                      data=sys.frame(sys.parent()), na.action=na.fail, work.dir=NULL)
UseMethod("bwSIMEX")

"bwSIMEX.default" <- 
  function (Y, W, method="HZ", sig, error="laplace", k_fold=2, B=10,
            h1=NULL, h2=NULL, length.h=10, Wdiff=NULL,
            data=sys.frame(sys.parent()),
            na.action=na.fail, 
            work.dir=NULL) {

    #########################################################################################
    # data structure
    #########################################################################################
    n = length(Y);
    dd = cbind(Y,W,Wdiff);
    newd = dd[sample(1:n, n),];
    Y = newd[,"Y"];
    W = newd[,"W"];
    
    ## Fourier transform settings
    m  = 2^16; m_mid = m/2+1;
    beta  = sqrt((2*pi)/m); 
    beta2  = (2*pi)/(m*beta); 
    input  = seq(-m*beta/2, m*beta/2-beta, by=beta)
    output= seq(-pi/beta, pi/beta-beta2, by=beta2)
    mconst= (-1)^(0:(m-1));
    
    ## Kernel for which CF is (1-t^2)^8 with Normal errors
    FKsup  = function(t) { ifelse( (t<=1 & t>=-1), (1-t^2)^8, 0) }
    FKoutput = rep(0, m);
    FKoutput[m_mid] = FKsup(output[m_mid]);
    i=1; indicator1 =1; indicator2 =1;
    while ( ( abs(indicator1)>1e-30 | abs(indicator2)>1e-30 ) & (i<(m/2)) ) {
      indicator1 = FKsup(output[m_mid-i]);
      indicator2 = FKsup(output[m_mid+i]);
      FKoutput[m_mid-i] = indicator1;
      FKoutput[m_mid+i] = indicator2;
      i=i+1;
    }
    ## inverse FFT to get K
    Kinput  = Re( (-1)^(0:(m-1))/beta*fft( (-1)^(0:(m-1))*FKoutput)/m );
    dt = 0.0001;
    tt = seq(-1,1,dt);
    
    #### SIMEX settings
    cumfold = seq(0, n, ceiling(n/k_fold));
    if(length(cumfold)==k_fold) cumfold=c(cumfold,n);
    Ws = matrix(0, n, B); Wss=Ws;
    for (i in 1:B){
      if(is.null(Wdiff)){
        if(error=="laplace"){
          Ws[,i] = W+.rlaplace(n,0,sig/sqrt(2));
          Wss[,i] = Ws[,i]+.rlaplace(n,0,sig/sqrt(2));
        }else{
          Ws[,i] = W+rnorm(n,0,sig);
          Wss[,i] = Ws[,i]+rnorm(n,0,sig);
        }
      }else{
        Wdiff = newd[,"Wdiff"];
        Ws[,i] = W + Wdiff[sample(1:n,n,replace=TRUE)];
        Wss[,i] = Ws[,i] + Wdiff[sample(1:n,n,replace=TRUE)];
      }
    }
    bw1=.dMISE(W, sig=sig, error=error);
    bw2=mean(apply(Ws, 2, .dMISE, sig=sig, error=error));
    if(is.null(h1)) h1 = seq(bw1-bw1/2, bw1+bw1/2, length.out=length.h)
    if(is.null(h2)) h2 = seq(bw2-bw2/2, bw2+bw2/2, length.out=length.h)
    
    ## density of W
    dd=density(W)
    dens=splinefun(dd$x, dd$y)
    fWhat=function(w) {
      resp=dens(w)
      resp[resp<0]=0
      resp
    }
    pW = fWhat(W);
    ## density of Ws
    dd = density(as.vector(Ws));
    dens = splinefun(dd$x, dd$y);
    fWhats = function(w) {resp = dens(w); resp[resp<0]=0; resp};
    pWs= matrix(1, n, B);
    for (b in 1:B){
      pWs[,b] = fWhats(Ws[,b]);
    }
    
    #########################################################################################
    # change working directory (if requested..)
    #########################################################################################
    if(!is.null(work.dir))
    {
      cat("\n Changing working directory to ",work.dir,"\n")
      old.dir <- getwd()  # by default work in current working directory
      setwd(work.dir)
    }
    
    #########################################################################################
    # calling the c++ code and # output
    #########################################################################################
    if(method=="HZ"){
      if(is.null(sig)) stop("please specify sig if non naive method is used")
      if(error=="laplace"){
        foo <- .Call("SIMEXnewLap", 
                     input_ = input, 
                     output_ = output, 
                     beta_ = beta, 
                     beta2_ = beta2, 
                     mconst_ = mconst, 
                     Kinput_ = Kinput, 
                     W_ = W, 
                     Y_ = Y, 
                     Ws_ = Ws, 
                     Wss_ = Wss,
                     h1_ = h1,
                     h2_ = h2,
                     sigU_ = sig, 
                     cumfold_ = cumfold,
                     pW_ = pW,
                     PWs_ = pWs,
                     PACKAGE = "lpme");
      } else {
        foo <- .Call("SIMEXnewGau", 
                     input_ = input, 
                     output_ = output, 
                     beta_ = beta, 
                     beta2_ = beta2, 
                     mconst_ = mconst, 
                     Kinput_ = Kinput, 
                     W_ = W, 
                     Y_ = Y, 
                     Ws_ = Ws, 
                     Wss_ = Wss,
                     h1_ = h1,
                     h2_ = h2,
                     sigU_ = sig, 
                     cumfold_ = cumfold,
                     pW_ = pW,
                     PWs_ = pWs,
                     PACKAGE = "lpme");
      } 
    } else if (method=="DFC") {
      if(is.null(sig)) stop("please specify sig if non naive method is used")
      if(error=="laplace"){
        foo <- .Call("SIMEXjasaLap", 
                     W_ = W, 
                     Y_ = Y, 
                     Ws_ = Ws, 
                     Wss_ = Wss,
                     h1_ = h1,
                     h2_ = h2,
                     sigU_ = sig, 
                     cumfold_ = cumfold,
                     pW_ = pW,
                     PWs_ = pWs,
                     dt_ = dt, 
                     t_ = tt,
                     PACKAGE = "lpme");
      } else {
        foo <- .Call("SIMEXjasaGau", 
                     W_ = W, 
                     Y_ = Y, 
                     Ws_ = Ws, 
                     Wss_ = Wss,
                     h1_ = h1,
                     h2_ = h2,
                     sigU_ = sig, 
                     cumfold_ = cumfold,
                     pW_ = pW,
                     PWs_ = pWs,
                     dt_ = dt, 
                     t_ = tt,
                     PACKAGE = "lpme");
      }
    } else {
      stop("please specify method to be HZ or DFC")
    }
    
    #########################################################################################
    # save state
    #########################################################################################
    if(!is.null(work.dir))
    {
      cat("\n\n Changing working directory back to ",old.dir,"\n")
      setwd(old.dir)
    }
    model.name <- "local polynomial estimation";
    
    #### Save to a list
    h1 = (foo$h1); CVh1 = foo$CVh1; 
    h2 = (foo$h2); CVh2 = foo$CVh2; 
    hwNEW = (h1[which.min(CVh1)])^2/h2[which.min(CVh2)];
    output <- list(modelname=model.name,
                   bw = hwNEW,
                   h1=foo$h1,
                   CVh1=foo$CVh1,
                   h2=foo$h2,
                   CVh2=foo$CVh2);
    class(output) <- c("bwSIMEX")
    output
  }