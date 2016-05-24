bw.dboot2 <- function(y,sig, h0='dboot1', error='normal',B=1000,grid=100,ub=2)
  {
    ################################################################### 
    suppL2 <- function(h,x0,fx0,y,sig){
      n=length(y);n2=length(fx0);
	  # implement the smoothing bootstrap approach
	  fy0 = density(y, bw="nrd0")
	  y = rnorm(n, sample(y, size = n, replace = TRUE), fy0$bw)
      if(length(sig)==1){
        fx1 = DeconPdf(y,sig,error='normal',bw=h,
          fft=TRUE,n=n2,from=min(x0),to=max(x0));
      }else{
        fx1 = DeconPdf(y,sig,x=x0,error='normal',bw=h);
      }
      return(mean((fx1$y-fx0)^2));
    }
    suppboot <- function(h,x0,fx0,y,B=1000,sig){
      hs = as.matrix(h,ncol=B);
      result = apply(hs,1,suppL2,x0=x0,fx0=fx0,y=y,sig=sig);
      return(mean(result));
    }
    hsupp <- function(y,sig,grid=100,ub=2,B=1000){
	  fnaive = 	DeconPdf(y,sig,error='normal',bw=h0,
        fft=TRUE,from=min(y),to=max(y));
      x0 = fnaive$x; fx0=fnaive$y;
      hk = matrix(seq(h0/ub,ub*h0,length=grid),ncol=1);
      bw = apply(hk,1,suppboot,x0=x0,fx0=fx0,y=y,B=B,sig=sig);
      return(hk[which(bw==min(bw))][1]);
    }
    ###################################################################
    gaussL2 <- function(h,x0,fx0,y,sig){
      n=length(y);n2=length(fx0);
	  # implement the smoothing bootstrap approach
	  fy0 = density(y, bw="nrd0")
	  y = rnorm(n, sample(y, size = n, replace = TRUE), fy0$bw)
      fx1 = DeconPdf(y,sig,error='laplacian',bw=h,fft=TRUE,
        n=n2,from=min(x0),to=max(x0));
      return(mean((fx1$y-fx0)^2));
    }
    gaussboot <- function(h,x0,fx0,y,B=1000,sig){
      hs = as.matrix(h,ncol=B);
      result = apply(hs,1,gaussL2,x0=x0,fx0=fx0,y=y,sig=sig);
      return(mean(result));
    }
    hgauss <- function(y,sig,grid=100,ub=2,B=1000){
	  fnaive = 	DeconPdf(y,sig,error='laplacian',bw=h0,
        fft=TRUE,from=min(y),to=max(y));
      x0 = fnaive$x; fx0=fnaive$y;
      hk = matrix(seq(h0/ub,ub*h0,length=grid),ncol=1);
      bw = apply(hk,1,gaussboot,x0=x0,fx0=fx0,y=y,B=B,sig=sig);
      return(hk[which(bw==min(bw))][1]);
    }
    ###################################################################
    homo = FALSE;
    if(length(sig)==1) homo=TRUE;
    if(!homo){
      if(length(y)!=length(sig))stop("Different length of 'y' and the SD(s)!");
      if(any(is.na(y))||any(is.na(sig))){
        sele = (!is.na(y))&(!is.na(sig));
        y=y[sele];sig=sig[sele];
      }
    }else{
      if(is.na(sig)||is.null(sig)) stop("SD(s) can not be empty!");
    }
    if(length(y)<3){stop("The sample size is too small!");}else{n=length(y);}
    if(!is.numeric(grid)) grid=100;
    if(!is.numeric(ub)) ub=2;
    if(ub<0) ub=abs(ub);
    if(ub<1) up=1/ub;
    if(ub==1){warning("Invalid upper boundary of searching grid!");ub=2}
    if(homo){
      type=switch(substr(tolower(error),1,3),
        nor='norm',lap='lap');
    }else{
      type=switch(substr(tolower(error),1,3),
        nor='hnorm',lap='hlap');
    }
	if (is.character(h0)) {
        h0 <- switch(tolower(h0), dnrd = bw.dnrd(y, sig, error = error), 
            dmise = bw.dmise(y, sig, error = error), 
			dboot1 = bw.dboot1(y, sig, error = error), 
			stop("Unknown bandwidth rule!!!"))
    }
    result=switch(type,
      'norm' = hsupp(y,sig=sig,grid=grid,ub=ub,B=B),
      'lap' = hgauss(y,sig=sig,grid=grid,ub=ub,B=B),
      'hnorm' = print("Heteroscedastic normal errors are not supported yet!"),
      'hlap' = print("Heteroscedastic Laplacian errors are not supported yet!"),
      stop("This error type is not supported yet!")
      );
    return(result);
  }
