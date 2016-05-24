DeconNpr <- 
function(y, sig, z, x, error = "normal", bw = "dboot1", adjust = 1,
         n = 512, from, to, cut =3, na.rm = FALSE,  
         grid=100,ub=2,...) 
{
  CheckValidity <- function(y,na.rm){
    if (!is.numeric(y)) 
      stop("argument 'y' must be numeric")
    name <- deparse(substitute(y))
    y <- as.vector(y)
    y.na <- is.na(y)
    if (any(y.na)) {
      if (na.rm)  y <- y[!y.na]
      else stop("'y' contains missing values")
    }
    y.finite <- is.finite(y)
    if (any(!y.finite)) {
      y <- y[y.finite]
    }
    ny <- length(y)
    list(y=y,ny=ny,name=name);
  }
  if (length(list(...)) > 0) 
    warning("non-matched further arguments are disregarded")
  error = match.arg(tolower(error),c("laplacian","snormal","normal"))
  
  kernel="normal";
  if(error=="normal")kernel="support";

  if(length(sig)==1){
    error1 = switch(error,
      laplacian="laplacian",
      snormal="normal",
      normal="normal");
    error2 = switch(error,
      laplacian="laplacian",
      snormal="snormal",
      normal="normal");
  }else{
    error1 = switch(error,
      laplacian="laplacian",
      snormal="normal",
      normal="normal");
    error2 = switch(error,
      laplacian="hlaplacian",
      snormal="hsnormal",
      normal="hnormal");
    }
  yout = CheckValidity(y);
  y=yout$y; ny=yout$ny; name=yout$name; N = ny;
  zout = CheckValidity(z);
  z=zout$y; nz=zout$ny;
  if(ny!=nz)stop("'z' and 'y' have different lengths!");
  if(any(sig<=0))stop("Standard deviations should be positive.");
  if(error=="hlaplacian"&length(sig)!=length(y))
    stop("'sig' and 'y' have different length.");
  if(!is.numeric(n)|n<3)
    stop(paste("n must be a positive integer power of 2 within the range 2^2 <= n <= 2^21." ));
  if (is.character(bw)) {
    if (ny < 2) 
      stop("need at least 2 points to select a bandwidth automatically")
    bw <- switch(tolower(bw),
	 			 dnrd = bw.dnrd(y,sig,error=error1),
	 			 dmise = bw.dmise(y,sig,error=error1,kernel=kernel,grid=grid,ub=ub),
    			 dboot1 = bw.dboot1(y,sig,error=error1,grid=grid,ub=ub),
    			 dboot2 = bw.dboot2(y,sig,error=error1,grid=grid,ub=ub),
                 stop("Unknown bandwidth rule!!!"))
  }
  if (!is.finite(bw)) 
    stop("non-finite 'bw'")
  bw <- adjust * bw
  if (bw <= 0) 
    stop("'bw' is not positive.")
  if(missing(x)){
    if (missing(from)) 
      from <- min(y) - cut * bw
    if (missing(to)) 
      to <- max(y) + cut * bw
    if (!is.finite(from)) 
      stop("non-finite 'from'")
    if (!is.finite(to)) 
      stop("non-finite 'to'")
    if(from>=to){
      stop("'from' is not smaller than 'to'!");
    }else{
      x=seq(from,to,length=n);
    };
  }else{
    x=sort(as.vector(x));
    xout = CheckValidity(x);
    x=xout$y; n=xout$ny;
  }
  if(substr(tolower(error2),1,3)=='sno'){
    if(any(sig>bw)){
      warning("The small normal error model is not appropriate!");
      sele = sig>bw;
      sig[sele] = .95*bw[sele];
    }
  }
  fhat=switch(substr(tolower(error2),1,3),
    nor=.C("NPRSupport",as.double(y),as.integer(length(y)),
      as.double(z),y=as.double(x),as.integer(length(x)),
      as.double(bw),as.double(sig)),
    sno=.C("NPRGauss",as.double(y),as.integer(length(y)),
      as.double(z),y=as.double(x),as.integer(length(x)),
      as.double(bw),as.double(sig)),
    lap=.C("NPRLaplace",as.double(y),as.integer(length(y)),
      as.double(z),y=as.double(x),as.integer(length(x)),
      as.double(bw),as.double(sig)),
    hno=.C("NPRHSupport",as.double(y),as.integer(length(y)),
      as.double(z),y=as.double(x),as.integer(length(x)),
      as.double(bw),as.double(sig)),
    stop("The specified  error type is not supported!"));
  return(structure(list(x = x,y = fhat$y,bw = bw,n = N,
                        call = match.call(), data.name = name,
                        has.na = FALSE), class = "Decon"))
}


