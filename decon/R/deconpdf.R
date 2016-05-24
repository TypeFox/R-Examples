DeconPdf <- 
function(y, sig, x, error = "normal", bw = "dboot1", adjust = 1,
         fft=FALSE, n = 512, from, to, cut =3, na.rm = FALSE, 
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
  error = match.arg(tolower(error),c("laplacian","snormal","normal","hlaplacian"))
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
    fft=FALSE;
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
    if(fft) warning("FFT method is not applicable on user-defined grids!");
    fft=FALSE;
  }
  if(substr(tolower(error2),1,3)=='sno'){
    if(any(sig>bw))stop("The small normal error model is not appropriate!");
  }
  
  if(fft){
    if(n<32)stop("FFT is not suggested if the number of points to be evaluated is too small (n<32).");
    if(n>2^21)stop("Number of points to be evaluated is too large (>2^21)");
    n2 = 2^ceiling(log2(n));
    if(n2 != n) stop("n must be a positive integer power of 2 within the range 2^2 <= n <= 2^21." );
    fhat=switch(substr(tolower(error2),1,3),
      sno=.Fortran("FFTGauss",as.double(y),as.integer(length(y)),
        as.double(min(x)),as.double(max(x)),as.double(bw),as.double(sig),
        as.double(x), y=as.double(x), as.integer(length(x)),PACKAGE='decon'),
      lap=.Fortran("FFTLaplace",as.double(y),as.integer(length(y)),
        as.double(min(x)),as.double(max(x)),as.double(bw),as.double(sig),
        as.double(x), y=as.double(x), as.integer(length(x)),PACKAGE='decon'),
      nor=.Fortran("FFTSupport",as.double(y),as.integer(length(y)),
        as.double(min(x)),as.double(max(x)),as.double(bw),as.double(sig),
        as.double(x), y=as.double(x), as.integer(length(x)),PACKAGE='decon'),
      stop("The specified  error type is not supported!"));
  }else{
    fhat=switch(substr(tolower(error2),1,3),
      nor=.C("DKESupport",as.double(y),as.integer(length(y)),
        y=as.double(x),as.integer(length(x)),
        as.double(bw),as.double(sig),as.integer(0),PACKAGE='decon'),
      sno=.C("DKEGauss",as.double(y),as.integer(length(y)),
        y=as.double(x),as.integer(length(x)),
        as.double(bw),as.double(sig),as.integer(0),PACKAGE='decon'),
      lap=.C("DKELaplace",as.double(y),as.integer(length(y)),
        y=as.double(x),as.integer(length(x)),
        as.double(bw),as.double(sig),as.integer(0),PACKAGE='decon'),
      hno=.C("DKEHSupport",as.double(y),as.integer(length(y)),
        y=as.double(x),as.integer(length(x)),
        as.double(bw),as.double(sig),as.integer(0),PACKAGE='decon'),
      stop("The specified  error type is not supported!"));
  }
  return(structure(list(x = x,y = fhat$y,bw = bw,n = N,
                        call = match.call(), data.name = name,
                        has.na = FALSE), class = "Decon"))
}

print.Decon <- function (x, digits = NULL, ...) 
{
    cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name, 
        " (", x$n, " obs.);", "\tBandwidth 'bw' = ", formatC(x$bw, 
            digits = digits), "\n\n", sep = "")
    print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
    invisible(x)
}

plot.Decon  <- 
function (x, main = NULL, xlab = NULL, ylab = "Density", type = "l", 
    zero.line = TRUE, ...) 
{
    if (is.null(xlab)) 
        xlab <- paste("N =", x$n, "  Bandwidth =", formatC(x$bw))
    if (is.null(main)) 
        main <- deparse(x$call)
    plot.default(x, main = main, xlab = xlab, ylab = ylab, type = type, 
        ...)
    if (zero.line) 
        abline(h = 0, lwd = 0.1, col = "gray")
    invisible(NULL)
}
