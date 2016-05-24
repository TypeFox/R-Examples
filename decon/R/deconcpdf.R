DeconCPdf <- 
function(y, sig, y0, error = "normal", bw1 = "dboot1", bw2="nrd0",adjust = 1,
         fft=FALSE, n = 512, from, to, cut =3, na.rm = FALSE, 
         grid=100,ub=2,tol=0, ...) 
  {
	if (length(sig)!=1) 
	  stop("Conditional density estimation is only the case of homoscedastic error!")
	if ((tol < 0) | (tol > 0.05))
		  stop("The tolerance value has to be greater than or equal to zero and less than 0.05!")
	type=switch(substr(tolower(error),1,3),
      nor = 1,lap = 2,stop("This error type is not supported yet!"));
	name <- deparse(substitute(y))
	## estimate the dke
	N <- length(y)
	fx <- DeconPdf(y,sig, error = error, bw = bw1, adjust = adjust,
	         fft=fft, n = n, from=from, to=to, cut = cut, na.rm = na.rm, 
	         grid=grid, ub=ub,...);
	## esitmate the mariginal density of y
    if (is.character(bw2)) {
        if (N < 2) 
            stop("need at least 2 points to select a bandwidth automatically")
        bw2 <- switch(tolower(bw2), nrd0 = bw.nrd0(y), nrd = bw.nrd(y), 
            ucv = bw.ucv(y), bcv = bw.bcv(y), sj = , `sj-ste` = bw.SJ(y, 
                method = "ste"), `sj-dpi` = bw.SJ(y, method = "dpi"), 
            stop("unknown bandwidth rule"))
    }
    if (!is.finite(bw2)) 
        stop("non-finite 'bw'")
    bw2 <- adjust * bw2
	fy <- density(y, from=y0-cut*bw2,to=y0+cut*bw2,n=1001,bw=bw2)
	fy0 <- fy$y[501]
	dlap <- function(x, sig) {exp(-1*abs(x)/sig)/(2*sig)}
	if (type==1) {fy.x <- dnorm(y0-fx$x, sd=sig)
		}else{
		if (type==2) {fy.x <- dlap(y0-fx$x, sig=sig)}	
		}
	if (tol>0) {fy.x[fy.x<tol] <- tol}
	fx.y <- fy.x*fx$y/fy0
    structure(list(x = fx$x, y = fx.y, bw= paste(round(fx$bw,5),"and", round(bw2,5), sep=" "),
        n = N, call = match.call(), data.name = name, has.na = FALSE), 
        class = "Decon")
  }

