diffractogram<-function(data, tau=2.5,gam=1, scl.factor=1.2, maxwdth=5, intnum=0, alpha=0.1, maxiter1=500, maxiter=10000, hmax=5, maxsolutions=3,heterosk=TRUE,baselim=c(0.05,5),dispers=1) {
  #Abfrage äquidistante Winkelwerte für 2Theta
  if(any(diff(data[,1], 1, 2) > 10^-10)) 
    stop("2 theta values have to be equidistant")
  #Abfrage ob NAs vorhanden
  if(any(is.na(data))) 
    stop("data contains NAs")
  #Abfrage ob negative Werte vorhanden
  if(any(data[,2] < 0)) 
    stop("negative counts are not permitted")
  #Tau muss im Interall [2,4] liegen
  if(!(tau >= 2 && tau <= 4)) 
    stop("tau must be between [2,4]")
  # baselim[1] muss im Intervall (0,1) liegen, baselim[2] > 0
  if(!(baselim[1] > 0 && baselim[1] < 1)) 
    stop("baselim[1] must be between (0,1)")
	if(baselim[2] < 0)
    stop("baselim[2] musst be greater zero")
  # alpha darf nur aus der Menge {0.01,0.02,...,0.2} stammen
  if(alpha > 0.01 && alpha < 0.2) 
      alpha = round(alpha, digits = 2)  
  if(!any(alpha == (1:20)/100)) {
    alpha <- 0.1 
    warning("alpha must be element of {0.01,0.02,...,0.2}! alpha is set to 0.1")
  }     
  basl<-baselinefit(data, tau,gam, scl.factor, maxwdth)
	pks=NULL
	if (is.list(basl)) pks<-pkdecomp(basl,intnum, alpha, maxiter1, maxiter, hmax, maxsolutions,heterosk,baselim,dispers)
	list(basl=basl,pks=pks)
}
 
