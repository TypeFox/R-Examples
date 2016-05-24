### P7 Kernels

`p7` <-
function(t,mu,a,m) {
	(1+((t-mu)^2)/(m*a^2))^(-m)
	}

### Transformation


`partrans` <-
function(par,d,n) {
	partr<-rep(0,length(par))
	k<-(length(par)-2)/4
	partr[1]<-log((par[1]+d)/(d-par[1]))
	partr[2]<-par[2]
	lokvek<-1
	for (i in 1:k) lokvek<-c(lokvek,par[(i-1)*4+5])
	lokvek<-c(lokvek,n)
	lokvek<-log(diff(lokvek[2:length(lokvek)])/diff(lokvek[1:(length(lokvek)-1)]))

	for (i in 1:k) {
		partr[(i-1)*4+3]<-log(par[(i-1)*4+3])								#Höhe
		partr[(i-1)*4+4]<-log((par[(i-1)*4+4]-1)/(1000-par[(i-1)*4+4]))					#Shape
		partr[(i-1)*4+5]<-lokvek[i]									#Lok
		partr[(i-1)*4+6]<-log(par[(i-1)*4+6])								#Breite
	}
		partr
}

`partrans2` <-
function(par,d,n) {
	partr<-rep(0,length(par))
	k<-(length(par)-2)/4
	partr[1]<-log((par[1]+d)/(d-par[1]))
	partr[2]<-par[2]

	for (i in 1:k) {
		partr[(i-1)*4+3]<-log(par[(i-1)*4+3])								#H�he
		partr[(i-1)*4+4]<-log((par[(i-1)*4+4]-1)/(1000-par[(i-1)*4+4]))				#Shape
		partr[(i-1)*4+5]<-par[(i-1)*4+5]								#Lok
		partr[(i-1)*4+6]<-log(par[(i-1)*4+6])								#Breite
	}
		partr
}


`parbtrans` <-
function(partr,d,n) {
	par<-rep(0,length(partr))
	k<-(length(partr)-2)/4
	par[1]<-d*(exp(partr[1])-1)/(exp(partr[1])+1)	
	par[2]<-partr[2]
	lokvek<-NULL
	z<-1
	p<-1
	for (i in 1:k) {
		lokvek<-c(lokvek,exp(partr[(i-1)*4+5]))
		p<-p*lokvek[i]
		z<-z+p
		}
	lokvek2<-(n-1)/z
	lokvek3<-1+lokvek2
	if (k>1) for (i in 2:k) {
		lokvek2<-c(lokvek2,lokvek2[i-1]*lokvek[i-1])
		lokvek3<-c(lokvek3,lokvek3[i-1]+lokvek2[i])
	}
	for (i in 1:k) {
		par[(i-1)*4+3]<-exp(partr[(i-1)*4+3])
		par[(i-1)*4+4]<-(1000*exp(partr[(i-1)*4+4])+1)/(exp(partr[(i-1)*4+4])+1)		#shape
		par[(i-1)*4+5]<-lokvek3[i]
		par[(i-1)*4+6]<-exp(partr[(i-1)*4+6])
	}
		par
}


partransoutp<-function(par,left,h) {
	height<-par[1]
	m<-par[2]
	mu<-par[3]
	a<-par[4]
	lok<-left+(mu-1)*h
	fwhm<-2*a*sqrt( m*(2^(1/m)-1))*h
	intens<-beta(m-0.5,0.5)*sqrt(m)*a*height*h
	c(lok,height,intens,fwhm,m)
}




### Funktionen f?r MRB 

mr<-function(resid) .C("multires",as.double(resid),as.double(1),as.integer(length(resid)),PACKAGE="diffractometry")[[2]]

`mrqsim` <-
function(n,alpha,rep=10000) {
empvert<-rep(0,rep)
	for(i in 1:rep) {
	test<-rnorm(n,mean=0,sd=1)
	empvert[i]<-.C("multires",as.double(test),as.double(1),as.integer(length(test)),PACKAGE="diffractometry")[[2]]
	}
	thresh<-quantile(empvert,prob=1-alpha)
	thresh
}

`mrquant` <-
function(n,q,erg) {

	erg$tab[erg$sizevek==n,erg$quantvek==q]
	
}




### Anpassung f?r 1 Interval mit gegebener Anz. Kerne


`pkdecompint` <-
function(baslfit, intnum, k, thresh=0, alpha=0.1, heterosk=TRUE, maxiter=10000, dispers=1, baselim=c(0.05,5)){
  # Parameters: beta_0,beta_1,mass_1,m_1,mu_1,a_1,...,mass_k,m_k,mu_k,a_k
 	# Check whether there is at least one peak 
	if(length(baslfit$npks) == 0)
    stop("no peaks are found in data")  
  # Check whether parameter intnum is a reasonable number of peaks
  if(intnum > length(baslfit$npks))
    stop(c("value intnum must be between 1 and ", length(baslfit$npks)))
	y<-baslfit$pks[baslfit$indlsep[intnum]:baslfit$indrsep[intnum]]		# Daten werden extrahiert
	x<-baslfit$x[baslfit$indlsep[intnum]:baslfit$indrsep[intnum]]		# Daten werden extrahiert
	gwidth<-diff(x)[1]
	gewichte<-rep(1,length(y))																							# Gewichte werden mi 1 initialisiert (f�r den Fall w=F)
if (heterosk==FALSE) {			
	sig<-baslfit$spl$sigma
	stand2<-sqrt(baslfit$pmg$fn[baslfit$indlsep[intnum]:baslfit$indrsep[intnum]])
	for (i in 1:length(stand2)) if (stand2[i]<sig) stand2[i]<-sig		       						
	gewichte<-1/stand2
	}

if (heterosk==TRUE) {
	stand2<-baslfit$pmg$scl[baslfit$indlsep[intnum]:baslfit$indrsep[intnum]]
	gewichte<-1/stand2
	dispers<-1
}

basl<-(baslfit$baseline$basisl[baslfit$indlsep[intnum]:baslfit$indrsep[intnum]])		# Basislinie

versch<-baselim[1]*mean(basl)

																	bslope<-baselim[2]*gwidth
# Maximal allowed shift of baseline 

n<-length(y)

if (thresh==0) {
	thresh<-mrquant(n,(1-alpha),mr)
	if (length(thresh)!=1) thresh<-mrqsim(n,alpha)															}						# Falls keine Schranke angegeben, wird simuliert

auswertung<-function(para) {																										# Hier wird der C-Code aufgerufen, um f�r gegebene Parameter
	erg<-.C("p7fit",as.double(y),double(n),double(n),as.double(para),double(4*k+2),double(1),as.double(versch),	# die Anpassung und RSS auszurechnen
	as.integer(n),as.integer(k),as.integer(1),as.double(gewichte),PACKAGE="diffractometry")
	
	list(fit=erg[[2]], resid=erg[[3]], rss=erg[[6]])
	}

ausw2<-function(para){																															# Das selbe, es wird nur die RSS ausgeben (f�r die Optimierung)
	auswertung(para)[[3]]
	}

	lower<-partrans2( c( -0.75*versch,-bslope, rep(c(0.1,1.00001, log(1/(n-k-1)) ,1),k) )  ,versch,n)    # Schranken f�r den Suchraum
 	upper<-partrans2( c( 0.75*versch,bslope, rep(c(1.4*max(y),500, log(n-k-1) ,n/2),k) ) ,versch,n)    


fit.not.ok<-1
l<-0
best.par<-rep(0,4*k+2)
best.rss<-Inf

while(fit.not.ok==1 && l<maxiter) {																		# Suche solange, bis MRB erf�llt oder maxiter erreicht

l<-l+1

	candid<-runif(4*k+2,min=lower,max=upper) 			
	erg<-try(optim(par=candid,fn=ausw2, method="BFGS"))

	if (!inherits(erg,"try-error")) {

	if (erg$value < best.rss) {
	best.rss<-erg$value
	best.par<-erg$par

	erg2<-auswertung(erg$par)
	erg.fit<-erg2$fit
	erg.resid<-erg2$resid
	erg.rss<-erg2$rss
	erg.par<-erg$par

	if (.C("multires",as.double(erg.resid),as.double(1),as.integer(n),PACKAGE="diffractometry" )[[2]]<=dispers*thresh) fit.not.ok<-0

	}
	}

cat(paste("Interval: ",intnum,", Number of Peaks: ",k,", Iteration: ",l,"\n",sep=""))
}


parb<-parbtrans(erg.par,versch,n)


parmat<-NULL
for (g in 1:k) parmat<-rbind(parmat,partransoutp(parb[(3+(g-1)*4):(6+(g-1)*4)] ,baslfit$x[baslfit$indlsep[intnum]], gwidth) )
colnames(parmat)<-c("loc","height","intens","FWHM","m")


parbl<-c(parb[1],parb[2]/gwidth)

fitpk<-NULL
for (g in 1:k) fitpk<-rbind(fitpk,parb[(3+(g-1)*4)]*p7((1:n),parb[(5+(g-1)*4)],parb[(6+(g-1)*4)],parb[(4+(g-1)*4)]) )

baslchg<-parb[1]+(1:n)*parb[2]

list(intnr=intnum, x=x, y=y, fit=erg.fit, fitpk=fitpk, basl=basl, baslchg=baslchg, rss=erg.rss, num.ker=k, par=parb, parbl=parbl, parpks=parmat, accept=(fit.not.ok==0), alpha=alpha, thresh=thresh)
}

### Function to fit peaks for whole data set


`pkdecomp` <-
function(baslfit,intnum=0, alpha=0.1, maxiter1=500, maxiter=10000, hmax=5, maxsolutions=3,heterosk=TRUE,baselim=c(0.05,5),dispers=1) {
  # baselim[1] muss im Intervall (0,1) liegen, baselim[2] > 0
  if(!(baselim[1] > 0 && baselim[1] < 1)) {
    stop("baselim[1] must be between (0,1)")
  }
	if(baselim[2] < 0) {
    stop("baselim[2] must be greater than zero")
  }
  # alpha must be in {0.01,0.02,...,0.2}
  if(alpha > 0.01 && alpha < 0.2) 
      alpha = round(alpha, digits = 2)  
  if(!any(alpha == (1:20)/100)) {
    alpha <- 0.1 
    warning("alpha must be element of {0.01,0.02,...,0.2}! alpha is set to 0.1")
  }   
	geserg<-NULL
	j<-0
	 	# chekc whether there is at least one peak
	if(length(baslfit$npks) == 0)
    stop("no peaks are found in data")  
  # check whether intnum is a reasonable number of peaks
  if(intnum > length(baslfit$npks))
    stop(c("value intnum must be between 1 and ", length(baslfit$npks)))
	if (intnum==0) intnum<-1:length(baslfit$indlsep)
	if ((heterosk==FALSE) && (dispers == 0)) dispers<-mad(calcdisp(baslfit))
	for (z in intnum) {
	h=baslfit$npks[z]
	accept<-F
	thresh<-0
	while (accept==FALSE & h<=hmax) {	
		j<-j+1                                            
		if (h==1) tfout<-pkdecompint(baslfit,z,h,thresh,alpha,heterosk=heterosk,baselim=baselim,dispers=dispers,maxiter=maxiter1)
		if (h>1) tfout<-pkdecompint(baslfit,z,h,thresh,heterosk=heterosk,baselim=baselim,dispers=dispers,alpha,maxiter)	
		geserg[[j]]<-tfout
		h<-h+1
		accept<-tfout$accept
		thresh<-tfout$thresh
		}
		
		if (maxsolutions>1 & accept==TRUE) {		#prodcue maxsolutions solutions, if more than one is requested
		for (m in 2:maxsolutions) {
			if (h==1) tfout<-pkdecompint(baslfit,z,(h-1),thresh,alpha,heterosk=heterosk,baselim=baselim,dispers=dispers,maxiter=maxiter1)
			if (h>1) tfout<-pkdecompint(baslfit,z,(h-1),thresh,alpha,heterosk=heterosk,baselim=baselim,dispers=dispers,maxiter)	
			j<-j+1
			geserg[[j]]<-tfout
			}
		}
		
	}
geserg
}


### Calculate additional dispersion parameter if heterosk=FALSE and disp=0 

`calcdisp` <-
function(daten) {
	intnum<-1:length(daten$y)
	ind2<-NULL
	for (i in 1:length(daten$indlsep)) ind2<-c(ind2,(daten$indlsep[i]:daten$indrsep[i]))
	intnum<-intnum[-ind2]
	((daten$y-daten$baseline$basisl)/sqrt(daten$baseline$basisl))[intnum]	
}

.onLoad <- function(lib, pkg) {
  if(version$major<2)
    stop("This version for R 2.00 or later")
  library.dynam("diffractometry", pkg, lib)
} 

