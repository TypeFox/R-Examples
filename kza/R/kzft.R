#===========================================================================#
# kzft.R				                                                            #
# Copyright (C) 2016 Brian Close                                            #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#
.packageName <- "kza"

coeff<-function(m, k){as.vector(polynom::polynomial(rep(1,m))^k/m^k)}

kzft <- function(x, f=0, m=1, k=1) 
{
	 	h <- floor  ((m-1)/2)
		t  <- ceiling((m-1)/2)
		x <- append(x, rep(NA,t*k), after = length(x) )
		x <- append(x, rep(NA,h*k), after = 0)
		
  	z <- as.complex(x)
		for(i in 1:k) {
  	  s     <- rep(0, length(z))
  	  count <- rep(0, length(z))
  	  
  	  count <- count + !is.na(z)
  	  z[is.na(z)] <- 0
  	  s <- as.complex(s) + z
  	  
  	  for(j in 1:h) {
  	      # add offset
  	      z   <- c(rep(NA, j), x[1:(length(x)-j)])
  	      z <- z * exp(complex(real = 0, imaginary = 2*pi*f*j))
  	      count <- count + !is.na(z)
  	      z[is.na(z)] <- 0
  	      s <- s + z
  	  }
			for(j in 1:t) {
  	      z   <- c(x[(j+1):length(x)], rep(NA, j))
  	      z <- z * exp(complex(real = 0, imaginary = -2*pi*f*j))
  	      count <- count + !is.na(z)
  	      z[is.na(z)] <- 0
  	      s <- s + z
  	  }
  	  z <- s/count
  	  x<-z
  	}
	  z <- z[(t*k+1):(length(z)-(h*k))]
	  return(z);
}	

kzp <- function(y, m=length(y), k=1)
{
	writeLines(c("This may take some time to complete."))
	M<-(m-1)*k+1
	n=length(y)
	z<-matrix(unlist(lapply(0:(m-1), function(i) kzft(y,f=i/m,m,k=k))),nrow=n, ncol=m)

	d<-apply(z,2,function(z) {(abs(z)^2)*M})
	a<-colMeans(d,na.rm=TRUE)
	a<-a[1:(m/2)]

	structure(list(
		periodogram = a,
		window=m,
		k=k,
		var=var(y),
		smooth_periodogram=NULL,
		smooth_method=NULL,
		call=match.call()
            ),
        class = "kzp")
}

kztp <- function(x, m, k, box=c(0,0.5,0,0.5))
{
    x<-as.vector(x)
    if (((m-1)*k+1)>length(x)) stop("The value of (m-1)*k+1 needs to be less then the length of the input data x")

	writeLines(c("This may take some time to complete."))
	M<-(m-1)*k+1
	n=length(x)
	z<-matrix(unlist(lapply(0:(m-1), function(i) kzft(x,f=i/m,m,k=k))),nrow=n, ncol=m)

	T<-dim(z)[1]
	m<-dim(z)[2]
	rp1=box[1]; rp2=box[2]; cp1=box[3]; cp2=box[4];
	rm1<-max(round(m*rp1),1)
    rm2<-round(m*rp2)
    cm1<-max(round(m*cp1),1)
    cm2<-round(m*cp2)
    delta.rm<-rm2-rm1+1
    delta.cm<-cm2-cm1+1
	zp<-array(NA,dim=c(delta.rm,delta.cm,T))
    for ( i in (1:delta.rm) ) for ( j in (1:delta.cm) ){
    	 zp[i,j,]<-z[,i+rm1-1]*z[,j+cm1-1]*Conj(z[,i+j+rm1+cm1-2])*m^2
    }               
    d<-rowMeans(zp,dims=2)
    return(d)
}

periodogram <- function(y) {
	fourier<-fft(y)
	
	magnitude<-Mod(fourier)
	 
	# extract the phase which is atan(Im(fourier)/Re(fourier))
	phase<-Arg(fourier)
	
	# select only first half of vectors
	magnitude_firsthalf <- magnitude[1:(length(magnitude)/2)]
	phase_firsthalf<-phase[1:(length(magnitude)/2)]
	 
	# generate x-axis
	x.axis <- 1:length(magnitude_firsthalf)/length(magnitude)
	
	p<-cbind(x.axis, magnitude_firsthalf)
	return(p)
}

plot.kzp <- function(x, ...)
{
	if (is.null(x$smooth_periodogram)) dz<-x$periodogram else dz<-x$smooth_periodogram
	omega<-(0:(length(x$periodogram)-1))/x$window
	plot(omega, dz-mean(dz), type="l", xlab="Frequency", ylab="")
}

summary.kzp <- function(object, digits = getOption("digits"), top=1, ...)
{
	cat(" Call:\n ")
	dput(object$call, control=NULL)

	M=object$window
	if (is.null(object$smooth_periodogram)) {	d<-object$periodogram } else { d<-object$smooth_periodogram }
	
	mlist<-rep(0,top)
	for (i in 1:top) {
		mlist[i]<-which.max(d)
		d[which.max(d)]=NA			
	}

   cat("\n Frequencies of interest:\n")
   print((mlist-1)/M, digits=digits, ...)

    cat("\n Periods of interest:\n")
    print(M/(mlist-1), digits=digits, ...)
    invisible(object)
}

transfer_function <- function(m, k, lamda=seq(-0.5,0.5,by=0.01), omega=0)
{
      lamda<-lamda*2*pi
      omega<-omega*2*pi

      N<-length(lamda)
      tf<-array(0,dim=c(N, m))

      for ( j in (1:m) ){
         tf[,j]<-exp(1i*(lamda-omega)*j)
      }

      tf1<-rep(0,N)
      for ( i in (1:N) ){
         tf1[i]<-sum(tf[i,])
      }

      tf2<-(1/m)*tf1
      tf2<-abs(tf2)^k
      return(tf2)
}

nonlinearity.kzp<-function(x)
{
    n<-length(x)

    s<-rep(0,n)
    for (t in 2:(n-1)) {
        s[t]<-abs(x[t+1]-2*x[t]+x[t-1])
    }

    sq<-array(0, dim=c(n, n))

    for (i in (1:n)) for (j in (2:n)) {
        sq[i,j]<-sum(s[(max(1,(i-j+1))):(min(n,(i+j-1)))])
    }
    return(list(total=sum(s), matrix=sq))
}

variation.kzp<-function(x)
{
    n<-length(x)
    s=c(diff(x)^2,0)

    q<-array(0, dim=c(n, n))
    for (i in (1:n)) for (j in (2:n)) {
        q[i,j]<-sum(s[(max(1,(i-j+1))):(min(n,(i+j-2)))])
    }
    return(list(total=sum(s), matrix=q))
}

smooth.kzp<-function(object, log=TRUE, smooth_level=0.05, method = "DZ")
{
	if (class(object)!='kzp') stop ("Object type needs to be kzp.")
    n<-length(object$periodogram)
    spg<-rep(0,n)
    m<-rep(0,n)

	if (log==TRUE) p=log(object$periodogram) else p=object$periodogram
    if (method == "DZ") q<-variation.kzp(p)
    else if (method == "NZ") q<-nonlinearity.kzp(p)

    cc<-smooth_level*q$total

    for ( i in (1:n) ) {
        m[i]<-sum(q$matrix[i,1:n]<=cc)

        spg[i]<-mean(p[(max(1,(i-m[i]+1))):(min(n,(i+m[i]-1)))])
    }
    
    object$smooth_periodogram<-spg
    object$smooth_method=method
    return(object)
}
