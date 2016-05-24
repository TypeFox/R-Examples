flpkspace <-
function(theta,x,t,w,
			k	=2,
			m1=as.integer(nrow(x)/2),
			m2=as.integer(nrow(x)-1),
			etas.params,
			etas.l,
			mh=1,
#			kern.var=FALSE,
			iprint=FALSE,
			indweight=TRUE)     {
		
			kern.var=FALSE
# function called (for optimization) by flp.etas.nlm
# essentially is an interface to the FORTRAN  subroutine deltafl1kspace
# in this version only for k=2  (x,y)

iprint=FALSE
if (iprint) cat("flpkspace theta", theta,"\n")

cat("-")	
		n	=as.integer(nrow(x))
		dens	=array(0,n)
		 s	=matrix(c(1,0,0,1),2,2)
 rangex=t(apply(x,2,range))
indanis	=as.integer(kern.var)
nh	=k*(1+indanis) 
h.start=numeric(0)
h.start[1:nh]=as.double(exp(theta[1:nh]))

	ris.fl<-.Fortran("deltafl1kspacevar",
		x	=as.double(x),
		t	=as.double(t),
		w	=as.double(w),
		n	=as.integer(n),
		k	=as.integer(k),
		m1	=as.integer(m1),
		m2	=as.integer(m2),
		nh	=as.integer(nh),
		rangex	=as.double(rangex),
		h	=as.double(h.start),
		hdef	=as.double(h.start),
		dens	=as.double(dens),
		integr	=as.double(dens),
		delta	=as.double(dens),
		indanis	=as.integer(indanis),
		expweight=as.double(-0.2),
		indweight=as.integer(indweight),
		allocationerr=as.integer(0),
		NAOK=TRUE
				)
		if (ris.fl$allocationerr)cat("memory allocation failed in flp step","\n")		
		lambda	=etas.params[1]
		val	=-sum(log(lambda*ris.fl$dens[m1:m2]+etas.l[(m1+1):(m2+1)])
		-lambda*ris.fl$integr[m1:m2])
if (iprint) {
cat("Function Value ",val,"\n")
cat("  lambda: ",lambda,"\n")
cat("summary(ris.fl$dens[m1:m2])","\n")
cat(summary(ris.fl$dens[m1:m2]),"\n")
cat("summary(etas.l[(m1+1):(m2+1)])","\n")
cat(summary(etas.l[(m1+1):(m2+1)]),"\n")
cat("summary(ris.fl$integr[m1:m2])","\n")
cat(summary(ris.fl$integr[m1:m2]),"\n")

            }
                attr(val, "dens") <- ris.fl$dens
                attr(val, "delta") <- ris.fl$delta
                attr(val, "integr") <- ris.fl$integr
                attr(val, "hdef") <- ris.fl$hdef
                attr(val, "h") <- ris.fl$h
		return(val)
 }
