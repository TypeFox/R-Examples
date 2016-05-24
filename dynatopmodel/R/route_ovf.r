# use the eigenvector approach to solving a system of linear first order ODEs
dist.eigen <- function(groups, A, ex, vof, dt, ichan=1, eig=NULL, tm=NULL, max.t=NULL)
{
	# A is the scaled weighting matrix for this discretisation
	ex.hru <- ex  # ex[-ichan]
	if(any(ex.hru > 0, na.rm=T))
	{
        # total excess storage
      #  ex.hru <- ex.hru

		# vof a vector of overflow flow
		# relate storage to discharge by q = vs
		# then s' = t(q) %*% w-q = v*(t(w)%*% s-s)
	#	ex.hru[ichan]<- 0
#		eig <- get_routing_eig(A, vof)
		ex.dtt <- ode(y=ex.hru,
                      times=c(0, dt), func=dsdt, parms=list(A=A, vof=vof))
		ex.hru <- ex.dtt[2,-1]  # exclude first column with times and take result at time dt, in final row


#		Av <- vof * (A - identity.matrix(length(vof)))
#		eig <- eigen(Av)

		# scale to give rate of change of storage in tewrms of storgage

	#	vof[ichan] <- 0 #vof[-ichan]
#		vof <- matrix(rep(vof), col=nrow(w))
# 		A <- t(whru) -identity.matrix(nrow(whru))
#
# 		A <- vof*A
# 		eig <- eigen(A)
		# eigenvalues
#		lambda <- eig$values
#		u <- eig$vectors

		# the generral asolution to this system is
		# s = c1*exp(l1*t)u + c2*exp(l2*t*u)+...

		# where the ci are arbitrary constants


		# Assume that ex is the initial state of the system so, setting t=0
		# ex = c[i]u[i]=t(c)%*%u=t(u)%c

		# solve for coefficients using initial conditions
#		ci <- solve(u, ex)

		# factors at t=dt
#		cdt <- ci * exp(lambda*dt)

#		ex.hru <- u %*% cdt

		# convert back to specific storage

		# the channel storages are unchanged
#		ex[-ichan] <- ex.hru
# specific storages
  #  ex.hru<- 0
		ex <- ex.hru/groups$area


	}

	# convert storage back to equivalent flow on return
	return(ex)
}

dsdt <- function(t, st, parms, ...)
{
	A <- parms$A
	vof <- parms$vof
    # convert to a discharge
	Av <- vof * (A - identity.matrix(length(vof)))
  dsdt <- vof*st %*% t(A)-vof*st
  # Av %*% s
	return(list(dsdt))

}



# show_ovf <- function(groups, tm, max.t, lambda, u, ci)
# {
# 	tms <- seq(0, max.t, length.out=6)
#
# 	qovf <- t(sapply(tms,
# 									 function(dt)
# 									 {
# 									 	cdt <-  ci * exp(lambda*dt)
# 									 	ex.hru <- u %*% cdt
# 									 }
# 	))
#
# 	s.ovf <- t(apply(round(qovf, 3), MARGIN=1, function(x)x/groups$area))
#
#
# 	s.ovf <- cbind(s.ovf, rowSums(qovf)/sum(groups$area))
# 	groups$id <- paste0("HRU", groups$id)
# 	groups[1,"id"]<- "Reach 1"
# 	colnames(s.ovf)<- c(groups$id, "Overall")
#
# 	n.breaks <- 50
# 	cols <- colorRampPalette(c("white", "blue"))(n.breaks)
# 	cols.hru <- cut(s.ovf, breaks=n.breaks, labels=F)
#
# 	its <- xts(order.by=tm+tms*3600, matrix(cols.hru, nrow=length(tms)))
# 	ts.ex <- xts(round(s.ovf*1000, 2), order.by=index(its))
#
# 	hru <- disc$hru #  drn*100+
# 	drn <- gBuffer(gwy$drn, w=3)
#
# 	write.zoo(ts.ex, "e:/junk/dtm/ovf.tsv", sep="\t", quote=F)
#
#
# 	bvals <-seq(min(s.ovf), max(s.ovf), length.out=n.breaks)
# #	cols[unique(its)]
#
# 	bvals <- round(1000*bvals[unique(its)], 2)
# #	hru[hru>200]<-1
# 	hru <- subs(hru, data.frame(unique(hru), 1:length(unique(hru))))
# 	sel <- extent(281400, 282000, 285650, 286082)
#
# browser()
#
# fn <- "e:/junk/dtm/ovf.jpg"
# #	jpeg(fn)
# #on.exit(dev.off())
#
# 	par(family="serif")
# 	layout(matrix(c(1:6, 7, 7), ncol=2, byrow=T), heights=c(0.5,0.5,0.5,0.1))
#
# par(mar=c(2,2.5,3,0))
# 	for(i in 1:length(tms))
# 	{
#
#
# 		sp::plot(hru, col=cols[its[i,-1]], ext=sel, legend=F, main=index(its)[i], cex.main=0.75)
# 		sp::plot(drn, col=cols[its[i,1]], add=T, border=NA)
# 		raster::contour(gwy$dem, add=T, col=make.transparent("brown"), nlevels=20)
#
#
# 	#	readline("")
#
# 	}

# plot.new()
#
# 	legend(x="bottomleft", fill=cols[sort(unique(its))], legend=sort(bvals), horiz=T)
#
# 	dev.off()
# }

# method of solution by eigen values and vectors
get_routing_eig <- function(A, vof, ichan=1)
{

	# scale to give rate of change of storage in terms of storage
	nhru <- nrow(A)

	#	vof[ichan] <- 0 #vof[-ichan]
	vof.mat <- matrix(rep(vof, nhru), nrow=nhru, byrow=T)
	Av <- A - diag(vof)

	eig <- eigen(Av)

	return(eig)

}
