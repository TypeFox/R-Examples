

neg3t = function(x){
	
	n <- dim(x)[1]
	m <- dim(x)[2]
	k <- dim(x)[3]
	
#c1 <- classcrit(x)	
	c2 <- classcrit(x[n:1,,]) + classcrit(x[,m:1,]) + classcrit(x[,,k:1])
	
	return(c2/3/iccrit(x))
	
}



JBCI <- function(x,r = 1){
	stopifnot(length(dim(x)) == 3)
	n <- dim(x)[1]
	m <- dim(x)[2]
	k <- dim(x)[3]
	
#c1 <- classcrit(x)	
	c2 <- classcrit(x[n:1,,]) + classcrit(x[,m:1,]) + classcrit(x[,,k:1])
	ix <- iccrit(apply(x,r,sum))
	x2 <- apply(x,c(1:3)[-r],sum)
	
	return(c2/ix/( classcrit(x2)+2*classcrit(x2,FALSE) )*sum(x)^2)
	
}


neg3jb = function(x,r = 1){
	
	n <- dim(x)[1]
	m <- dim(x)[2]
	k <- dim(x)[3]
	
	x2 <- apply(x,c(1:3)[-r],sum)
	v1<-3*iccrit(x2)/( classcrit(x2)+2*classcrit(x2,FALSE) )*neg3t(x)
	
	return(v1)
	
}


CBCI = function(x,r = 1, joint.order = FALSE){
	stopifnot(length(dim(x)) == 3)
	n <- dim(x)[1]
	m <- dim(x)[2]
	k <- dim(x)[3]
	if(joint.order){
		x2 <- optile(x,iter=100,method="joint")
	}else{
		x2 <- x
	}
	xz <- apply(x2,c(c(1:3)[-r][1],r),sum)
	yz <- apply(x2,c(c(1:3)[-r][2],r),sum)
	
	bxz <- classcrit(xz)
	byz <- classcrit(yz)
	
	nxz <- bxz + classcrit(xz, FALSE)
	nyz <- byz + classcrit(yz, FALSE)
	
	iz <- iccrit(apply(x2,r,sum))
	
	v1 <- (nxz*nyz)/iz - (bxz*byz)/iz
	
	#v1<- neg3t(x)*3*iccrit(x)/v1
	v1 <- BCC(x)/v1
	return(v1)
	
}



# CBCI = function(x,vars = c(1,2), joint.order = FALSE){
	# stopifnot(length(vars) == 2)
	# nd <- length(dim(x))
	# stopifnot(all(vars < nd+1))
	
	# n <- dim(x)[1]
	# m <- dim(x)[2]
	# k <- dim(x)[3]
	# if(joint.order){
		# x2 <- optile(x,iter=100,method="joint")
	# }else{
		# x2 <- x
	# }
	# xz <- apply(x2,c(c(1:3)[-r][1],r),sum)
	# yz <- apply(x2,c(c(1:3)[-r][2],r),sum)
	
	# bxz <- classcrit(xz)
	# byz <- classcrit(yz)
	
	# nxz <- bxz + classcrit(xz, FALSE)
	# nyz <- byz + classcrit(yz, FALSE)
	
	# iz <- iccrit(apply(x2,r,sum))
	
	# v1 <- (nxz*nyz)/iz - (bxz*byz)/iz
	
	# #v1<- neg3t(x)*3*iccrit(x)/v1
	# v1 <- BCC(x)/v1
	# return(v1)
	
# }

