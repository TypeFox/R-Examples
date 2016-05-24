logLikStrat <-
function(x11,x14,x24,y,beta){
	x11 <- as.matrix(x11)
	x14=as.matrix(x14)
	x24=as.matrix(x24)
    x24 <- as.matrix(x24)
	n <- dim(x11)[1]
	n11 <- dim(x11)[2]
	n14 <- dim(x14)[2]
	n24 <- dim(x24)[2]
	b11 <- as.matrix(beta[1:n11])
	#b13 <- as.matrix(beta[(n11+1)])
	b14 <- as.matrix(beta[(n11+1):(n11+1+n14-1)])
	b24 <- as.matrix(beta[-c(1:(n11+1+n14-1))])
	u24 <- x24%*%b24
	p2 <- pnorm((u24))
	u11 <- x11%*%b11
	#u13 <- b13*rep(1,n)
	u14 <- x14%*%b14
	u12 <- p2*u14 #+ (1-p2)*u13
	p1 <- pnorm((u12-u11))
	
	# double-check
	vec1 <- which(y==1); vec4 <- which(y==4); vec3 <- which(y==3) # which obs belong where?
	
	
	#vec1 <- which(y[,1]==0)
	#vec3 <- which(y[,2]==0)
	#vec4 <- which(y[,2]==1)
	chunk1 <- log(1-p1)[vec1]
	chunk3 <- log(p1*(1-p2))[vec3]
	chunk4 <- log(p1*p2)[vec4]
	ll <- sum(chunk1) + sum(chunk3) + sum(chunk4) 
	return(-ll)	
	}
