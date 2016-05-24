chrom.isoreg <-
function(A, tresh=0.2)
{
	A <- as.matrix(A)
	if(all(is.na(A))) return(A)
	if(ncol(A)==0) return(A)
	#max.ind <- apply(abs(A),2,which.max)
	#max.sign <- sweep(A,2,max.ind)
	#mul.sign <- diag(as.matrix(apply(as.matrix(A[max.ind,]),2,sign)))
	#A <- sweep(as.matrix(A),2,mul.sign,"*")
	#x <- A[,2]
		
	Ar <- apply(as.matrix(A),2,function(x)
	{
		x[x<0] <- 0
		max.x <- max(x, na.rm=T)
		x <- normalize(x)
		x.rigth <- x[which.max(x):1]
		if(which.max(x)<length(x)) {x.left <- x[(which.max(x)+1):length(x)]}else{x.left <- NA}
		
		if(length(x.rigth)>2)
		{
			flag.isoreg <- F
			for(i in 2:length(x.rigth))
			{
				if(x.rigth[i]<tresh) flag.isoreg=T
				if(flag.isoreg) 
					if(x.rigth[i]>x.rigth[i-1]) x.rigth[i] <- x.rigth[i-1]
			}
		}
		if(length(x.left)>2)
		{
			flag.isoreg <- F
			for(i in 2:length(x.left))
			{
				if(x.left[i]<tresh) flag.isoreg=T
				if(flag.isoreg) 
					if(x.left[i]>x.left[i-1]) x.left[i] <- x.left[i-1]
			}
		}
		x <- c(x.rigth[length(x.rigth):1],x.left)
		del.na <- which(is.na(x)==T)
		if(length(del.na)!=0) x <- x[-del.na]
		#length(x)
		#plot(x, type="l")
		x <- x*max.x
		x
	})
	
	xmax.vect <- apply(Ar,2,max)
	Ar.diff <- normalize(apply(Ar,2,diff))
	Ar.n <- apply(Ar.diff,2, function(x){	
		x[which(abs(x[1:which.max(x)])<0.01)] <- 0
		x[(which.min(x)-1)+ which(abs(x[which.min(x):length(x)])<0.01)] <- 0
		x <- cumsum(c(0,x))
		x[x<0] <- 0
		x[normalize(x)<0.01] <- 0
		x
	})
	
	Ar.n <- sweep(normalize(Ar.n),2,xmax.vect,"*")
	Ar.n	
}
