# gauss_dist <- function(vector, pos, sigma, amp)
# {
	# z.par <- (vector-pos)/sigma
	# gaussian.function <- amp*exp(-(z.par^2)/2)	
	# gaussian.function
# }

# normalize <- function(x)
# {
	# x[is.na(x)] <- 0
	# if(is.matrix(x)==T) norm.x <- sweep(x, 2, apply(x, 2, function(k) max(k, na.rm=T)), "/")
	# if(is.matrix(x)==F) norm.x <- x/max(x, na.rm=T)
	# norm.x[is.na(norm.x)] <- 0
	# norm.x
# }

is.even <- function(x) x %% 2 == 0

paste.sp <- function(x,y) {paste(x,y, sep=",")}

avoid.processing <- function(sampleRD)
{
	#if(min(sampleRD@avoid.processing.mz)<sampleRD@min.mz) stop(paste("Error: Minimum Mz adquired is higher than Avoid Processing Mz. Please change Avoid Processing Mz parameter in SetAlPar() function, whith at least from Mz number",sampleRD@min.mz))
	avoid.mz <- sampleRD@avoid.processing.mz - (sampleRD@min.mz - 1)
	UnderMZ <- which(sampleRD@avoid.processing.mz<sampleRD@min.mz)
	AboveMZ <- which(sampleRD@avoid.processing.mz>sampleRD@max.mz)
	out.mz <- unique(c(UnderMZ, AboveMZ))
	if(length(out.mz)!=0) avoid.mz <- avoid.mz[-out.mz]

	sampleRD@data[,avoid.mz] <- 0
	sampleRD
}


get.nrowcol <- function(vect.ind,ncl,nrw)
{
	mat.indexes <- apply(as.matrix(vect.ind),1, function(position){
		numcl <- as.integer(position/nrw) + 1
		numrw <- position - (numcl-1)*nrw 
		c(numcl,numrw)	
		})
	t(mat.indexes)
}

sparse.to.vector <- function(sparse.profile)
{	
	splitted.profile.list <- strsplit(as.character(sparse.profile),split=" ")[[1]]
	splitted.profile.matrix <- apply(as.matrix(splitted.profile.list),1,function(c.bin){strsplit(c.bin,split=",")[[1]]})
	
	output <- list(time=as.numeric(splitted.profile.matrix[1,]), int=as.numeric(splitted.profile.matrix[2,]))
	output
}

get.max.sign <- function(mat) {
	apply(mat,2,function(X){if(max(X)>abs(min(X))){return(1)}else{return(-1)}})	
}

lagFFT <- function(x,y, corr.coef=F)
{
	#lag value: the 'lag' positions that vector 'y' has to be displaced to be aligned with 'x'
	x[is.na(x)] <- 0; y[is.na(y)] <- 0
	X <- fft(x)
	Y <- fft(y)
	lag.1 <-  which.max(abs(fft(X*Conj(Y),inverse=T))) - 1
	lag.2 <-  which.max(abs(fft(Y*Conj(X),inverse=T))) - 1	
	if(abs(lag.1)>abs(lag.2)) lag <- -lag.2
	if(abs(lag.1)<abs(lag.2)) lag <- lag.1
	if(abs(lag.2)==abs(lag.1)) lag <- lag.1*(-1)
	#cat("L1:",lag.1, "L2:", lag.2, "\n")
	##Correlation coefficients for x,y
	if(corr.coef==T) abs(fft(X*Conj(Y),inverse=T))/(sqrt(sum(x^2)*sum(y^2))*length(x))
	lag
}

log.error <- function(w,x,y) {sum(log(abs(x-w*y)+1))}

break.vector <- function(x)
{
	x.vect <- x
	x.vect.o <- x.vect
	x.vect[normalize(x.vect)<0.01] <- 0
	x.vect[x.vect!=0] <- 1
	x.vect[x.vect==0] <- -1
	sqnc <- cumprod(x.vect)*x.vect
	 sqnc.groups <- c(1,diff(sqnc)!=0)
	 split.grouping <- cumsum(sqnc.groups)		
	split.indexes <- unique(split.grouping[which(duplicated(split.grouping)==T)])
	x.vect.s <- matrix(0, nrow=length(x.vect.o), ncol=length(split.indexes))
	x.vect.s <- apply(as.matrix(split.indexes),1,function(i) { x.out <- x.vect.o*0
		x.out[which(split.grouping==i)] <- x.vect.o[which(split.grouping==i)]
		x.out})
	x.vect.s
}
