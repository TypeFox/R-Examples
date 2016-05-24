MRwarp <-
function(Xdata,Ydata,chain=400,thin=10,burnin=200,kernel.s,components=1,selection="FIXED",shr=0.3,
outputfit=1,alpha=0.1)
{
select=0
model.selection=0

S <- dim(Xdata)[1]
N <- dim(Xdata)[2]
X <- matrix(0,1,S*N)
Y <- matrix(0,1,S*N)
Q <-1
for (i in 1:S)
{
X[((i-1)*N+1):(i*N)] <- Xdata[i,]
Y[((i-1)*N+1):(i*N)] <- Ydata[i,]
}


# input error messages
#-----------------------
if (dim(Xdata)[1]==1){stop("input: X matrix must have at least two rows")}
if (dim(Ydata)[1]==1){stop("input: Y matrix must have at least two rows")}
if (dim(Ydata)[1]!=dim(Xdata)[1]){stop("input: Y matrix must have the same dimension as the X matrix")}
if (dim(Ydata)[2]!=dim(Xdata)[2]){stop("input: Y matrix must have the same dimension as the X matrix")}
if (burnin < 0 || burnin >= chain){stop("input: burnin must be a positieve integer, strictly smaller than the chainsize")}
if (as.integer(chain)!=chain || chain <=0){stop("input: chains must be a strictly positive integer")}
if (as.integer(thin)!=thin || thin <=0){stop("input: thin must be a strictly positive integer")}
if (shr < 0 ){stop("input: shrinkage parameter must be strictly positive")}



if (as.integer(length(kernel.s)/3) != length(kernel.s)/3){stop("input:  length of vector with starting values for the kernel parameters not a multiple of 3")}

if (selection == "FIXED")
	{
	if (as.integer(components) != components || components <= 0) {stop("input: Number of warping components must be a strictly positive integer")}
	}

# for storing the chains:
#-------------------------
Chain.warp <- matrix(0,1,(chain -burnin) * (3*Q + Q*(S-1)+1))
Chain.logsigma <- matrix(0,1,(3*Q + Q*(S-1)+1))
Chain.shift <- matrix(0,1,(chain -burnin) *(S-1))
Chain.amp <- matrix(0,1,(chain -burnin)* length(kernel.s)/3)
Chain.ker <- matrix(0,1,(chain -burnin) * length(kernel.s))
Namp <- length(kernel.s)/3
iternr <- 0

# The C function call
#---------------------------------------
stop=0

if (selection == "STEP"){selectionm = as.integer(c(1,1)) }else{
	selectionm = as.integer(c(2,components))	}

ll <- min(X) 
ul <- max(X)
lla <- min(Y) 
ula <- max(Y)

priors <- matrix(0,Q*2*(3 + (S-1))+1,1)

q = seq(-3,3,by=0.01)
pnorm <- pnorm(q, mean = 0, sd = 1)
ans.prev <- 0
ans <- 0

if (outputfit==1)
	{
	dev.new(width=5,height=5)
	ylim = c(lla,ula)
	xlim = c(ll,ul)	
	plot(Xdata[1,],Ydata[1,],ylim=ylim,xlim=xlim,main="Original data",xlab="",ylab="")
	lines(Xdata[1,],Ydata[1,])
	for (j in 2:S)
		{
		points(Xdata[j,],Ydata[j,],col=j)
		lines(Xdata[j,],Ydata[j,],col=j)
		}
	}

# priors first component

f <- 1:14
b <- 1:15
indexf <- Q + (f-1)*Q
indexb <- Q + (b-1)*Q
breaksv <- seq(from = ll-(ul-ll)/10, to = ul+(ul-ll)/10, length.out=15)
priors[indexf]<-(chain-burnin)/14
priors[14*(3 + (S-1))*Q+indexb]<-breaksv
priors[14*Q + indexf]<-(chain-burnin)/14
priors[14*(3 + (S-1))*Q+15*Q+indexb]<-breaksv
priors[14*Q*2 + indexf]<-(chain-burnin)/14
priors[14*(3 + (S-1))*Q+15*Q*2+indexb]<-breaksv	
for (j in 1:(S-1))
	{
	breaksv <- seq(from = -1, to = 1, length.out=15)
	priors[14*(Q*3+(Q-1)*(S-1)+(j-1))+indexf]<-(chain-burnin)/14
	priors[14*(Q*(3+S-1))+15*(Q*3+Q*(j-1))+indexb]<-breaksv
	}

while (stop==0)
	{
	cat("start C loop")
	ans <- .C("interfaceRC",as.double(X),as.double(Y),as.integer(S),as.integer(N), as.integer(chain), as.integer(thin),as.integer(burnin),as.integer(Namp),as.double(kernel.s), as.integer(Q),as.double(shr),as.double(priors),as.double(pnorm),as.double(Chain.warp),as.double(Chain.shift),as.double(Chain.amp),as.double(Chain.ker),as.double(Chain.logsigma),as.integer(iternr),PACKAGE="MRwarping")
	cat("C loop finished")	

	# organizing output
	kernels <- t(matrix(ans[[17]],length(kernel.s),chain-burnin))
	warping <- t(matrix(ans[[14]],Q*3+Q*(S-1)+1,chain-burnin))
	shift <- t(matrix(ans[[15]],(S-1),chain-burnin))
	ampl <- t(matrix(ans[[16]],length(kernel.s)/3,chain-burnin))
	logsigma <- t(matrix(ans[[18]],Q*3+Q*(S-1)+1,chain-burnin))
	iter <- ans[[19]]

	if (outputfit==1)
		{
		# shift
		WX <- Xdata
		for (j in 1:(S-1))
			{
			WX[j,] <- WX[j,]+shift[iter,j]
				}
		WX[S,] <- WX[S,] - sum(shift[iter,])

		#warping
	
		A <- warping[iter,1:Q]
		Ll <- warping[iter,(Q+1):(2*Q)]
		Ul <- warping[iter,(2*Q+1):(3*Q)]
		Lambda <- matrix(0,S,Q)
		for (j in 1:(S-1))
			{
			for (jj in 1:Q)
				{
				Lambda[j,jj] <-  warping[iter,(3*Q)+(jj-1)*(S-1)+j]
				}
			}
		if (Q==1){Lambda[S]=-sum(Lambda)}else{
			Lambda[S,]=-colSums(Lambda[1:(S-1),])}
			
		for (j in 1:S)
			{
			WX[j,] <- warp(A,Lambda[j,],A-Ll,Ul-A,WX[j,])
			}

		dev.new(width=5,height=5)
		ylim = c(lla,ula)
		xlim = c(ll,ul)	
		plot(WX[1,],Ydata[1,],ylim=ylim,xlim=xlim,main=paste("Shift and ",Q," component(s)"),xlab="",ylab="")
		lines(WX[1,],Ydata[1,])
		for (j in 2:(S))
			{
			points(WX[j,],Ydata[j,],col=j)
			lines(WX[j,],Ydata[j,],col=j)
			}
		}
	model.selection=0
	a.int <- boa.hpd(warping[,Q], alpha)
	ll.int <- boa.hpd(warping[,2*Q], alpha)
	ul.int <- boa.hpd(warping[,3*Q], alpha)
	L.int <- matrix(0,S,2)
	sum = 0
	for (j in 1:(S-1))
		{
		L.int[j,] <- boa.hpd(warping[,3*Q + (S-1)*(Q-1)+j], alpha)
		if (0 > L.int[j,1] && 0 < L.int[j,2]){sum=sum+1}
		}
	L.int[S,] <- boa.hpd(-rowSums(warping[,(3*Q + (S-1)*(Q-1)+1):(3*Q + (S-1)*(Q-1)+S-1)]), alpha)
	if (0 > L.int[S,1] && 0 < L.int[S,2]){sum=sum+1}
	if  (ll.int[2] > ul.int[1]) {model.selection=1}else{
		if  (sum == S) {model.selection=1}
		}
	if (model.selection==1)
		{
		cat(paste("Bayesian model selection criterion suggests not to include this component"))
		select=1
		}else{cat(paste("Bayesian model selection criterion includes this component"))}
	answer="n"
	if (selectionm[1] == 1)
		{ 
		stepbystep <- function(){readline("Do you want to continue and add a component? (y/n)")} 
		error<- function(){readline("Error reading input. Please answer with 'y' or 'n' only.")} 
		answer <- stepbystep() 
		while (answer != "y" && answer != "n") {answer <- error()}
		if (answer == "n")
			{
			stop = 1
			cat(paste("program stopped after ",Q," components"))
			}
		}

		if (selectionm[1]==2 && Q == components)
			{
			stop=1
			}
		if ( (selectionm[1]==2 && Q+1 <= components)|| answer == "y") #add new component in next model
			{
			Q <- Q+1

			Chain.warp <- matrix(0,1,(chain -burnin) * (3*Q + Q*(S-1)+1))
			Chain.logsigma <- matrix(0,1,(3*Q + Q*(S-1)+1))
			Chain.shift <- matrix(0,1,(chain -burnin) *(S-1))
			Chain.amp <- matrix(0,1,(chain -burnin)* length(kernel.s)/3)
			Chain.ker <- matrix(0,1,(chain -burnin) * length(kernel.s))

			#transfer starting values for next model
			Chain.warp[1:(3*(Q-1)+ (Q-1)*(S-1)+1)] <- warping[iter,] 
			Chain.shift[1:(S-1)] <- shift[iter,] 
			Chain.ker[1:(length(kernel.s))] <- kernels[iter,] 
			Chain.amp[1:(length(kernel.s)/3)] <- ampl[iter,] 
			Chain.logsigma[1:(3*(Q-1) + (Q-1)*(S-1)+1)] <- logsigma[iter,] 

			# new priors
			priors <- matrix(0,14*(3 + (S-1))*Q+15*(3 + (S-1))*Q)
			for (i in 1:(Q-1))
				{
				f <- 1:14
				b <- 1:15
				indexf <- i + (f-1)*Q
				indexb <- i + (b-1)*Q
				# a
				breaksv <- seq(from = min(warping[,i]), to = max(warping[,i]), length.out=15)
				priors[indexf]<-hist(warping[,i], breaks=breaksv,plot=FALSE)$counts
				priors[14*(3 + (S-1))*Q+indexb]<-breaksv
				
				# ll
				breaksv <- seq(from = min(warping[,Q-1+i]), to = max(warping[,Q-1+i]), length.out=15)
				priors[14*Q+indexf]<-hist(warping[,Q-1+i], breaks=breaksv,plot=FALSE)$counts
				priors[14*(3 + (S-1))*Q+15*Q+indexb]<-breaksv
				
				# ul
				breaksv <- seq(from = min(warping[,2*(Q-1)+i]), to = max(warping[,2*(Q-1)+i]), length.out=15)
				priors[14*Q*2+indexf]<-hist(warping[,2*(Q-1)+i], breaks=breaksv,plot=FALSE)$counts
				priors[14*(3 + (S-1))*Q+15*Q*2+indexb]<-breaksv
		
				# intensities
				for (j in 1:(S-1))
					{
					indexf <- (i-1)*(S-1) + (f-1)*Q*(S-1) +(j)
					indexb <- (i-1)*(S-1) + (b-1)*Q*(S-1) +(j)
					breaksv <- seq(from = min(warping[,3*(Q-1)+(S-1)*(i-1)+j]), to = max(warping[,3*(Q-1)+(S-1)*(i-1)+j]), length.out=15)
					priors[14*(Q*3)+indexf]<-hist(warping[,3*(Q-1)+(S-1)*(i-1)+j], breaks=breaksv,plot=FALSE)$counts
					priors[14*(Q*(3+S-1))+15*(Q*3)+indexb]<-breaksv
					}
				}
			# new component
			indexf <- Q + (f-1)*Q
			indexb <- Q + (b-1)*Q
			breaksv <- seq(from = ll-(ul-ll)/10, to = ul+(ul-ll)/10, length.out=15)
			priors[indexf]<-(chain-burnin)/14
			priors[14*(3 + (S-1))*Q+indexb]<-breaksv

			priors[14*Q + indexf]<-(chain-burnin)/14
			priors[14*(3 + (S-1))*Q+15*Q+indexb]<-breaksv

			priors[14*Q*2 + indexf]<-(chain-burnin)/14
			priors[14*(3 + (S-1))*Q+15*Q*2+indexb]<-breaksv	

			for (j in 1:(S-1))
					{
					indexf <- (Q-1)*(S-1) + (f-1)*Q*(S-1) +(j)
					indexb <- (Q-1)*(S-1) + (b-1)*Q*(S-1) +(j)
					breaksv <- seq(from = -1, to = 1, length.out=15)
					priors[14*(Q*3)+indexf]<-(chain-burnin)/14
					priors[14*(Q*(3+S-1))+15*(Q*3)+indexb]<-breaksv
					}

			# for final output
			kernels.prev <- kernels
			warping.prev <- warping
			shift.prev <- shift
			iter.prev <- iter
			}

	}

if (Q==1)
{
kernels.prev <- kernels
warping.prev <- warping
shift.prev <- shift
iter.prev <- iter
}

# unload .dll files
#----------------
#dyn.unload("interfaceRC_warp.dll")
#----------------

lambdas.output <- matrix(0,chain-burnin,S*Q)
for (i in 1:Q)
	{
	lambdas.output[,(S*(i-1)+1):(S*(i-1)+S-1)] <- warping[,(3*Q + (S-1)*(i-1)+1):(3*Q + (S-1)*(i-1)+S-1)]
	lambdas.output[,S*(i-1)+S] <- -rowSums(warping[,(3*Q + (S-1)*(Q-1)+1):(3*Q + (S-1)*(Q-1)+S-1)])
	}

warpingpars = list(lower=warping[,(Q+1):(2*Q)],A=warping[,1:Q],upper=warping[,(2*Q+1):(3*Q)],Intensities=lambdas.output)

last=list(shift= cbind(shift,-rowSums(shift[,])),warping= warpingpars,kernels=kernels,error.variance=warping[,3*Q+(S-1)*Q+1],max.post.dens=iter )
if (Q!=1)
{Q <- Q-1}
lambdas.output <- matrix(0,chain-burnin,S*Q)
for (i in 1:(Q))
	{
	lambdas.output[,(S*(i-1)+1):(S*(i-1)+S-1)] <- warping.prev[,(3*Q + (S-1)*(i-1)+1):(3*Q + (S-1)*(i-1)+S-1)]
	lambdas.output[,S*(i-1)+S] <- -rowSums(warping.prev[,(3*Q + (S-1)*(Q-1)+1):(3*Q + (S-1)*(Q-1)+S-1)])
	}

warpingpars = list(lower=warping.prev[,(Q+1):(2*Q)],A=warping.prev[,1:Q],upper=warping.prev[,(2*Q+1):(3*Q)],Intensities=lambdas.output)

previous=list(shift= cbind(shift.prev,-rowSums(shift.prev[,])),warping= warpingpars,kernels=kernels.prev,error.variance=warping.prev[,3*Q+(S-1)*Q+1],max.post.dens=iter.prev )



list(last=last,previous=previous)
}
