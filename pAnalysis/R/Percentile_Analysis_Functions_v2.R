##############################################################################################################################
#
# The R-squared calculation based on two numeric vectors of equal length
#
#
R2 <- function(x, y){
	if(length(x)!=length(y)){stop("r.sq measure: length of x must equal length of y")}
	xh <- x-mean(x)
	yh <- y-mean(y)
	
	num <- sum(xh*yh)^2
	den <- sum(xh^2)*sum(yh^2)
	
	R2 <- num/den
	return(R2)
}


##############################################################################################################################
#
# simple function to capitalize first letters of words for use in titles
#
#
cap1 <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s,1,1)), substring(s, 2),
          sep="", collapse=" ")
}


##############################################################################################################################
#
# generate a data frame with possible values of R2 going from 0 to 1 with corresponding 
#   probabaility density and cumulative density functions given a particular number of degrees
#	of freedom (dof) and a particular noise distribution function (dist).
#	By default, use a million samples (order=6), normal distribution with mean=0 and sd=1, and 
#   return values with 3 decimal places.
#
#
pcdfs <- function(dof, order=6, ndecimals=3, dist='normal', par1=0, par2=1){

if(order>=7){stop(paste("order is too large (10M) -- calculation time too long. Make order<7. Fractions OK."))}
N=round(10^order,0)
bw=1/10^ndecimals #bin width

dN = dof*N
rnums <- rep(0,dN)
R2c <- rep(0,N)

if(				dist=="normal")		{	rnums <- rnorm(	n=dN, mean=par1, sd=par2)
	} else if(	dist=="uniform")	{	rnums <- runif(	n=dN, min=par1, max=par2)
	} else if(	dist=="lognormal")	{	rnums <- rlnorm(n=dN, meanlog=par1, sdlog=par2)
	} else if(	dist=="poisson")	{	rnums <- rpois(	n=dN, lambda=par1)
	} else if(	dist=="binomial")	{	rnums <- rbinom(n=dN, size=par1, prob=par2)}

#define R2 components x and y (y is just e, residuals, random numbers)				
x  <- seq(1,dof)
xb <- mean(x)
xd <- x-xb
xd <- t(matrix(rep(xd, N), nrow=dof, ncol=N))

e  <- matrix(rnums, nrow=N, ncol=dof)	#make a matrix of N rows by n columns
eb <- rowSums(e)/dof					#get means for each row
ed <- e-eb								#get delta

#calculate R2 numerator
n1 = xd*ed
n1s = rowSums(n1)
num = n1s*n1s							#numerator

#calculate R2 denominator
d1 = xd^2
d2 = ed^2
d1s = rowSums(d1)
d2s = rowSums(d2)
den = d1s*d2s							#denominator

#calculate R2 (this is an array of R2 calculations based on noise)
R2c = num/den							#R2

#separate R2s into bins of width bw (histogram R2) to get pdf.  sum pdf to get cdf.
br = seq(0, 1, by=bw)
R2 = br[2:length(br)]					#R2 values (bins)
R2h = hist(R2c, breaks=br, plot=F)		#histogrammed R2
pdf <- R2h$counts/sum(R2h$counts)		#prob density
cdf <- cumsum(pdf)						#cumulative probability density
R2df <- data.frame(R2, pdf, cdf)		#

return(R2df)
}


##############################################################################################################################
#
#	determine the baseline noise level (R2p) for a corresponding number of degrees of freedom(dof) and noise percentile(pct)
#
#
R2p <- function(dof, pct=0.95, ndecimals=3,...){
	cdf <- pcdfs(dof, ndecimals=ndecimals,...)[,c(1,3)]		#just R2 and cdf columns
	R2p <- cdf$R2[cdf[,2]>=pct][1]		#R2p is the value of cdf just at the point where it's just >= p
	R2p <- R2p + rnorm(1)*10^(-(ndecimals+2))  #add a small random number to remove any binning errors.
	R2p <- round(R2p,ndecimals)
	return(R2p)
}


##############################################################################################################################
#
#	For a given value of R2, dof and pct, determine the noise-normalized, dof-independent, 
#      distribution-independent, R2 equivalent:  R2k
#
R2k <- function(R2, dof, pct=0.95, ndecimals=3,...){
	r2p <- R2p(dof=dof, pct=pct, ndecimals=ndecimals,...)								#find R2p
	
	r2k <- (R2-r2p)/(1-r2p + 0.00000000001)			#rescale R2 to where it lies between 1 and baseline noise; add smidge incase r2p=1
	
	#make R2k consistent with the number of decimal places in R2p.
	fl <- floor(r2p)								
	nd=rep(0,length(R2))
	if(r2p-fl>0){
		nd <- nchar(sapply(strsplit(as.character(r2p), ".",fixed=T), "[[", 2))} 
	r2k <- round(r2k, ndecimals)
	
	return(r2k)
}



##############################################################################################################################
#
# Generate a percentile analysis table listing baseline R2ps for various degrees of freedom and percentiles. 
#		Any measured R2 falling below these values (for the corresponding dof and pct)
#			1) are indistinguishable from noise and
#			2) will yield a negative R2k
#			3) should be discarded
#
#		This will run a call to R2p for each combination of dof and pct
#		A dof list of one dof (say,10), will take one pct just under a minute, each additional pct adds an equivalent amount (approx).
#		60 dof will take one pct about 5.75 min.
#
#		See plotR2Equiv for a plot of all R2 equal to a particular measure of R2 (with a certain dof and pct).
#
R2pTable <- function(doflist=NULL, pctlist=NULL, order=4, ndecimals=2,...){
	if(is.null(doflist)){doflist=c(4,8,16,32,64,128)}
	if(is.null(pctlist)){pctlist=c(0.7,0.9,0.95,0.99)}
	nds <- length(doflist)	#need test here for pos integer
	nps <- length(pctlist)	#need test here for nums >0 and <1
	rownams = as.character(doflist)
	colnams = as.character(pctlist)
	
	shell <- matrix(nrow=nds, ncol=nps)
	r2ptab <- matrix(mapply(function(x,i,j) R2p(doflist[i],pctlist[j], order=order, ndecimals=ndecimals,...), shell,row(shell),col(shell)), nrow=nds, ncol=nps)

	r2ptab <- as.data.frame(r2ptab)
	colnames(r2ptab) <- colnams
	rownames(r2ptab) <- rownams
	return(r2ptab)
}



##############################################################################################################################
#
# Plot the probability density function for a given number of degrees of freedom and noise distribution function
#
#
plotpdf <- function(dof, order=4, dist='normal',...){

df <- pcdfs(dof=dof,order=order,dist=dist,...)
N = 10^order
dist2 <- sapply(dist, cap1)
mxy = max(df$pdf)

plot <- ggplot(df) + 
		geom_point(aes(R2, pdf),size=1) +
		ggtitle(paste("Probability Density Function")) +
		ylim(0,mxy) +
		xlab(expression(R^2)) + 
		ylab("Probability Density") +
		ggtitle(paste("Probability Density Function")) +
		geom_text(aes(x=0.95,y=0.9*mxy,label=paste("Noise Distribution:",dist2,
													"\nDegrees of Freedom:",dof,
													"\nNumber of  Samples:",floor(N))),size=3,hjust=1)


return(plot)
}


##############################################################################################################################
#
# Plot the cumulative probability density function (cdf) for a given number of degrees of freedom and noise distribution function
#
#
plotcdf <- function(dof, order=4, dist='normal',...){  		#need to explicitly state distribiution here in order to get it into the plot title

r2cdf <- pcdfs(dof=dof,order=order,dist=dist,...)
cdf <- NULL													#see http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when.  Need this to eliminate a note during R CMD check
N = 10^order
dist2 <- sapply(dist, cap1)
mxy <- max(r2cdf$cdf)
plot <- ggplot(r2cdf) + 
		geom_point(aes(R2, cdf),size=1) +
		ylim(0,mxy) + 
		xlab(expression(R^2)) + 
		ylab("Cumulative Probability") +
		ggtitle(paste("Cumulative Probability Density Function")) +
		geom_text(aes(x=0.95,y=0.3*mxy,label=paste("Noise Distribution:",dist2,
													"\nDegrees of Freedom:",dof,
													"\nNumber of  Samples:",floor(N))),size=3,hjust=1)

return(plot)
}


##############################################################################################################################
#
# Plot R2p for a list of percentiles (not greater than 5 in the list)
#
#
plotR2p <- function(doflist=c(2:30), pctlist=c(0.95), order=4, ndecimals=3, ...){
	if(length(pctlist)>5){stop(paste("Too many percentiles to calculate", length(pctlist)))}

	doflist <- doflist[doflist>1] 
	doflim <- min(30, length(doflist))
	doflist <- doflist[1:doflim]
	
	pctlim <- min(5,length(pctlist))
	pctlist <- pctlist[1:pctlim]
	pctlist <- formatC(as.numeric(pctlist),width=(ndecimals+1),format='f',digits=ndecimals,flag='0')
	doflength <- length(doflist)
	pctlength <- length(pctlist)
	
	mcolor <- c("black", "blue", "red", "green", "darkgreen")
	sizes <- c(3.6, 3.2, 2.8, 2.4, 2.0)/2  #make the points of the first plots larger so they can be seen
	#sizes <- c(1, 0.8, 0.6, 0.4, 0.2)   #use these for geom_path instead of geom_point
	r2pdf <- R2pTable(doflist=doflist, pctlist=pctlist, order=order,...)
	mxy <- 0.9*max(r2pdf[,1])
	N = 10^order
	plt <- ggplot(r2pdf)

	if(pctlength>=1){plt <- plt + 
		geom_point(aes(as.numeric(row.names(r2pdf)), r2pdf[,1]), color=mcolor[1], size=sizes[1]) +
		geom_text(aes(x=max(doflist), y=mxy-0.00, label=paste0("p = ",pctlist[1])), color=mcolor[1], hjust=1, size=4)}
	if(pctlength>=2){plt <- plt + 
		geom_point(aes(as.numeric(row.names(r2pdf)), r2pdf[,2]), color=mcolor[2], size=sizes[2]) +
		geom_text(aes(x=max(doflist), y=mxy-0.05, label=paste0("p = ",pctlist[2])), color=mcolor[2], hjust=1, size=4)}
	if(pctlength>=3){plt <- plt + 
		geom_point(aes(as.numeric(row.names(r2pdf)), r2pdf[,3]), color=mcolor[3], size=sizes[3]) +
		geom_text(aes(x=max(doflist), y=mxy-0.10, label=paste0("p = ",pctlist[3])), color=mcolor[3], hjust=1, size=4)}
	if(pctlength>=4){plt <- plt + 
		geom_point(aes(as.numeric(row.names(r2pdf)), r2pdf[,4]), color=mcolor[4], size=sizes[4]) +
		geom_text(aes(x=max(doflist), y=mxy-0.15, label=paste0("p = ",pctlist[4])), color=mcolor[4], hjust=1, size=4)} 
	if(pctlength>=5){plt <- plt + 
		geom_point(aes(as.numeric(row.names(r2pdf)), r2pdf[,5]), color=mcolor[5], size=sizes[5]) +
		geom_text(aes(x=max(doflist), y=mxy-0.20, label=paste0("p = ",pctlist[5])), color=mcolor[5], hjust=1, size=4)}
	
	plt <- plt + 
	ggtitle("R2 Baseline Noise Level (R2p) \nfor Various Noise Percentiles (p)") +
	xlab("Degrees of Freedom") +
	ylab(expression(R^2)) +
	geom_text(aes(x=max(doflist), y=mxy-0.25, label=paste0("Number of Samples:",N)), color='black', hjust=1, size=3)

	
	return(plt)
}




##############################################################################################################################
#
# Plot R2k for a single R2 across a range of dofs
#
plotR2k <- function(R2, doflist=c(2:30), pct=0.95, order=4, ndecimals=3,...){

	pct <- pct[1]										#ensure only one pct is used
	df <- R2pTable(doflist=doflist, pctlist=pct, ndecimals=ndecimals, order=order,...)

	df$R2k <- NA

	n <- nrow(df)

	for(i in 1:n){df$R2k[i] <- R2k(R2, dof=as.numeric(row.names(df)[i]), pct=pct, ndecimals=ndecimals, order=order,...)}
	
	maxx=max(doflist)
	plot <- ggplot(df) + 
		geom_point(aes(as.numeric(row.names(df)),df[,1]),color='red') + 		#column df[,1] is named for the pct used, which can change every time.
		geom_point(aes(as.numeric(row.names(df)),R2k),color='blue',na.rm=T) + 
		geom_hline(aes(yintercept=R2),color='black') + 
		ylim(0,1) +
		xlab("Degrees of Freedom") +
		ylab("R2") +
		ggtitle("R2k for a Given Baseline Noise Level (R2p) and \na Constant Measured R2") +
		geom_text(aes(x=maxx, y=0.17),		label=paste0("Baseline Noise Level\np = ",pct), 	color='red',  hjust=1) + 
		geom_text(aes(x=maxx, y=(R2-0.1)),	label=paste0("R2k"), 					color='blue', hjust=1) +
		geom_text(aes(x=maxx, y=(R2+0.05)),	label=paste0("Measured R2 = ",R2), 		color='black',hjust=1)
	
	return(plot)
	
}

##############################################################################################################################
#
# Plot the R2s for a range of dofs that are equivalent to a single measured R2
#		Plot the measured value of R2 (green asterisk)
#		Plot the noise level (R2p) (red)
#		Shade the area shoing improved R2 measures
#		Plot the R2k equivalent curve (black)
#		if desired, plot the noise level that equals R2
#
#
plotR2Equiv <- function(R2, dof, pct=0.95, order=4, plot_pctr2=F,...){

	mcolor <- c("red", "blue", "forestgreen", "slategray4", "gray20", "black")


	# get the pcdf for this dof
	df <- pcdfs(dof=dof, order=order,...)
	
	doflist = c(2:30)
	pctlist = c(pct)
	
	
	# this will add a plot of points that follow the curve where the pct equals R2 -- just to show the probability of noise for this R2.
	#plot_pctr2 = F
	# using that df, find the closest df$R2 (aka pct) to the given R2, or pct_r2
	if(plot_pctr2){
		pct_R2 <- df$cdf[df$R2>=R2][1]
		pctlist=c(pctlist,pct_R2)			#take out if not plotting pct_R2.  Comprising pct and pct_R2
		}

	doflength = length(doflist)
	pctlength = length(pctlist)

	ptable <- R2pTable(doflist=doflist,pctlist=pctlist,order=order,...)
	
	r2p <- ptable[(dof-1),1]
	r2k <- R2k(R2,dof=dof,...)
	f = (R2-r2p)/(1-r2p)
	ptable$R2Equiv <- f*(1-ptable[,1]) + ptable[,1]   	#this should be column 3 if pctr2 is being plotted and 2 if not

	tx = max(doflist[doflength])

	if(length(ptable)==3){ptable <- ptable[c(1,3,2)]}	#ensure proper order of columns;  if pct_r2=T, then swap 2<->3 to keep R2Equiv in pos 2

	plt <- ggplot(ptable) +
			geom_point(aes(as.numeric(row.names(ptable)),ptable[,1]),color=mcolor[1],size=2,na.rm=T) +
			geom_point(aes(as.numeric(row.names(ptable)),ptable[,2]),color=mcolor[6],size=2,na.rm=T) +
			geom_ribbon(aes(x=as.numeric(row.names(ptable)), ymin=ptable[,2], ymax=1),fill=mcolor[4],alpha=0.3,na.rm=T) +
			geom_point(data=data.frame(R2,dof), aes(dof,R2),shape=8, color=mcolor[3],size=5,na.rm=T) + 
			ggtitle("Noise Baseline and R2 Equivalent (R2k)") +
			xlab("Degrees of Freedom") +
			ylab(expression(R^2)) + 
			geom_text(x=tx, y=0.80, label=paste0("R2 = ",R2,"   dof = ",dof), color=mcolor[3], hjust=1, size=4) +
			geom_text(x=tx, y=0.75, label=paste0("R2k = ",r2k), color=mcolor[6], hjust=1, size=4) +
			geom_text(x=tx, y=0.70, label=paste0("Noise Baseline: R2p = ",pctlist[1]), color=mcolor[1], hjust=1,size=4) + 
			geom_text(x=tx, y=0.65, label=paste0("Improved R2 Measure"), color=mcolor[4], hjust=1, size=4)

	#if pct_r2 is T, plot the noise level where pct=R2
	if(plot_pctr2){plt <- plt + geom_point(aes(as.numeric(row.names(ptable)),ptable[,3]),color=mcolor[2],size=2) +
			geom_text(x=tx, y=0.60, label=paste0("R2 Noise Percentile = ",pctlist[2]), color=mcolor[2], hjust=1,size=4,na.rm=T)}
				
	#if R2 is in the noise, show the improved but still noisy R2 values in a black ribbon.
	if(any(ptable[,2]<ptable[,1])){ plt <- plt +
			geom_ribbon(aes(x=as.numeric(row.names(ptable)), ymin=ptable[,2], ymax=ptable[,1]),fill=mcolor[5],alpha=0.7,na.rm=T) +
			geom_text(x=tx, y=0.55, label=paste0("Unacceptable Noise"), color=mcolor[5], hjust=1, size=4) +
			geom_point(data=data.frame(R2,dof), aes(dof,R2),size=4,shape=8, color=mcolor[3],na.rm=T) }
				
			
	return(plt)
}

