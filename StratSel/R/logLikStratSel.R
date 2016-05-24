##
## "logLikStratSel" is function which computes the log-likelihood of an QRE agent error model which
## takes selection into account (Leemann, 2013):
##
##       1
##      /\
##     /  \
##    /    \ 2
##   u11   /\
##        /  \
##       /    \
##     u13    u14
##     0      u24
##
## whereas there are three possible outcomes (1,3, and 4). Note, to achieve identification u23 is
## set equal to zero (u23=0). Finally, the model implementation here differs slightley with respect
## to the assumed error variance. While Signorino and Kenkel and Signorino rely on \sqrt(2) I just
## use 1 which does not affect any substantive result (predicted probabilities, z-values, etc.).
## Unlike the Games package (Kenkel and Signorino, 2013) this returns the sum of the logLik and not
## the individual contributions.
##
## This model estimates the correlation of the random components of player 1 and player 2




logLikStratSel <- function(x11,x14,x24,y,beta){
	# NOTE: parametrization is now intuitive (first 1, then 2) unlike Greene/Wawro
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
	b24 <- as.matrix(beta[(n11+1+n14):(n11+1+n14+n24-1)])
	para <- as.matrix(beta[-c(1:(n11+1+n14+n24-1))])
	rho <- 2*(1/(1+exp(-para)))-1
	xb2.sp <- x24%*%b24												# latent D.utility for player 2 for a_4
	p.p <- pnorm(xb2.sp/(sqrt(1)))									# p(Y4=1|Y1=0)
	u1.34 <- -x11%*%b11 + p.p*(x14%*%b14) #+ (1-p.p)*rep(b13,n)		# latent D.utility for player 1 for a_2
	part1 <- cbind(-xb2.sp,u1.34); part2 <- cbind(xb2.sp,u1.34)		# combining for mvnorm
	covar1 <- matrix(c(1,-rho,-rho,1),2,2)				# negative rho
	covar2 <- matrix(c(1,rho,rho,1),2,2)				# positive rho
	prob1 <- pbivnorm(x=cbind(-xb2.sp, u1.34), rho=-rho)
	#prob2 <- apply(part2,1,pmnorm, mean=rep(0,2),varcov=covar2)	# y_1=1; y_2=1
	prob2 <- pbivnorm(x=cbind(xb2.sp, u1.34), rho=rho)
	prob3 <- 1-pnorm(u1.34)										# y_1=0; y_2=?
	#Y1 <- y[,1]; Y2 <- y[,2]
	Y <- y
	#vec3 <- which(Y1==0); vec2 <- which(Y2==1); vec1 <- which(Y2==0) # which obs belong where?
	vec3 <- which(Y==1); vec2 <- which(Y==4); vec1 <- which(Y==3) # which obs belong where?
	chunk1 <- sum(log(prob1[vec1]))								# for each y group, the
	chunk2 <- sum(log(prob2[vec2]))								# partial LL contrinution
	chunk3 <- sum(log(prob3[vec3]))								# it makes
	ll <- chunk1 + chunk2 + chunk3
	return(-ll)
	}
