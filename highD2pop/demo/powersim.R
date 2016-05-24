
##### Set simulation parameters here. 

simparms <- list	(	
					n = c(30,45),
					m = c(45,60),
					p = 100,
					beta = c(0,.1,.2,.4,.6,.8,.9,1),
					delta = .5,
					muX = function(p) return(rep(0,p)),
					muY = function(p,b,delta) { 
												times1 <- b * p
												times2 <- (1-b)*p+1e-7
												mu <- c(rep(delta,times1),rep(0,times2))
												return(mu)
											},
					commoncov = TRUE,
					VarScaleY = 2,
					dep = c('IND','ARMA','LR'),
					ARMAparms = list(coefs=list(ma=c(.2,.3) , ar=c(.4,-.1))),
					LRparm = .75,	
					S = 25,
					#innov = function(n,...) rdblepareto(n,16.5,8),
					#innov = function(n,...) rnorm(n,0,1),
					innov = function(n,...) rgammashift(n,4,2),
					heteroscedastic=FALSE,
					het.diag = function(p) return(diag(1/2 + rexp(p,2)))
					)

					
testparms <- list(	r = 15,
			smoother = c("parzen","trapezoid"),
			ntoorderminus=c(0,1,2))
			
							
library(highD2pop)


highD2popBetaPower.sim <- function 	(	
										ChenQin = FALSE,
										SK	= FALSE,
										CLX = FALSE,
										GCT = TRUE,
										simparms,
										testparms
										)
{
	# Hmmm, I may just have to redefine the function if I change my simparms...
	
	output<-list
	output<-list(output,output,output,output,output,output,output,output)  # l
	output<-list(output,output,output)	# k
	output<-list(output)	# j
	output<-list(output,output) 	# i
	
	for(i in 1:length(simparms$n)){
	for(j in 1:length(simparms$p)){
	for(k in 1:length(simparms$dep)){
	for(l in 1:length(simparms$beta)){
		
		i<-1
		j<-1
		k<-1
		l<-3

		DATA <-	build2popData(	
								n = simparms$n[i],
								m = simparms$m[i],
								p = simparms$p,
								muX = simparms$muX(simparms$p),
								muY = simparms$muY(simparms$p,simparms$beta[l],simparms$delta),
								commoncov = simparms$commoncov,
								VarScaleY = simparms$VarScaleY,
								dep = simparms$dep[k],
								ARMAparms = simparms$ARMAparms,
								LRparm = simparms$LRparm,
								S = simparms$S,
								innov = simparms$innov,
								heteroscedastic=simparms$heteroscedastic,
								het.diag = simparms$het.diag(simparms$p)
		  					)
		 						
		  						
		if(SK==TRUE) SK.results <- SK.sim(DATA) else SK.results=NULL
		if(ChenQin==TRUE) ChenQin.results <- ChenQin.sim(DATA) else ChenQin.results=NULL
		if(CLX==TRUE) CLX.results <- CLX.sim.Covtest(DATA) else CLX.results=NULL
		if(GCT==TRUE) GCT.results11 <- GCT.sim(DATA,r = testparms$r[1],smoother=testparms$smoother[1],ntoorderminus=testparms$ntoorderminus[1]) else GCT.results1=NULL
		if(GCT==TRUE) GCT.results12 <- GCT.sim(DATA,r = testparms$r[1],smoother=testparms$smoother[1],ntoorderminus=testparms$ntoorderminus[2]) else GCT.results2=NULL
		if(GCT==TRUE) GCT.results13 <- GCT.sim(DATA,r = testparms$r[1],smoother=testparms$smoother[1],ntoorderminus=testparms$ntoorderminus[3]) else GCT.results3=NULL
		if(GCT==TRUE) GCT.results21 <- GCT.sim(DATA,r = testparms$r[1],smoother=testparms$smoother[2],ntoorderminus=testparms$ntoorderminus[1]) else GCT.results1=NULL
		if(GCT==TRUE) GCT.results22 <- GCT.sim(DATA,r = testparms$r[1],smoother=testparms$smoother[2],ntoorderminus=testparms$ntoorderminus[2]) else GCT.results2=NULL
		if(GCT==TRUE) GCT.results23 <- GCT.sim(DATA,r = testparms$r[1],smoother=testparms$smoother[2],ntoorderminus=testparms$ntoorderminus[3]) else GCT.results3=NULL

		
		output[[i]][[j]][[k]][[l]] <- list(	SK.results = SK.results, 
										ChenQin.results = ChenQin.results, 
										CLX.results = CLX.results,
										GCT.results11 = GCT.results11,
										GCT.results12 = GCT.results12,
										GCT.results13 = GCT.results13,										
										GCT.results21 = GCT.results21,
										GCT.results22 = GCT.results22,
										GCT.results23 = GCT.results23,
										n = simparms$n[i],
										m = simparms$m[i],
										p = simparms$p,
										muX = simparms$muX(simparms$p),
										muY = simparms$muY(simparms$p,simparms$beta[l],simparms$delta),
										commoncov = simparms$commoncov,
										VarScaleY = simparms$VarScaleY,
										dep = simparms$dep[k],
										ARMAparms = simparms$ARMAparms,
										LRparm = simparms$LRparm,
										S = simparms$S,
										innov = simparms$innov,
										heteroscedastic=simparms$heteroscedastic,
										het.diag = simparms$het.diag(simparms$p)
										)										
			print(paste(i,j,k,l))
										
	}						
	}
	}
	}
	
	return(output)
	
}

output <- highD2popBetaPower.sim(SK=TRUE,ChenQin=TRUE,CLX=TRUE,GCT=TRUE,simparms=simparms,testparms=testparms)


##### Plot the results.


par(mfrow = c(2,3))

par(oma = c(2.5,2.5,4,1),
	mar = c(0,0,0,0),
	lwd = 1,
	cex.axis = .8,
	cex.lab = .85,
	tck = -.02)
	
Pzen0 <- Pzen2 <- ChQ <- SK <- CLX <- matrix(0,length(output),length(output[[1]][[1]][[1]]))


for( i in 1:length(output))
	for(j in 1:3)
	{
		for(l in 1:length(output[[1]][[1]][[1]]))
			{
				Pzen0[i,l] <- mean(output[[i]][[1]][[j]][[l]]$GCT.results11[,2] < .05) # parzen, r = 20
				Pzen2[i,l] <- mean(output[[i]][[1]][[j]][[l]]$GCT.results13[,2] < .05) # parzen, r = 20
				ChQ[i,l] <- mean(output[[i]][[1]][[j]][[l]]$ChenQin.results[,2] < .05) # Ch-Q
				SK[i,l] <- mean(output[[i]][[1]][[j]][[l]]$SK.results[,2] < .05) # SK
				CLX[i,l] <- mean(output[[i]][[1]][[j]][[l]]$CLX.results[,2] < .05) # CLX
			}

			plot(	Pzen0[i,],
					ylim = c(0,1),
					xlim = c(1,8),
					yaxt = "n",
					ylab = "",
					type = "l",
					xaxt = "n",
					xlab = "",
					lwd = 1)	

			lines(Pzen0[i,],lwd=1)
			lines(Pzen2[i,],lty=5,lwd=1)
			lines(ChQ[i,],lty=3,lwd=1)
			lines(SK[i,],lty=2,lwd=1)
			lines(CLX[i,],lty=4,lwd=1)
		
				
	if(i==1 & j==2){
					
		x.pos <- grconvertX(.5,from="nfc",to="user")
		y.pos <- grconvertY(1.125,from="nfc",to="user")
									
		legend(	x=x.pos,
				y=y.pos,
				lwd=c(1,1,1,1,1),
				lty=c(1,5,3,2,4),
				legend=c("mod-p GCT","lg-p GCT","Ch-Q","SK","CLX"),
				bty="n",
				seg.len=3,
				horiz = TRUE,
				xpd=NA,
				xjust = .5,
				cex=1.25)
	}
				
	if(i==1)
	{
		
		mtext(side=3,text=output[[1]][[1]][[j]][[1]]$dep,line=-1.5,cex=.85,adj=.1)
		text(x=1,y=.5,"(n,m) = (45,60)",srt=90,cex=.85)
		
	}
	
	if(i==2)
	{	
		axis(1,at=1:8, as.character(simparms$beta),padj=-1.5)
		text(x=1,y=.5,"(n,m) = (90,120)",srt=90,cex=.85)
	}
	if(j==1) 	axis(2,padj=1)
	
	abline(h=.05,lwd=.5)
	text(x=7.5,y=.025, expression(alpha==.05),cex=.6)
	
}						
								
	text(x=7.5,y=.025, expression(alpha==.05),cex=.6)
	
	mtext(outer=TRUE,side=1,expression(beta), line=1.5,cex=1)
	mtext(outer=TRUE,side=2,"Power",line=1.25,cex=1)
	mtext(outer=TRUE,side=3,text="Power curves under skewed innovations and equal covariance matrices",line=2.5,cex=1)


