#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#	TITLE: 	R - Meta-Analysis	- Features relating to further evidence			#
#	DATE:		23/02/2010										      #
#	R script for new graphical framework									#
#___________________________________________________________________________________________#

############################### Description of possible arguments ##########################################################
#SS - effect estimates of the current studies													#
#seSS - standard errors of the current studies													#
#sig.level - significance level															#
#method - "fixed" or "random"																#
#ylim - c(y1,y2) limits of the y axis														#
#xlim - c(x1,x2) limits of the x axis														#
#contour.points - number of points for additional features - more means a smoother contour but takes longer to compute	#
#isquared - e.g. c(0.4,0.5,0.6,0.7) values defining i-squared contours	- type element>1 for current isquared value	#
#tausquared - e.g. c(1,2,4,6) values defining tau-squared contours										#
#contour - TRUE/FALSE for displaying the significance contour											#
#legend - TRUE/FALSE for displaying key/legend													#
#expxticks - custom ticks for the x axis on a exponential scale										#
#xticks - custom ticks for the x axis														#
#yticks - custom ticks for the y axis														#
#zero - value for the null effect (usually 0, even when expxticks is used for odds ratios)						#
#xlab - label for the x axis																#
#ylab - label for the y axis																#
#plot.zero - TRUE/FALSE plot the null effect vertical line as defined by the arguament 'zero'					#
#plot.summ - TRUE/FALSE plot the pooled effect vertical line of current studies								#
#legendpos - position of the upper left hand corner on the x/y axis scale									#
#summ - summary diamond including pooled effect and confidence interval (significance level as defined by sig.level		#
#xpoints/ypoints - add an extra point(s) in the plot just to show as an example								#
#points - whether the study points should be displayed at all (TRUE default)								#
#pred.interval - displays predictive interval allong with the summary diamond.								#
#rand.load - show percentage of computations complete when the random effects contours are calculated				#
############################################################################################################################


##### START OF EXTFUNNEL FUNCTION
extfunnel <- function(SS, seSS, method, 
	sig.level=0.05, contour=FALSE, isquared=NULL, tausquared=NULL, contour.points=200, summ=FALSE,
	summ.pos=0, pred.interval=FALSE, plot.zero=FALSE, plot.summ=FALSE, ylim=NULL, xlim=NULL, legend=TRUE,
	expxticks=NULL, xticks=NULL, yticks=NULL, zero=0, xlab=NULL, ylab=NULL, rand.load=10,
	legendpos=c(xlim[2]+0.05*(xlim[2]-xlim[1]),ylim[2]), xpoints=NULL, ypoints=NULL, points=TRUE) {


############################### Description of possible arguments
#check for correct input later
if (!is.null(isquared) & !is.null(tausquared)) stop("Tau^2 and I^2 contours cannot be plotted on the same graph")

##### Checks of required packages
#rmeta used for the function meta.summaries
#require("rmeta") || stop("`rmeta' package not found")

#Converts the significance level into 'z' from the normal distribution (i.e. 0.05 -> 1.96)
ci <- qnorm(1-((sig.level)/2))

#Number of studies in meta-analysis for use later
length <- length(SS)
length.vector <- c(rep(1,length))


##Calculates the summary effect estimate
meta <- meta.summaries(SS, seSS, method=method, conf.level=(1-sig.level))
tau2 <- meta$tau2



#####CURRENT WEIGHTINGS of studies

if (method=="random") {
   df <- NROW(SS) - 1
   #df2 = the degrees of freedom given a prospective trial is also included
   df2= df+1
   size <- 1/((seSS^2)+tau2)
   }

else size <- 1/(seSS^2)

######Defines the limits of the y axis so that by default the features are intelligently displayed

sediff <- max(seSS) - min(seSS)

if (!is.null(ylim)) if (ylim[1]<ylim[2]) ylim <- rev(ylim)

if (is.null(ylim)) {
   ylim <- c(max(seSS) + 0.20*sediff, min(seSS) - 0.25*sediff)
   if (ylim[2]<0) ylim[2] <- 0
   }

axisdiff <- ylim[2] - ylim[1]

#####CALUCLATES THE X-AXIS LIMITS
SSdiff <- max(SS) - min(SS)

if (is.null(xlim)) {
   xlim <- c(min(SS) - 0.2*SSdiff, max(SS) + 0.2*SSdiff)
   }

#####Creates a vector of weights/sizes so that various contours can be produced
if (contour || !is.null(isquared) || !is.null(tausquared)) {
   cSS <- seq(xlim[1], xlim[2], length.out=contour.points)
   csize <- seq(ylim[1], ylim[2], length.out=contour.points)
   csize[csize<=0] <- 0.0000001*min(seSS)
   for (k in 2:length(csize)) if (csize[k]==0 & csize[k-1]==0) csize[k] <- NA
   csize <- csize[!is.na(csize)]
   }


#####HETEROGENEITY CONTOURS

#####HETEROGENEITY CONTOURS I^2

if (!is.null(isquared)) {

   #calculate current isquared
   Q <- sum((SS^2)/(seSS^2)) - (sum(SS/(seSS^2))^2/sum(1/(seSS^2)))
   isq <- 100*round((Q-length(SS)+1)/Q,digits=3)
   if (isq<0) isq <- 0
   isquared2 <- isquared
   isquared2[isquared2<0] <- isq
   isquared2 <- sort(isquared2)/100

   isqcontourSS <- matrix(rep(1,times=length(isquared)*length(csize)*2),length(csize))
   eval1 <- sum((SS^2)/(seSS^2))
   eval2 <- sum(SS/(seSS^2))
   eval3 <- sum(1/(seSS^2))

   #Two I^2 contours since size -> SS is not a 1 to 1 function
   isqcontourSS1 <- c(rep(1,times=length(csize)))
   isqcontourSS2 <- c(rep(1,times=length(csize)))
   vwt <- 1/(csize^2)

   for (j in 1:length(isquared2)) {

	const <- length/(1-isquared2[j])

	for (i in 1:length(csize))  {
	   if (!is.na(csize[i])) {
		term2 <- ( ((eval3 + vwt[i])/(eval3*vwt[i]))*((eval2^2/(eval3 + vwt[i])) + const - eval1) ) + (eval2/eval3)^2
		if (term2<0) { 
		   isqcontourSS1[i] <- NA    
		   isqcontourSS2[i] <- NA
		   }
		else {
		   term1 <- eval2/eval3
		   isqcontourSS1[i] <- term1+sqrt(term2)
		   isqcontourSS2[i] <- term1-sqrt(term2)

		   #for testing purposes...
		   #Q1 <- eval1 + (isqcontourSS1[i]^2)*vwt[i] - ((eval2 + (isqcontourSS1[i]*vwt[i]))^2/(eval3+vwt[i]))
		   #Q2 <- eval1 + (isqcontourSS2[i]^2)*vwt[i] - ((eval2 + (isqcontourSS2[i]*vwt[i]))^2/(eval3+vwt[i]))
		   #isq1 <- (Q1-length)/Q1
		   #isq2 <- (Q2-length)/Q2
		   }
	      }
	   }

   Q1 <- eval1 + (isqcontourSS1^2)*vwt - ( (eval2 + (isqcontourSS1*vwt))^2/(eval3+vwt) )
   Q2 <- eval1 + (isqcontourSS2^2)*vwt - ( (eval2 + (isqcontourSS2*vwt))^2/(eval3+vwt) )
   isq1 <- (Q1-length(SS))/Q1
   isq2 <- (Q2-length(SS))/Q2
	

      isqcontourSS[,j] <- isqcontourSS1
      isqcontourSS[,j+length(isquared2)] <- isqcontourSS2

	#Displays warning messages when some of the I^2 contours dont show - one study can only change I^2 by limited amount
	isqcontourperc <- round((length(isqcontourSS1[is.na(isqcontourSS1)])/length(isqcontourSS1))*100, digits=2)
	if (isqcontourperc>0) 
		cat("I^2=",isquared2[j]*100," could not be displayed","\n")

	}

###calculates the minimum Isq
#calculates the minimum I^2 which is calculated from the pooled effect of a fixed effects meta analysis
metafixed <- meta.summaries(SS, seSS, method="fixed", conf.level=(1-sig.level))

#the number 1 used is purely arbitrary - same minimum Isq no matter what s.e.
Qmin <- eval1 + (metafixed$summary^2) - ( (eval2 + (metafixed$summary))^2/(eval3+1) )
isqmin <- 100*round((Qmin-length(SS))/Qmin,digits=5)

cat("Minimum I^2 = ", isqmin ,"\n")

}

#####HETEROGENEITY CONTOURS Tau^2
if (!is.null(tausquared)) {

   tausquared2 <- tausquared
   tausquared2[tausquared2<0] <- round(tau2,digits=3)
   tausquared2 <- sort(tausquared2)

   tausqcontourSS <- matrix(rep(1,times=length(tausquared2)*length(csize)*2),length(csize))

   #Two Tau^2 contours since size -> SS is not a 1 to 1 function
   tausqcontourSS1 <- c(rep(1,times=length(csize)))
   tausqcontourSS2 <- c(rep(1,times=length(csize)))

   evalv=sum((1/(seSS^2))^2)
   evalx=sum(SS*(1/(seSS^2)))
   evaly=sum(1/(seSS^2))
   evalz=sum((1/(seSS^2))*(SS)^2)
   vwt <- 1/(csize^2)

   for (j in 1:length(tausquared2)) {

   	for (i in 1:length(csize)) {
	if (!is.na(csize[i])) {
	SStest1 <- evalx/evaly
	SStest2 <- ((evaly + vwt[i])/(evaly*vwt[i])) * 
	   (evaly + vwt[i] -((evalv + vwt[i]^2)/(evaly + vwt[i]))) * tausquared2[j]
	SStest3 <- ((evaly + vwt[i])/(evaly*vwt[i])) * (evalz - (evalx^2/(evaly + vwt[i])) - length(SS))

	eval.sqrt <- SStest2 + (SStest1^2) - SStest3
	
	if (eval.sqrt>=0) {
	   tausqcontourSS1[i] <- SStest1 + sqrt(SStest2 + (SStest1^2) - SStest3)
	   tausqcontourSS2[i] <- SStest1 - sqrt(SStest2 + (SStest1^2) - SStest3)
	   }
	else {
	   tausqcontourSS1[i] <- NA
	   tausqcontourSS2[i] <- NA
	   }

	}
	}

   tausqcontourSS[,j] <- tausqcontourSS1
   tausqcontourSS[,j+length(tausquared)] <- tausqcontourSS2

   #Displays warning messages when some of the I^2 contours dont show - one study can only change I^2 by limited amount
   tausqcontourperc <- round((length(tausqcontourSS1[is.na(tausqcontourSS1)])/length(tausqcontourSS1))*100, digits=2)
   if (tausqcontourperc>0) 
	cat(tausqcontourperc, "% of Tau^2=",tausquared2[j]," could not be displayed","\n")

   }

###calculates the minimum Tau sq
#calculates the minimum tau^2 which is calculated from the pooled effect of a fixed effects meta analysis
metafixed <- meta.summaries(SS, seSS, method="fixed", conf.level=(1-sig.level))

#the number 1 used is purely arbitrary - same minimum Isq no matter what s.e.
metafixed2 <- meta.summaries(c(SS,metafixed$summary), c(seSS,(0.0000001*min(seSS^2))^0.5), 
   method="fixed", conf.level=(1-sig.level))
tausqmin <- metafixed2$tau2
cat("Minimum Tau^2 = ", tausqmin ,"\n")

}


#####Creating the contours
if (contour) {

if (method=="fixed")  {
   vwt <- 1/(csize^2)
   c1SS <- (1/vwt)*(zero*(sum(size) + vwt) - sum(size*SS) +  ci * (sum(size)+vwt)^0.5)
   c2SS <- (1/vwt)*(zero*(sum(size) + vwt) - sum(size*SS) -  ci * (sum(size)+vwt)^0.5)
   }

if (method=="random")  {

   matcont <- matrix(rep(NA,times=length(cSS)*length(csize)),nrow=length(csize))	
   overmax <- 0

   for (i in 1: length(csize))  {

	#display percentage of computation complete
	if (rand.load>0) {
   	   roundi<-i/rand.load
	   flush.console()
   	   if (roundi==round(roundi,0)) {
   	      perc_complete <- (i/contour.points)*100
	      cat(perc_complete, "%")
	   }
	   else cat(".")
	}

	if ( !is.na(csize[i]) ) {

	   for (j in 1:length(cSS))  {

	      metacont <- meta.summaries(c(SS, cSS[j]), c(seSS, csize[i]), method=method, conf.level=(1-sig.level))
	      lc <- metacont$summary - ci*metacont$se.summary
	      uc <- metacont$summary + ci*metacont$se.summary

 	      if (lc < zero & uc < zero) matcont[i,j] <- 1
	      if (lc < zero & uc > zero) matcont[i,j] <- 0
	      if (lc > zero & uc > zero) matcont[i,j] <- 2

	      }	   
	
	   }

	else matcont[i,] <- 3

	}

   }

}


#####Defining Axis labels
if (is.null(xlab)) xlab <- "Effect"
if (is.null(ylab)) ylab <- "Standard Error"

#####Summary diamond
if (summ) {
   xsumm <- c(meta$summary - ci*meta$se.summary, meta$summary, meta$summary + ci*meta$se.summary, meta$summary)
   ysumm <- c(ylim[2]-0.10*axisdiff+summ.pos,ylim[2]-0.07*axisdiff+summ.pos,ylim[2]-0.10*axisdiff+summ.pos,ylim[2]-0.13*axisdiff+summ.pos)

   if (pred.interval) {	
	predint1 <- meta$summary - qt(p=sig.level,df=length-2)*(meta$tau2+(meta$se.summary^2))^0.5
	predint2 <- meta$summary + qt(p=sig.level,df=length-2)*(meta$tau2+(meta$se.summary^2))^0.5
	}
}


######CREATING THE LEGEND
#creates an inital matrix (which will be modified) containing all information needed to define the legend if TRUE

legendmat1 <- c(rep(NA,times=13))
legendmat2 <- c(rep(NA,times=13))
legendmat3 <- c(rep(NA,times=13))
legendmat4 <- c(rep(NA,times=13))

if (contour) {
#1

   legendmat1[1] <- paste("Non Sig Effect (",sig.level*100,"% level)", sep="")
   legendmat2[1] <- "black"
   legendmat3[1] <- 0
   legendmat4[1] <- 0

#2
   legendmat1[2] <- paste("Sig Effect>NULL (",sig.level*100,"% level)", sep="")
   legendmat2[2] <- "gray72"
   legendmat3[2] <- 0
   legendmat4[2] <- 15

#3
   legendmat1[3] <- paste("Sig Effect<NULL (",sig.level*100,"% level)", sep="")
   legendmat2[3] <- "gray91"
   legendmat3[3] <- 0
   legendmat4[3] <- 15
   }

#4
if (plot.zero) {
   legendmat1[4] <- "Null Effect"
   legendmat2[4] <- "lightgrey"
   legendmat3[4] <- 1
   legendmat4[4] <- 46
   }

#5
if (plot.summ) {
   legendmat1[5] <- "Pooled effect"
   legendmat2[5] <- "slategrey"
   legendmat3[5] <- 1
   legendmat4[5] <- 46
   }

#6
if (pred.interval) {
   legendmat1[6] <- paste((1-sig.level)*100,"% Pred Interval", sep="")
   legendmat2[6] <- "black"
   legendmat3[6] <- 1
   legendmat4[6] <- 46
   }

#7-12
if (!is.null(isquared) || !is.null(tausquared)) {
   legendmat1[7] <- ""
   legendmat2[7] <- "white"
   legendmat3[7] <- 0
   legendmat4[7] <- 46
   if (!is.null(tausquared)) legendmat1[8] <- "Tau^2 Contours:"
   else legendmat1[8] <- "I^2 Contours:"
   legendmat2[8] <- "white"
   legendmat3[8] <- 0
   legendmat4[8] <- 46

   if (!is.null(isquared)) het <- paste(isquared2*100,"%")
   if (!is.null(tausquared)) het <- tausquared2

     for (i in 1:4) {
	if (length(het)>i-1)  {
	   legendmat1[i+8] <- het[i]
	   if (!is.null(isquared)) { 
		if (het[i]==paste(isq,"%")) {
		   legendmat1[i+8] <- paste(het[i],"(current)")
		   legendmat2[i+8] <- "black"
		   legendmat3[i+8] <- 1
		   }
		if (het[i]!=paste(isq,"%")) {
		   if (i==1) legendmat2[i+8] <- "sienna1"
	   	   if (i==2) legendmat2[i+8] <- "olivedrab"
	   	   if (i==3) legendmat2[i+8] <- "deepskyblue2"
	   	   if (i==4) legendmat2[i+8] <- "firebrick2"
	   	   if (i==1) legendmat3[i+8] <- 5
	   	   if (i==2) legendmat3[i+8] <- 2
	   	   if (i==3) legendmat3[i+8] <- 4
	   	   if (i==4) legendmat3[i+8] <- 3
	   	   if (i==1) legendmat4[i+8] <- 46
	   	   if (i==2) legendmat4[i+8] <- 46
	   	   if (i==3) legendmat4[i+8] <- 46
	   	   if (i==4) legendmat4[i+8] <- 46
		   }
		}

	   if (!is.null(tausquared)) {
		if (het[i]==round(tau2,digits=3)) {
		   legendmat1[i+8] <- paste(het[i],"(current)")
		   legendmat2[i+8] <- "black"
		   legendmat3[i+8] <- 1
		   }
		if (het[i]!=round(tau2,digits=3)) {
	   	   if (i==1) legendmat2[i+8] <- "sienna1"
	   	   if (i==2) legendmat2[i+8] <- "olivedrab"
	   	   if (i==3) legendmat2[i+8] <- "deepskyblue2"
	   	   if (i==4) legendmat2[i+8] <- "firebrick2"
	   	   if (i==1) legendmat3[i+8] <- 5
	   	   if (i==2) legendmat3[i+8] <- 2
	   	   if (i==3) legendmat3[i+8] <- 4
	  	   if (i==4) legendmat3[i+8] <- 3
	  	   if (i==1) legendmat4[i+8] <- 46
	  	   if (i==2) legendmat4[i+8] <- 46
	   	   if (i==3) legendmat4[i+8] <- 46
	   	   if (i==4) legendmat4[i+8] <- 46
               }
		}
	   }
	}
   }

#summary diamond
if (summ) {
   legendmat1[13] <- "Pooled result"
   legendmat2[13] <- "lavenderblush4"
   legendmat3[13] <- 0
   legendmat4[13] <- 18
   }

#####CREATING THE PLOT
if (legend) {
   dev.new(width=15, height=11)
   par(mai = c(1, .8, .2, 3))
   par(plt=c(0.0841, 0.6967, 0.13, 0.95))
}
else {
   dev.new(width=11, height=11)
   par(mai = c(1, .8, .2, .2))
   par(plt=c(0.13, 0.95, 0.13, 0.95))
}



#default axis is supressed if x axis values specified in options (exp or otherwise)
if (is.null(expxticks) & is.null(xticks)) xaxis <- "s" else xaxis <- "n"

#default axis is supressed if y axis values specified in options
if (is.null(yticks)) yaxis <- "s" else yaxis <- "n"

if (!points) cexpoints <- 0
else cexpoints <- 0.8

#The intial plot of individual points
plot(SS, seSS, ylim=ylim, xlim=xlim, xlab = xlab, ylab = ylab, pch=19, cex=cexpoints, col= "black", 
   xaxt=xaxis, yaxt=yaxis, xaxs="i", yaxs="i")

#The significance contours
if (contour & method=="random") {
   for (j in 1:(length(cSS)-1)) {
	for (i in 1:(length(csize)-1)) {
   	   if (matcont[i,j]==0 & matcont[i+1,j]==0 & matcont[i,j+1]==0 & matcont[i+1,j+1]==0)
		polygon(c(cSS[j],cSS[j],cSS[j+1],cSS[j+1]),c(csize[i],csize[i+1],csize[i+1],csize[i]), 
		   border="white", col = "white")
   	   if (matcont[i,j]==1 & matcont[i+1,j]==1 & matcont[i,j+1]==1 & matcont[i+1,j+1]==1)
		polygon(c(cSS[j],cSS[j],cSS[j+1],cSS[j+1]),c(csize[i],csize[i+1],csize[i+1],csize[i]), 
		   border="gray91", col = "gray91")
   	   if (matcont[i,j]==2 & matcont[i+1,j]==2 & matcont[i,j+1]==2 & matcont[i+1,j+1]==2)
		polygon(c(cSS[j],cSS[j],cSS[j+1],cSS[j+1]),c(csize[i],csize[i+1],csize[i+1],csize[i]), 
		   border="gray72", col = "gray72")
	   }
	}
}

if (contour & method=="fixed") {
polygon(c(c1SS,rev(c2SS)), c(csize, rev(csize)), border="white", col = "white")
polygon(c(c2SS,xlim[1], xlim[1]), c(csize, ylim[2], ylim[1]), border="gray91", col = "gray91")
polygon(c(c1SS,xlim[2], xlim[2]), c(csize, ylim[2], ylim[1]), border="gray72", col = "gray72")

}

#Summary Diamond
if (summ) {
   if (pred.interval) {
	segments(x0=predint1, y0=ysumm[1], x1=predint2, y1=ysumm[1])
	}

   polygon(xsumm,ysumm, border="black", col = "lavenderblush4")
   }

#Pooled effect line
if (plot.summ) abline(v = meta$summary, col="slategrey")

#Null vertical line
if (plot.zero) abline(v = zero, col="lightgrey", lty=1)

#Adding x axis ticks for exponential effects
if (!is.null(expxticks)) axis(1, at=log(expxticks), labels=expxticks)
#Adding x axis ticks for specified (non-exponential effects)
if (!is.null(xticks)) axis(1, at=xticks, labels=xticks)
#Adding y axis ticks for pre-specified standard errors
if (!is.null(yticks)) axis(2, at=yticks, labels=yticks)

##########Heterogeneity contours
###Isquared
if (!is.null(isquared)) {
   colour <- legendmat2[9:12]
   colour <- colour[!is.na(colour)]
   lty <- legendmat3[9:12]
   lty <- lty[!is.na(lty)]

   for (i in 1:(length(isquared2))) {
      points(isqcontourSS[!is.na(isqcontourSS[,i]),i], csize[!is.na(isqcontourSS[,i])], 
	   type="l", col=colour[i], lty=lty[i], lwd=2)
   }

   for (i in length(isquared2):(2*length(isquared2))) {
      points(isqcontourSS[!is.na(isqcontourSS[,i]),i], csize[!is.na(isqcontourSS[,i])], 
	   type="l", col=colour[i-length(isquared2)], lty=lty[i-length(isquared2)], lwd=2)
   }
}

###Tausquared
if (!is.null(tausquared)) {
   colour <- legendmat2[9:12]
   colour <- colour[!is.na(colour)]
   lty <- legendmat3[9:12]
   lty <- lty[!is.na(lty)]

   for (i in 1:(length(tausquared2))) {
      points(tausqcontourSS[!is.na(tausqcontourSS[,i]),i], csize[!is.na(tausqcontourSS[,i])], 
	   type="l", col=colour[i], lty=lty[i], lwd=2)
   }

   for (i in length(tausquared2):(2*length(tausquared2))) {
      points(tausqcontourSS[!is.na(tausqcontourSS[,i]),i], csize[!is.na(tausqcontourSS[,i])], 
	   type="l", col=colour[i-length(tausquared2)], lty=lty[i-length(tausquared2)], lwd=2)
   }
}


#Features may have ovelayed axes, this part redraws axes
box()

#Overlay points
points(SS, seSS, ylim=ylim, xlim=xlim, xlab = xlab, ylab = ylab, pch=19, cex=cexpoints, col= "black", xaxt=xaxis)

#Plots the legend
if (legend) {
   par(xpd=TRUE)
   legendmat1 <- legendmat1[!is.na(legendmat1)]
   legendmat2 <- legendmat2[!is.na(legendmat2)]
   legendmat3 <- legendmat3[!is.na(legendmat3)]
   legendmat4 <- legendmat4[!is.na(legendmat4)]
   if (length(legendmat1)>0)  legend(legendpos[1], legendpos[2],legendmat1, col = legendmat2, lty = legendmat3, 
	pch=legendmat4, pt.cex=2.5)
   }

if (!is.null(xpoints) & !is.null(ypoints)) {
points(xpoints, ypoints, pch=21, cex=4, col= "black")
points(xpoints, ypoints, pch=4, cex=3, col= "black")
}

}
