# fclass.R
# 
# function class for easy plotting
#
# Note 060909: From 2.10 on this class is obsolete, we use spatstat's fv-class!
#
# Author: Tuomas Rajala <tarajala@maths.jyu.fi>
###############################################################################


#segfcl<-function(x)
#{
#	class(x)<-'segfcl'
#}
#
#print.segfcl<-function(x,...)
#{
#	p<-paste("(",x$gtype)
#	if(length(x$par)< 2){if(x$par!=0)p<-paste(p,", par=",x$par)}
#	else p<-paste(p,", range (",paste(range(x$par),collapse=","),") ",sep="")
#	
#	p<-paste(p,")")
#	cat(paste(x$note2, p,"\n"))
#	if(x$note!="")cat(paste("Note:",x$note,"\n"))
##	cat(names(x)[1]);cat("=\n")
##	print(x[[1]])
#	
#	cat(spatialsegregation_function_explanations[[names(x)[[1]]]])
#	
#}
#
#
#plot.segfcl<-function(x, ...)
#{
#	par<-x[['par']]
#	f<-x[[1]]
#	plot(par, f, ...)
#	p<-x[['poisson']]
#	if(length(p)<2) abline(h=p,col=2,lty=2)
#	else lines(par,p,col=2,lty=2)
#}
#
#
##########
#spatialsegregation_function_explanations<-list(
#		M="$M : means over typewise mingling values\n$Mw: weighted means, relative abundances as weights\n$sd : sd of typewise values\n$par : parameters used\n$typewise : typewise mingling values\n$poisson : Value for Poisson\n",
#		S="$S : spatial Simpson index values\n$par : parameters used\n$typewise : typewise values of which S is the sum\n$poisson : Value for Poisson\n",
#		H="$H : spatial Shannon index values\n$par : parameters used\n$typewise : pi_tau values i.e. local relative frequencies\n$global : Global entropy\n$poisson : Value for Poisson\n",
#		I="$I : means over typewise ISAR values\n$Iw: weighted means, relative abundances as weights\n$sd : sd of typewise values\n$par : parameters used\n$typewise : typewise ISAR values\n$poisson : Values for Poisson\n"
#		)
#
########################################################################################
