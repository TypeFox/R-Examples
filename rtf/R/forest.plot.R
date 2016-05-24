#########################################################################/**
# @RdocFunction rtf.forest.plot
#
# @title "Get an RTF encoded forest plot"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{x}{x (e.g. hazard ratio).}
# 	\item{min}{Minimum whisker (e.g. lower bound of 95\% hazard ratio CI).}
# 	\item{max}{Maximum whisker (e.g. upper bound of 95\% hazard ratio CI).}
# 	\item{xlim}{A @vector specifying the x limits.}
# 	\item{width}{Plot width in inches.}
# 	\item{height}{Plot height in inches.}
# 	\item{cex}{A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.}
# 	\item{lwd}{Line width.}
# 	\item{res}{Output resolution in dots per inch.}
# }
#
# \details{
# 	Create a forest plot and convert PNG to RTF code.  This is useful for
# 	embedding into a data frame of hazard ratios and then writing an
# 	RTF output file.  See the example below for usage.
# }
#
# \examples{
# \dontrun{
# tab<-data.frame(
#	Label=c("Test1","Test2","Test3"),
# 	HR=c(1,2,0.45),
# 	Lower.CI=c(0.5,1.1,0.25),
# 	Upper.CI=c(2,3.5,0.9), 
# 	stringsAsFactors=FALSE, 
# 	check.names=FALSE)
#
# # create forest plots by row
# forest.plot.args<-list(xlim=c(0.1,5),width=3.0,height=0.3,cex=1,lwd=0.75,res=300)
# tab$"HR Plot (log scale)"<-mapply(rtf.forest.plot,tab$HR,tab$Lower.CI,tab$Upper.CI,
#			MoreArgs=forest.plot.args)
#
# # rbind the x-scale to the table in the plot column
# xscale<-rtf.forest.plot.xscale(xlim=c(0.1,5),width=3.0,height=0.3,cex=1,
# 			lwd=0.75,res=300)
#
# tab<-data.frame(lapply(tab, as.character), 
# 			stringsAsFactors=FALSE, 
#			check.names=FALSE)
#
# tab<-rbind(tab,list("","","","",xscale))
#
# # write the RTF output
# rtf<-RTF("test_rtf.forest.plot.doc",width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
# addTable(rtf,tab,col.widths=c(0.75,0.75,0.75,0.75,3))
# done(rtf)
# }
# }
#*/#########################################################################
rtf.forest.plot<-function(x=1.25,min=0.5,max=2,xlim=c(0.1,12),
	width=3.0,height=0.3,cex=1,lwd=0.75,res=300) {
	
	tmp.file<-tempfile("temp_forest_plot")
	ret<-.rtf.plot(.forest.plot,tmp.file=tmp.file,width=width,height=height,res=res,x=x,min=min,max=max,xlim=xlim,cex=cex,lwd=lwd)
	if(file.exists(tmp.file) ) {
		unlink(tmp.file)
	}
	ret
}

rtf.forest.plot.xscale<-function(xlim=c(0.1,12),
	width=3.0,height=0.3,cex=1,lwd=0.75,res=300) {
	tmp.file<-tempfile("temp_forest_xscale")
	ret<-.rtf.plot(.forest.plot.scale,tmp.file=tmp.file,width=width,height=height,res=res,xlim=xlim,cex=cex,lwd=lwd)
	if(file.exists(tmp.file) ) {
		unlink(tmp.file)
	}
	ret
}

######################################################################################

.forest.plot<-function(x=1.25,min=0.5,max=2,xlim=c(0.1,12),cex=1,lwd=0.75) {
	par(oma=c(0,0,0,0),mar=c(0,0,0,0))
	plot(c(1,1),type="n",axes=FALSE,xlab="",ylab="",xlim=xlim,ylim=c(1,1),log="x",xaxt="n",yaxt="n",bty="l")

	abline(v=1,lwd=0.5,col="black") # tick marks at '1' behind the plot glyphs

	if(!is.na(min)) { arrows(x,1,min,1,angle=90,lwd=lwd,length=0.05) }
	if(!is.na(max)) { arrows(x,1,max,1,angle=90,lwd=lwd,length=0.05) }	
	if(!is.na(x)) { points(x,1,pch=21,bg="black",cex=cex) }
}

.forest.plot.scale<-function(xlim=c(0.1,12),cex=1,lwd=0.75,res=300) {
	par(oma=c(0,0,0,0),mar=c(0,0,0,0))
	plot(c(1,1),type="n",axes=FALSE,xlab="",ylab="",xlim=xlim,ylim=c(0.95,1.05),log="x",xaxt="n",yaxt="n",bty="l")
	ticks<-axTicks(1)
	y.max=1.025
	sapply(ticks,function(x,y.max){lines(x=c(x,x),y=c(1,y.max),lwd=lwd)},y.max=y.max)
	lines(c(ticks[1],ticks[length(ticks)]),c(y.max,y.max),lwd=lwd)
	text(ticks,1,ticks,pos=1,offset=0.25,xpd=NA,cex=cex)
}
