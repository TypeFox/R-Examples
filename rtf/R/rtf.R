##########################################################################################
## RTF Output Functions for R                                                            #
##                                                                                       # 
## Author: Michael E. Schaffer, Ph.D.                                                    #                                                                       #
##                                                                                       #
## Description:                                                                          #
## A set of R functions to output RTF files with high resolution graphics and tables.    #
## This is useful for reporting R results to an Microsoft Word-compatible report.  All   #
## graphics must be in a format supported by Word.  Currently the most compatible format #
## is PNG.                                                                               #
##                                                                                       #
## For details about the RTF format (a Microsoft format), see:                           #
## http://latex2rtf.sourceforge.net/rtfspec_7.html#rtfspec_paraforprop                   #
## http://www.pindari.com/rtf2.html                                                      #
##                                                                                       #
## For use as source file include: require(R.oo)                                         #
##########################################################################################



#########################################################################/**
# @RdocClass RTF
#
# @title "The RTF class"
#
# \description{
#	This is the class representing an RTF file output.
#
#	@classhierarchy
# }
#
# @synopsis
#
# \arguments{
# 	\item{file}{The path of the output file.}
# 	\item{width}{The width of the output page.}
# 	\item{height}{The width of the output page.}
# 	\item{omi}{A @vector representing the outer margins in inches (bottom, left, top, right).}
#	\item{font.size}{Default font size for the document in points.}
# 	\item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
# 	@allmethods
# }
#
# \examples{
# \dontrun{
# output<-"test_RTF-class.doc"
# png.res<-300
#
# rtf<-RTF(output,width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
# addHeader(rtf,title="Test",subtitle="2011-08-15\n")
# addPlot(rtf,plot.fun=plot,width=6,height=6,res=300, iris[,1],iris[,2])
# 
# # Try trellis plots
# if(require(lattice) & require(grid)) {
# 	# single page trellis objects
# 	addPageBreak(rtf, width=11,height=8.5,omi=c(0.5,0.5,0.5,0.5))
# 
# 	p <- histogram( ~ height | voice.part, data = singer, xlab="Height")
# 	addTrellisObject(rtf,trellis.object=p,width=10,height=7.5,res=png.res)
# 
# 	p <- densityplot( ~ height | voice.part, data = singer, xlab = "Height")
# 	addTrellisObject(rtf,trellis.object=p,width=10,height=7.5,res=png.res)
# 
# 	# multipage trellis object
# 	p2<-xyplot(uptake ~ conc | Plant, CO2, layout = c(2,2))
# 	addTrellisObject(rtf,trellis.object=p2,width=6,height=6,res=png.res)
# }
# 
# addPageBreak(rtf, width=6,height=10,omi=c(0.5,0.5,0.5,0.5))
# addTable(rtf,as.data.frame(head(iris)),font.size=10,row.names=FALSE,NA.string="-")
#
# addSessionInfo(rtf)
# 
# done(rtf)
# }
# }
#
# @author
#
# \seealso{
# 	@seeclass
# }
#*/#########################################################################
setConstructorS3("RTF", 
	function(file="",width=8.5,height=11,omi=c(1,1,1,1),font.size=10) {
	this <- extend(Object(), "RTF", 
		.rtf = .start.rtf(width,height,omi), 
		.file = file,
		.font.size = font.size,
		.indent = 0,
		.page.width = width,
		.page.height = height,
		.content.width = width - omi[2] - omi[4]
	);
	
	this;
});


###########################################################################/**
# @RdocMethod addTable
#
# @title "Insert a table into the RTF document"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{dat}{A matrix, data frame, or table.}
#	\item{col.widths}{A @vector of column widths in inches. \bold{optional}.}
#	\item{col.justify}{A single value or @vector of column justifications ('L', 'R', 'C', or 'J' for Left, Right, Center, and Justify, respectively). \bold{optional}.}
#	\item{header.col.justify}{A single value or @vector of table header column justifications ('L', 'R', 'C', or 'J' for Left, Right, Center, and Justify, respectively). \bold{optional}.}
#	\item{font.size}{Font size in points. \bold{optional}.}
#	\item{row.names}{Boolean argument to include row names in tables. \bold{optional}.}
#	\item{NA.string}{A character string to replace NA values in the table.}
#	\item{space.before}{Space before each row (in inches). \bold{optional}.}
#	\item{space.after}{Space after each row (in inches). \bold{optional}.}
# 	\item{...}{Not used.}
# }
#
# \examples{
# rtf<-RTF("test_addTable.doc",width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
# addTable(rtf,as.data.frame(head(iris)),font.size=10,row.names=FALSE,NA.string="-",
#          col.widths=rep(1,5))
# 
# tab<-table(iris$Species,floor(iris$Sepal.Length))
# names(dimnames(tab))<-c("Species","Sepal Length")
# addTable(rtf,tab,font.size=10,row.names=TRUE,NA.string="-",col.widths=c(1,rep(0.5,4)) )
#
# done(rtf)
# }
#
# @author
#
# \seealso{
# 	@seeclass
# }
# 
# @keyword -internal
#*/###########################################################################
setMethodS3("addTable", "RTF", function(this,dat,col.widths=NULL,col.justify=NULL,header.col.justify=NULL,font.size=NULL,row.names=FALSE,NA.string="-", space.before=NULL, space.after=NULL, ...) {
	if(is.null(font.size)) {
		font.size = this$.font.size  # default
	}
	
	this$.rtf <- paste(this$.rtf, .add.table(dat,col.widths,col.justify,header.col.justify,font.size,row.names,indent=this$.indent,NA.string=NA.string,max.table.width=this$.content.width, space.before=space.before, space.after=space.after, ...),sep="")
})

#########################################################################/**
# @RdocMethod view
#
# @title "View encoded RTF"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{...}{Not used.}
# }
#
# \value{
# 	Output the content of the object as RTF code.
# }
#
# @author
#
# \seealso{
# 	@seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("view", "RTF", function(this, ...) {
	print(this$.rtf)
	print(this$.file)
})

#########################################################################/**
# @RdocMethod done
#
# @title "Write and close the RTF output"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{...}{Not used.}
# }
#
# @author
#
# \seealso{
# 	@seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("done", "RTF", function(this, ...) {
	#this$.rtf <- paste(this$.rtf,.end.rtf(),sep="")
	#write(this$.rtf,this$.file)
	
	# write the file, but don't close it out so that we can continue to 
	# add to the object and write it out again.
	write(paste(this$.rtf,.end.rtf(),sep=""),this$.file)
})


#########################################################################/**
# @RdocMethod addTOC
#
# @title "Insert table of contents field"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("addTOC", "RTF", function(this,...) {
	toc<-"{\\field\\flddirty\\fldedit{\\*\\fldinst TOC f h}{\\fldrslt Update Field (right-click in MS Word) to show Table of Contents}}\\line\\line"
	
	if(!is.null(this$.font.size)) {
		font.size = this$.font.size  # default
	}
	
	this$.rtf <- paste(this$.rtf,.start.paragraph(indent=this$.indent,font.size=font.size),sep="")
	this$.rtf <- paste(this$.rtf,toc,sep="")
	this$.rtf <- paste(this$.rtf,.end.paragraph(),sep="")
})


#########################################################################/**
# @RdocMethod addHeader
#
# @title "Insert a header into the RTF document"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
#   \item{title}{Header title text.}
#   \item{subtitle}{Header subtitle text. \bold{optional}.}
#	\item{font.size}{Font size in points. \bold{optional}.}
#	\item{TOC.level}{Indent level for table of contents. \bold{optional}.}
# 	\item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("addHeader", "RTF", function(this, title,subtitle=NULL,font.size=NULL,TOC.level=NULL,...) {
	if(is.null(font.size)) {
		font.size = this$.font.size  # default
	}
	
	this$.rtf <- paste(this$.rtf,.add.header(title,subtitle=subtitle,indent=this$.indent,font.size=font.size,TOC.level=TOC.level),sep="")
})



#########################################################################/**
# @RdocMethod addText
#
# @title "Insert text into the RTF document"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{bold}{Bold text. \bold{optional}.}
# 	\item{italic}{Italic text. \bold{optional}.}
# 	\item{...}{Any number of character strings to concatenate.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("addText", "RTF", function(this, ..., bold=FALSE, italic=FALSE) {
	text<-paste(... , sep="")
	if(bold) {
		text<-paste("\\b ",text,"\\b0",sep="")
	}
	if(italic){
		text<-paste("\\i ",text,"\\i0",sep="")
	}
	this$.rtf <- paste(this$.rtf,.add.text(text),sep="")
})

#########################################################################/**
# @RdocMethod addParagraph
#
# @title "Insert a paragraph into the RTF document"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{...}{A character @vector of text to add.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("addParagraph", "RTF", function(this, ...) {
	text<-paste(... , sep="")
	
	if(!is.null(this$.font.size)) {
		font.size = this$.font.size  # default
	}
	
	this$.rtf <- paste(this$.rtf,.start.paragraph(indent=this$.indent,font.size=font.size),sep="")
	this$.rtf <- paste(this$.rtf,.add.text(text),sep="")
	this$.rtf <- paste(this$.rtf,.end.paragraph(),sep="")
})

#########################################################################/**
# @RdocMethod startParagraph
#
# @title "Start a new paragraph in the RTF document"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("startParagraph", "RTF", function(this, ...) {
	this$.rtf <- paste(this$.rtf,.start.paragraph(indent=this$.indent),sep="")
})

#########################################################################/**
# @RdocMethod endParagraph
#
# @title "End a paragraph in the RTF document"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("endParagraph", "RTF", function(this, ...) {
	this$.rtf <- paste(this$.rtf,.end.paragraph(),sep="")
})

###########################################################################/**
# @RdocMethod addPageBreak
#
# @title "Insert a page break into the RTF document optionally specifying new page settings"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{width}{New page width in inches. \bold{optional}.}
# 	\item{height}{New page height in inches. \bold{optional}.}
# 	\item{omi}{A @vector of page margins (botton, left, top, right) \bold{optional}.}
# 	\item{...}{Not used.}
# }
#
# \examples{
# rtf<-RTF("test_addPageBreak.doc",width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
# addPageBreak(rtf,width=11,height=8.5,omi=c(0.5,0.5,0.5,0.5))
# done(rtf)
# }
#
# @author
#
# \seealso{
# 	@seeclass
# }
# 
# @keyword -internal
#*/###########################################################################
setMethodS3("addPageBreak", "RTF", function(this, width=8.5,height=11,omi=c(1,1,1,1), ...) {
	this$.rtf <- paste(this$.rtf,.add.page.break(width=width,height=height,omi=omi),sep="")	
	this$.page.width = width
	this$.page.height = height
	this$.content.width = width - omi[2] - omi[4]
})

#########################################################################/**
# @RdocMethod addNewLine
#
# @title "Insert a new line into the RTF document"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{n}{Number of lines to add. Default is 1. \bold{optional}}
# 	\item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("addNewLine", "RTF", function(this, n=1, ...) {
	this$.rtf <- paste(this$.rtf,.add.newline(n=n,font.size=this$.font.size),sep="")
})

#########################################################################/**
# @RdocMethod increaseIndent
#
# @title "Increase RTF document indent"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("increaseIndent", "RTF", function(this, ...) {
	this$.indent <- this$.indent + 720 # 1/2" increments
})

#########################################################################/**
# @RdocMethod decreaseIndent
#
# @title "Decrease RTF document indent"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("decreaseIndent", "RTF", function(this, ...) {
	this$.indent <- max(0,this$.indent - 720) # 1/2" increments
})

#########################################################################/**
# @RdocMethod setFontSize
#
# @title "Set RTF document font size"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{font.size}{New default font size in points.}
# 	\item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("setFontSize", "RTF", function(this, font.size, ...) {
	this$.font.size <- font.size
})

#########################################################################/**
# @RdocMethod addPlot
#
# @title "Insert a plot into the RTF document"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{plot.fun}{Plot function.}
# 	\item{width}{Plot output width in inches.}
# 	\item{height}{Plot output height in inches.}
# 	\item{res}{Output resolution in dots per inch.}
# 	\item{...}{Arguments for \code{plot.fun}.}
# }
#
# \details{
# 	Plots are added to the document as PNG objects.  This function will work with all
#   base graphics methods for plotting.  For more sophisticated plots, you may need to
#   wrap your plot code into a function, and then pass a reference to that function to
#   this method.  The parameters for the plot method would then get passed in as '...'
#   above.
#
#   To output a ggplot2 plot, simply assign the plot to a variable.  Then use 'print'
#   as the plot function and pass in the plot variable assigned above. 
# }
#
# \examples{
# rtf<-RTF("test_addPlot.doc",width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
# addPlot(rtf,plot.fun=plot,width=6,height=6,res=300, iris[,1],iris[,2])
# done(rtf)
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("addPlot", "RTF", function(this,plot.fun=plot.fun,width=3.0,height=0.3,res=300, ...) {
	if(!is.null(this$.font.size)) {
		font.size = this$.font.size  # default
	}
	
	tmp.file<-tempfile("temp_rtf_plot")
	this$.rtf <- paste(this$.rtf,.start.paragraph(indent=this$.indent,font.size=font.size),sep="")
	this$.rtf <- paste(this$.rtf,.rtf.plot(plot.fun=plot.fun,tmp.file=tmp.file,width=width,height=height,res=res, ...),sep="")
	this$.rtf <- paste(this$.rtf,.end.paragraph(),sep="")
	
	if(file.exists(tmp.file) ) {
		unlink(tmp.file)
	}
})

#########################################################################/**
# @RdocMethod addPng
#
# @title "Insert an existing PNG image into the RTF document"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{file}{Image file path.}
# 	\item{width}{Plot output width in inches.}
# 	\item{height}{Plot output height in inches.}
# 	\item{...}{Not used.}
# }
#
# \details{
# 	Add existing PNG file to RTF document.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("addPng", "RTF", function(this,file,width=3.0,height=0.3, ...) {
	if(!is.null(this$.font.size)) {
		font.size = this$.font.size  # default
	}
	
	this$.rtf <- paste(this$.rtf,.start.paragraph(indent=this$.indent,font.size=font.size),sep="")
	this$.rtf <- paste(this$.rtf,.add.png(file,width=width,height=height,verbose=FALSE),sep="")
	this$.rtf <- paste(this$.rtf,.end.paragraph(),sep="")
})

#########################################################################/**
# @RdocMethod addTrellisObject
#
# @title "Insert a trellis plot object into the RTF document"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
# 	\item{trellis.object}{The trellis plot object.}
# 	\item{width}{Plot output width in inches.}
# 	\item{height}{Plot output height in inches.}
# 	\item{res}{Output resolution in dots per inch.}
# 	\item{rotate}{Object rotation in degrees. \bold{optional}.}
# 	\item{...}{Not used.}
# }
#
# \details{
# 	Plots are added to the document as PNG objects.  Multi-page trellis objects are 
#	automatically split across multiple pages in the RTF output file.  To rotate the
#   object to landscape orientation within the RTF output, use rotate=90.  When using 
#   rotation, width and height still refer to the unrotated plot dimensions and not the 
#   rotated output dimensions on the RTF page.  An alternative to rotating the plot is
#   to rotate the entire page using a call to addPageBreak with suitable page width and
#   height dimensions.
# }
#
# \examples{
# \dontrun{
# rtf<-RTF("test_addTrellisObject.doc",width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
# if(require(lattice) & require(grid)) {
# 	# multipage trellis object
# 	p2<-xyplot(uptake ~ conc | Plant, CO2, layout = c(2,2))
# 	addTrellisObject(rtf,trellis.object=p2,width=8,height=4,res=300, rotate=90)
# }
# done(rtf)
# }
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("addTrellisObject", "RTF", function(this,trellis.object,width=3.0,height=0.3,res=300, rotate=NULL, ...) {
	tmp.file<-tempfile("temp_rtf_trellis")
	this$.rtf <- paste(this$.rtf,.rtf.trellis.object(trellis.object=trellis.object,tmp.file=tmp.file,width=width,height=height,res=res,rotate=rotate),sep="")
	if(file.exists(tmp.file) ) {
		unlink(tmp.file)
	}
})


#########################################################################/**
# @RdocMethod addSessionInfo
#
# @title "Insert session information into the RTF document"
#
# \description{
#	@get "title".
# }
#
# @synopsis
#
# \arguments{
# 	\item{this}{An RTF object.}
#	\item{locale}{Output the locale.}
# 	\item{...}{Not used.}
# }
#
# \details{
# 	Exports session information to the RTF document in a similar
# }
#
# \examples{
# rtf<-RTF("test_addSessionInfo.doc",width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
# addSessionInfo(rtf)
# done(rtf)
# }
#
# @author
#
# \seealso{
#   @seeclass, \code{\link{sessionInfo}}.
# }
# 
# @keyword -internal
#*/#########################################################################
setMethodS3("addSessionInfo", "RTF", function(this, locale = TRUE, ...) {
	
	si<-sessionInfo()
	opkgver <- sapply(si$otherPkgs, function(x) x$Version)
    nspkgver <- sapply(si$loadedOnly, function(x) x$Version)
	
	startParagraph(this)
	addText(this,"\\s1 Session Information",bold=TRUE)
	endParagraph(this)
	
	startParagraph(this)
	addText(this,si$R.version$version.string,bold=TRUE,italic=TRUE)
	endParagraph(this)
	
	increaseIndent(this)
	startParagraph(this)
	addText(this,"Platform: ")
    addText(this,si$R.version$platform,italic=TRUE)
    
	if (locale) {
		addText(this,"\nLocale: ")
        addText(this,si$locale,italic=TRUE)
    }
    endParagraph(this)
    decreaseIndent(this)
    
    startParagraph(this)
    addText(this,"Packages",bold=TRUE,italic=TRUE)
    endParagraph(this)
    
    increaseIndent(this)
    startParagraph(this)
    addText(this,"Base: ")
	addText(this, paste(sort(si$basePkgs), collapse = ", "),italic=TRUE)
	
	if (length(opkgver)) {
        opkgver <- opkgver[sort(names(opkgver))]
        addText(this,"\nOther: ")
        vers<-paste("(v",opkgver,")",sep="")
        addText(this, paste(names(opkgver), vers, sep = " ", collapse = ", "),italic=TRUE)
    }
    
    if (length(nspkgver)) {
        nspkgver <- nspkgver[sort(names(nspkgver))]
    	addText(this,"\nLoaded (not attached): ")
    	vers<-paste("(v",nspkgver,")",sep="")
    	addText(this, paste(names(nspkgver), vers, sep = " ", collapse = ", "),italic=TRUE)
    }
    endParagraph(this)
    decreaseIndent(this)
    
    startParagraph(this)
    addText(this,"Session Details",bold=TRUE,italic=TRUE)
    endParagraph(this)
    increaseIndent(this)
    startParagraph(this)
    addText(this,"Working directory: ")
	addText(this, getwd(),italic=TRUE)
	addText(this,"\nOutput file: ")
	addText(this, this$.file,italic=TRUE)
	endParagraph(this)
    decreaseIndent(this)
})

######################################################################################

.start.rtf<-function(width=8.5,height=11,omi=c(1,1,1,1)) {
	paste("{\\rtf1\\ansi\n\\deff",.add.font.table(),.add.paper.size(width=width,height=height),"\n",.add.page.margins(omi),"\n",.add.page.numbers(),"\n",sep="")
}

.add.font.table<-function() {
	fonts<-character()
	fonts[1]<-"{\\f1\\fswiss\\fcharset0 Helvetica;}"
	fonts[2]<-"{\\f2\\ffroman\\charset0\\fprg2 Times New Roman;}"
	fonts[3]<-"{\\f3\\ffswiss\\charset0\\fprg2 Arial;}"
	fonts[4]<-"{\\f4\\fftech\\charset0\\fprg2 Symbol;}"
	fonts[5]<-"{\\f4\\ffroman\\charset0\\fprg2 Cambria;}"
	
	paste("{\\fonttbl",paste(fonts,collapse="\n"),"}",sep="\n")
}

.add.page.numbers<-function() {
	paste("\\titlepg\\headery720\\footery720","{\\footer {\\pard\\qc\\fi0\\li0 \\f2\\fs20 \\field{\\fldinst{page}} \\par}}",sep="\n")
}

.add.paper.size<-function(width=8.5,height=11) {
	paste("\\paperw",round(width*1440,0),"\\paperh",round(height*1440,0),"\\widowctrl\\ftnbj\\fet0\\sectd",if(width>height){"\\lndscpsxn"} else {""},"\\linex0",sep="")
}

.add.page.margins<-function(margins=c(1,1,1,1)) {
	paste("\\margl",round(margins[2]*1440,0),"\\margr",round(margins[4]*1440,0),"\\margt",margins[3]*1440,"\\margb",margins[1]*1440,sep="")
}

.add.header<-function(title,subtitle=NULL,indent=0,font.size=10,TOC.level=NULL) {
	if(is.null(subtitle)) {
		paste("{\\pard\\fi0\\li",indent,"\\f2\\fs",font.size*2,"\\b",.get.TOC.level(TOC.level)," ",.convert(title),"\\b0\\line\\par}\n",sep="")
	} else {
		paste("{\\pard\\fi0\\li",indent,"\\f2\\fs",font.size*2,"\\b",.get.TOC.level(TOC.level)," ",.convert(title),"\\b0\\par}\n{\\pard\\fi0\\f2\\fs",font.size*2," ",.convert(subtitle),"\\line\\par}\n",sep="")
	}
}

.get.TOC.level<-function(section.level) {
	ret<-""
	
	if(!is.null(section.level)) {
		ret<-paste("\\s",section.level,sep="")
	}
	
	ret
}

.start.paragraph<-function(indent=0,font.size=10) {
	paste("{\\pard\\fi0\\li",indent,"\\f2\\fs",font.size*2,"\n",sep="")
}

.add.text<-function(x) {
	paste(.convert(x),sep="")
}

.end.paragraph<-function() {
	paste("\\par}\n",sep="")
}

.end.rtf<-function() {
	paste("}",sep="")
}

.add.table.row<-function(col.data=c("c1","c2","c3"),col.widths=c(1.0,4.5,1.0),col.justify=NULL,font.size=10,last.row=FALSE,indent=0, border.top=FALSE, border.bottom=FALSE,space.before=NULL, space.after=NULL) {
	header<-paste("\\trowd\\trgaph100\\trleft",indent,sep="")  # trqc for centered
	
	if(length(col.data) != length(col.widths)) {
		stop(paste("The number of data columns (",length(col.data),") doesn't match the column widths (",length(col.widths),")!  Input data: ",col.data,sep=""))
	}
	
	justify<-vector()
	justify["L"]<-"\\ql"
	justify["R"]<-"\\qr"
	justify["C"]<-"\\qc"
	justify["J"]<-"\\qj"
	
	# Default: Left justify everything except numeric which are right justified
	justify.v<-rep(justify["L"],length(col.data))
	numeric.cols<-which(!is.na(suppressWarnings(as.numeric(col.data))))
	if(length(numeric.cols)>0) { justify.v[numeric.cols] <- justify["R"] }
	
	if(!is.null(col.justify)) {	
		if(length(col.justify)==1) {
			if(col.justify %in% names(justify)) {
				justify.v<-rep(justify[col.justify],length(col.data))
			} else {
				stop(paste("col.justify parameter not recognized: ",col.justify," (should be L, R, C, or J)",sep=""))
			}
		} else if(length(col.justify)==length(col.data)) {
			justify.v<-justify[col.justify]
		} else {
			stop(paste("The number of data columns (",length(col.data),") doesn't match the col.justify (",length(col.justify),") parameter!  Input data: ",paste(col.data,sep="",collapse=", "),sep=""))
		}
	}
	
	# TODO: allow decimal alignment: RTF code: \tqdec\tx810 (where the tx part is the tabbed position)
	
	btop<-""
	bbottom<-""
	
	if(border.top == TRUE) btop <- "\\clbrdrt\\brdrs\\brdrw15"
	if(last.row==TRUE | border.bottom==TRUE) bbottom <- "\\clbrdrb\\brdrs\\brdrw15"
	
	# Top vertical alignment for rows (clvertalt)
	cols.prefix<-paste("\\clvertalt\\clshdrawnil\\clwWidth",round(col.widths*1440,0),"\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph",btop,bbottom ,"\\cellx",c(1:length(col.widths)),"\n",sep="",collapse="")
	cols<-paste("\\pard",justify.v, .get.space.before.after(space.before, space.after),"\\widctlpar\\intbl\\fi0\\f2\\fs",font.size*2," ", .convert(col.data),"\\cell\n",sep="",collapse="")
	end.row<-"\\widctlpar\\intbl\\row\n"
	paste(header,cols.prefix,cols,end.row,sep="")
}

.get.space.before.after<-function(space.before=NULL,space.after=NULL) {
	ret<-""
	
	# \sbN -- N twips of extra (vertical) space before this paragraph (default: 0)
	# \saN -- N twips of extra (vertical) space after this paragraph (default: 0)
	
	if(!is.null(space.before)) {
		ret<-paste(ret,"\\sb",(space.before*1440),sep="")
	}
	
	if(!is.null(space.after)) {
		ret<-paste(ret,"\\sa",(space.after*1440),sep="")
	}
	
	ret
}

.add.merged.table.row<-function(col.data=c("c1","c2","c3"),col.widths=c(1.0,4.5,1.0),justify="LEFT",font.size=10,last.row=FALSE,indent=0, border.top=FALSE, border.bottom=FALSE) {
	header<-paste("\\trowd\\trgaph100\\trleft",indent,sep="")  # trqc for centered
	
	justify.q<-"\\ql"
	if(justify=="LEFT") justify.q<-"\\ql"
	if(justify=="RIGHT") justify.q<-"\\qr"
	if(justify=="CENTER") justify.q<-"\\qc"
	if(justify=="JUSTIFY") justify.q<-"\\qj"
	
	btop<-""
	bbottom<-""
	
	if(border.top == TRUE) btop <- "\\clbrdrt\\brdrs\\brdrw15"
	if(last.row==TRUE | border.bottom==TRUE) bbottom <- "\\clbrdrb\\brdrs\\brdrw15"
	
	merged<-c("","\\clmgf",rep("\\clmrg",length(col.data)-2))
	
	cols.prefix<-paste("\\clvertalc \\clshdrawnil \\clwWidth",round(col.widths*1440,0),"\\clftsWidth3 \\clheight260 \\clpadl100 \\clpadr100 \\gaph",btop," ",bbottom ,merged,"\\cellx",c(1:length(col.widths)),"\n",sep="",collapse="")
	cols<-paste("\\pard",justify.q,"\\widctlpar\\intbl\\fi0\\f2\\fs",font.size*2," ",.convert(col.data),"\\cell\n",sep="",collapse="")
	end.row<-"\\widctlpar\\intbl \\row \n\n"
	paste(header,cols.prefix,cols,end.row,sep="")
}

.add.table.header.row<-function(col.data=c("c1","c2","c3"),col.widths=c(1.0,4.5,1.0),col.justify=NULL,font.size=10,repeat.header=FALSE,indent=0) {
	header<-paste("\\trowd\\trgaph100\\trleft",indent,sep="")  # trqc for centered

		if(length(col.data) != length(col.widths)) {
		stop(paste("The number of data columns (",length(col.data),") doesn't match the column widths (",length(col.widths),")!  Input data: ",col.data,sep=""))
	}
	
	justify<-vector()
	justify["L"]<-"\\ql"
	justify["R"]<-"\\qr"
	justify["C"]<-"\\qc"
	justify["J"]<-"\\qj"
	
	# Default: Left justify everything
	justify.v<-rep(justify["L"],length(col.data))
	
	if(!is.null(col.justify)) {	
		if(length(col.justify)==1) {
			if(col.justify %in% names(justify)) {
				justify.v<-rep(justify[col.justify],length(col.data))
			} else {
				stop(paste("header.col.justify parameter not recognized: ",col.justify," (should be L, R, C, or J)",sep=""))
			}
		} else if(length(col.justify)==length(col.data)) {
			justify.v<-justify[col.justify]
		} else {
			stop(paste("The number of data columns (",length(col.data),") doesn't match the header.col.justify (",length(col.justify),") parameter!  Input data: ",paste(col.data,sep="",collapse=", "),sep=""))
		}
	}

	if(repeat.header==TRUE) header<-paste(header,"\\trhdr")
	
	# Bottom vertical alignment for headers (clvertalb)
	cols.prefix<-paste("\\clvertalb\\clshdrawnil\\clwWidth",round(col.widths*1440,0),"\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\clbrdrb\\brdrs\\brdrw15\\cellx",c(1:length(col.widths)),"\n",sep="",collapse="")
	cols<-paste("\\pard",justify.v,"\\widctlpar\\intbl\\fi0\\f2\\fs",font.size*2,"\\b ",.convert(col.data),"\\b0\\cell\n",sep="",collapse="")
	end.row<-"\\widctlpar\\intbl\\row\n\n"
	paste(header,cols.prefix,cols,end.row,sep="")
}


.add.table<-function(dat,col.widths=NULL,col.justify=NULL,header.col.justify=NULL,font.size=10,row.names=FALSE,indent=0,NA.string="-",max.table.width=NULL, space.before=NULL, space.after=NULL) {
	ret<-"{\\pard\n"
	
	if("table" %in% class(dat)) {
		if(length(dim(dat))==1) {
		
			varnames<-names(dimnames(dat))[1]
			nc<-2
			nr<-length(dimnames(dat)[[1]])
			
			if(is.null(col.widths)){ col.widths<-rep(6.5/nc,nc)}
			
			ret<-paste(ret,.add.table.header.row(c(names(dimnames(dat))[1]," "),col.widths,header.col.justify,font.size=font.size,repeat.header=TRUE,indent=indent),sep="")
			
			if(nrow(dat)>1) {
				for(i in 1:(nrow(dat)-1) ) {
					rn<-rownames(dat)[i]
					ret<-paste(ret,.add.table.row(c(rn,as.character(dat[i])),col.widths,col.justify,font.size=font.size,indent=indent, space.before=space.before, space.after=space.after),sep="")
				}
			}
			
			rn<-rownames(dat)[nrow(dat)]
			ret<-paste(ret,.add.table.row(c(rn,as.character(dat[nrow(dat)])),col.widths,col.justify,font.size=font.size,indent=indent,border.bottom=TRUE, space.before=space.before, space.after=space.after),sep="")
		
		} else if(length(dim(dat))==2) {
			varnames<-names(dimnames(dat))
			nc<-ncol(dat)+1
			nr<-nrow(dat)
			
			if(is.null(col.widths)){ col.widths<-rep(6.5/nc,nc)}
			
			# ret<-paste(ret,.add.table.header.row(c(" ",colnames(dat)),col.widths,font.size=font.size,repeat.header=TRUE,indent=indent),sep="")
			
			ret<-paste(ret,.add.merged.table.row(c(" ",paste("\\b ",varnames[2]," \\b0",sep=""),rep(" ",nc-2)),col.widths,font.size=font.size,indent=indent,border.top=TRUE),sep="")
			ret<-paste(ret,.add.table.row(c(paste("\\b ",varnames[1]," \\b0",sep=""),colnames(dat)),col.widths,col.justify,font.size=font.size,indent=indent,border.bottom=TRUE),sep="")
			
			if(nrow(dat)>1) {
				for(i in 1:(nrow(dat)-1) ) {
					rn<-rownames(dat)[i]
					ret<-paste(ret,.add.table.row(c(rn,as.character(dat[i,])),col.widths,col.justify,font.size=font.size,indent=indent, space.before=space.before, space.after=space.after),sep="")
				}
			}
			
			rn<-rownames(dat)[nrow(dat)]
			ret<-paste(ret,.add.table.row(c(rn,as.character(dat[nrow(dat),])),col.widths,col.justify,font.size=font.size,indent=indent,border.bottom=TRUE, space.before=space.before, space.after=space.after),sep="")
		
		} else {
			stop("Table dimensions can't be written")
		}
	
	} else if("xtab" %in% class(dat)) {
		
		nc<-ncol(dat$counts)+2
		nr<-nrow(dat$counts)
		
		if(is.null(col.widths)){ col.widths<-rep(6.5/nc,nc)}
		
		# ret<-paste(ret,.add.table.header.row(c(" ",colnames(dat)),col.widths,font.size=font.size,repeat.header=TRUE,indent=indent),sep="")
		
		ret<-paste(ret,.add.merged.table.row(c(" ",paste("\\b ",dat$varnames[2]," \\b0",sep=""),rep(" ",nc-2)),col.widths,font.size=font.size,indent=indent,border.top=TRUE),sep="")
		ret<-paste(ret,.add.table.row(c(paste("\\b ",dat$varnames[1]," \\b0",sep=""),colnames(dat$counts),"Total"),col.widths,col.justify,font.size=font.size,indent=indent,border.bottom=TRUE),sep="")
		grand.total<-sum(dat$col.margin)
		
		if(nrow(dat$counts)>1) {
			for(i in 1:(nrow(dat$counts)) ) {
				rn<-rownames(dat$counts)[i]
				ret<-paste(ret,.add.table.row(c(rn,as.character(dat$counts[i,]),paste( dat$row.margin[i]," (",sprintf("%0.1f",dat$row.margin[i]/grand.total*100),"%)" ,sep="")),col.widths,col.justify,font.size=font.size,indent=indent, space.before=space.before, space.after=space.after),sep="")
				ret<-paste(ret,.add.table.row(c(" ",paste("(",sprintf("%0.1f",dat$counts[i,]/dat$row.margin[i]*100),"% R)",sep="")," "),col.widths,col.justify,font.size=font.size,indent=indent, space.before=space.before, space.after=space.after),sep="")
				ret<-paste(ret,.add.table.row(c(" ",paste("(",sprintf("%0.1f",dat$counts[i,]/dat$col.margin*100),"% C)",sep="")," "),col.widths,col.justify,font.size=font.size,indent=indent, space.before=space.before, space.after=space.after),sep="")
				ret<-paste(ret,.add.table.row(rep(" ",nc),col.widths,col.justify,font.size=font.size,indent=indent, space.before=space.before, space.after=space.after),sep="")
			}
		}
		
		# Total rows
		ret<-paste(ret,.add.table.row(c("Total",paste(as.character(dat$col.margin),paste(" (",sprintf("%0.1f",dat$col.margin/grand.total*100),"%)",sep="")),as.character(grand.total)),col.widths,font.size=font.size,last.row=TRUE,indent=indent, space.before=space.before, space.after=space.after),sep="")
		
	# handle etables (etables are just matrices with a 'start cell' defined)
	} else if ("matrix" %in% class(dat) & !is.null(attributes(dat)$"start cell")) { 
		start.row<-attributes(dat)$"start cell"[1]
		
		# convert matrix to data frame
		dat<-as.data.frame(dat,stringsAsFactors=FALSE)
		dat[is.na(dat)] <- NA.string
		dat[dat=="NA"] <- NA.string
		
		# if no column widths are specified, then calculate optimal sizes that fit the page
		if(is.null(col.widths) & !is.null(max.table.width)) {
			col.widths <- .optimize.col.widths(dat,include.row.names=row.names,max.table.width=max.table.width,font.size=font.size)
		}
		
		# render the header
		nc<-ncol(dat)
		if(is.null(col.widths)){ col.widths<-rep(6.5/nc,nc)}
		
		if(nrow(dat)>1) {
			for(i in 1:(nrow(dat)-1)) {
				if(i<start.row) {
					border.top=FALSE
					border.bottom=FALSE
					if (i==1) { border.top=TRUE }
					if (i==(start.row-1)) { border.bottom=TRUE }
					ret<-paste(ret,.add.table.row(paste("\\b ",as.character(dat[i,])," \\b0",sep=""),col.widths,col.justify,font.size=font.size,indent=indent,border.top=border.top,border.bottom=border.bottom,space.before=space.before, space.after=space.after),sep="")

				} else {
					ret<-paste(ret,.add.table.row(as.character(dat[i,]),col.widths,col.justify,font.size=font.size,indent=indent, space.before=space.before, space.after=space.after),sep="")
				}
			}
			ret<-paste(ret,.add.table.row(as.character(dat[nrow(dat),]),col.widths,col.justify,font.size=font.size,last.row=TRUE,indent=indent, space.before=space.before, space.after=space.after),sep="")
		}
	
	} else if ("data.frame" %in% class(dat) || "matrix" %in% class(dat)) {
		# convert matrix to data frame
		if("matrix" %in% class(dat)) {
			dat<-as.data.frame(dat)
		}
		
		# convert factor columns in a data frame to characters
		rnames<-rownames(dat)
		is.na(dat) <- is.na(dat) # handle NaN values by converting to NAs since this throws an error: dat[is.nan(dat)] <- NA.string
		dat<-data.frame(lapply(dat,as.character),stringsAsFactors=FALSE,check.names=FALSE)
		dat[is.na(dat)] <- NA.string
		dat[dat=="NA"] <- NA.string
		rownames(dat)<-rnames
		
		# if no column widths are specified, then calculate optimal sizes that fit the page
		if(is.null(col.widths) & !is.null(max.table.width)) {
			col.widths <- .optimize.col.widths(dat,include.row.names=row.names,max.table.width=max.table.width,font.size=font.size)
		}
		
		# render the header
		nc<-ncol(dat)
		if(row.names==TRUE){ nc<-nc+1 }
		
		if(is.null(col.widths)){ col.widths<-rep(6.5/nc,nc)}
		
		
		if(row.names==TRUE){
			ret<-paste(ret,.add.table.header.row(c(" ",colnames(dat)),col.widths,header.col.justify,font.size=font.size,repeat.header=TRUE,indent=indent),sep="")
		} else {
			ret<-paste(ret,.add.table.header.row(colnames(dat),col.widths,header.col.justify,font.size=font.size,repeat.header=TRUE,indent=indent),sep="")
		}
		
		# render the rows
		if(nrow(dat)>1) {
			for(i in 1:(nrow(dat)-1)) {
				if(row.names==TRUE){
					rn<-rownames(dat)[i]
					ret<-paste(ret,.add.table.row(c(rn,as.character(dat[i,])),col.widths,col.justify,font.size=font.size,indent=indent, space.before=space.before, space.after=space.after),sep="")
				} else {
					ret<-paste(ret,.add.table.row(as.character(dat[i,]),col.widths,col.justify,font.size=font.size,indent=indent, space.before=space.before, space.after=space.after),sep="")
				}
			}
		}
		
		if(row.names==TRUE){
			rn<-rownames(dat)[nrow(dat)]
			ret<-paste(ret,.add.table.row(c(rn,as.character(dat[nrow(dat),])),col.widths,col.justify,font.size=font.size,last.row=TRUE,indent=indent, space.before=space.before, space.after=space.after),sep="")
		} else {
			ret<-paste(ret,.add.table.row(as.character(dat[nrow(dat),]),col.widths,col.justify,font.size=font.size,last.row=TRUE,indent=indent, space.before=space.before, space.after=space.after),sep="")
		}
		
	} else {
		warning("No suitable RTF converter for object class!")
	}
	
	ret<-paste(ret,"}\n\n",sep="")
	ret
}

.add.page.break<-function(width=8.5,height=11,omi=c(1,1,1,1)) {
	#	"\\pard {\\f1 \\sect } \\sectd \\lndscpsxn\\pgwsxn16840\\pghsxn11907\\left\\widctlpar\\fi0\\f2\\fs18 \\par"
	# previous: "\\pard {\\f1 \\column }\\left\\widctlpar\\fi0\\f2\\fs18 \\par"
	paste("\\pard{\\f1\\sect}\\sectd",.add.paper.size(width=width,height=height),.add.page.margins(omi),"\\left\\widctlpar\\fi0\\f2\\fs18",sep="")
}

.convert<-function(x) {
	# http://www.ssec.wisc.edu/~tomw/java/unicode.html
	#x<-gsubfn("\\u(\\d+)", .hex2dec, x, engine="R")     # format UTF-8 characters from hex to dec
	#x<-gsub("\\u(\\d+)","\\\\u\\1\\\\3",x)  # format UTF-8 characters from hex to dec
	
	x<-gsub("\\n"," \\\\line ",x)         # convert new line to RTF \line
	#x<-gsub("\\t"," \\\\tab ",x)         # convert tab to RTF \tab
	x<-gsub("<=","\\\\u8804\\\\3",x)      # convert <= to RTF symbol
	x<-gsub(">=","\\\\u8805\\\\3",x)      # convert >= to RTF symbol
	
# 	x<-gsub(":delta:","\\\\u0916\\\\3",x) # convert :delta: to uppercase Greek delta
# 	
# 	x<-gsub("&alpha;","\\\\u0945\\\\3",x) # convert &alpha; to lowercase Greek alpha
# 	x<-gsub("&beta;","\\\\u0946\\\\3",x)  # convert &beta; to lowercase Greek beta
# 	x<-gsub("&gamma;","\\\\u0947\\\\3",x) # convert &gamma; to lowercase Greek gamma
# 	x<-gsub("&delta;","\\\\u0948\\\\3",x) # convert &delta; to lowercase Greek delta
# 	x<-gsub("&epsilon;","\\\\u0949\\\\3",x) # convert &epsilon; to lowercase Greek epsilon
# 	x<-gsub("&theta;","\\\\u0952\\\\3",x) # convert &theta; to lowercase Greek theta
# 	x<-gsub("&kappa;","\\\\u0954\\\\3",x) # convert &kappa; to lowercase Greek kappa
# 	x<-gsub("&lambda;","\\\\u0955\\\\3",x) # convert &lambda; to lowercase Greek lambda
# 	x<-gsub("&mu;","\\\\u0956\\\\3",x)    # convert &mu; to lowercase Greek lambda
# 	
# 	x<-gsub("&Alpha;","\\\\u0913\\\\3",x) # convert &Alpha; to uppercase Greek alpha
# 	x<-gsub("&Beta;","\\\\u0914\\\\3",x)  # convert &Beta; to uppercase Greek beta
# 	x<-gsub("&Gamma;","\\\\u0915\\\\3",x) # convert &Gamma; to uppercase Greek gamma
# 	x<-gsub("&Delta;","\\\\u0916\\\\3",x) # convert &Delta; to uppercase Greek delta
# 	x<-gsub("&Epsilon;","\\\\u0917\\\\3",x) # convert &Epsilon; to uppercase Greek epsilon
# 	x<-gsub("&Theta;","\\\\u0920\\\\3",x) # convert &Theta; to uppercase Greek theta
# 	x<-gsub("&Kappa;","\\\\u0922\\\\3",x) # convert &Kappa; to lowercase Greek kappa
# 	x<-gsub("&Lambda;","\\\\u0923\\\\3",x) # convert &Lambda; to lowercase Greek lambda
# 	x<-gsub("&Mu;","\\\\u0924\\\\3",x)     # convert &Mu; to lowercase Greek lambda
	
	# convert HTML characters
	x<-gsub("&gt;",">",x)
	x<-gsub("&lt;","<",x)
	
	# convert uppercase and lowercase Greek letters
	x<-gsub("&Alpha;","\\\\u0913\\\\3",x)
	x<-gsub("&Beta;","\\\\u0914\\\\3",x)
	x<-gsub("&Gamma;","\\\\u0915\\\\3",x)
	x<-gsub("&Delta;","\\\\u0916\\\\3",x)
	x<-gsub("&Epsilon;","\\\\u0917\\\\3",x)
	x<-gsub("&Zeta;","\\\\u0918\\\\3",x)
	x<-gsub("&Eta;","\\\\u0919\\\\3",x)
	x<-gsub("&Theta;","\\\\u0920\\\\3",x)
	x<-gsub("&Iota;","\\\\u0921\\\\3",x)
	x<-gsub("&Kappa;","\\\\u0922\\\\3",x)
	x<-gsub("&Lambda;","\\\\u0923\\\\3",x)
	x<-gsub("&Mu;","\\\\u0924\\\\3",x)
	x<-gsub("&Nu;","\\\\u0925\\\\3",x)
	x<-gsub("&Xi;","\\\\u0926\\\\3",x)
	x<-gsub("&Omicron;","\\\\u0927\\\\3",x)
	x<-gsub("&Pi;","\\\\u0928\\\\3",x)
	x<-gsub("&Rho;","\\\\u0929\\\\3",x)
	x<-gsub("&Sigma;","\\\\u0931\\\\3",x)
	x<-gsub("&Tau;","\\\\u0932\\\\3",x)
	x<-gsub("&Upsilon;","\\\\u0933\\\\3",x)
	x<-gsub("&Phi;","\\\\u0934\\\\3",x)
	x<-gsub("&Chi;","\\\\u0935\\\\3",x)
	x<-gsub("&Psi;","\\\\u0936\\\\3",x)
	x<-gsub("&Omega;","\\\\u0937\\\\3",x)
	x<-gsub("&alpha;","\\\\u0945\\\\3",x)
	x<-gsub("&beta;","\\\\u0946\\\\3",x)
	x<-gsub("&gamma;","\\\\u0947\\\\3",x)
	x<-gsub("&delta;","\\\\u0948\\\\3",x)
	x<-gsub("&epsilon;","\\\\u0949\\\\3",x)
	x<-gsub("&zeta;","\\\\u0950\\\\3",x)
	x<-gsub("&eta;","\\\\u0951\\\\3",x)
	x<-gsub("&theta;","\\\\u0952\\\\3",x)
	x<-gsub("&iota;","\\\\u0953\\\\3",x)
	x<-gsub("&kappa;","\\\\u0954\\\\3",x)
	x<-gsub("&lambda;","\\\\u0955\\\\3",x)
	x<-gsub("&mu;","\\\\u0956\\\\3",x)
	x<-gsub("&nu;","\\\\u0957\\\\3",x)
	x<-gsub("&xi;","\\\\u0958\\\\3",x)
	x<-gsub("&omicron;","\\\\u0959\\\\3",x)
	x<-gsub("&pi;","\\\\u0960\\\\3",x)
	x<-gsub("&rho;","\\\\u0961\\\\3",x)
	x<-gsub("&sigmaf;","\\\\u0962\\\\3",x)
	x<-gsub("&sigma;","\\\\u0963\\\\3",x)
	x<-gsub("&tau;","\\\\u0964\\\\3",x)
	x<-gsub("&upsilon;","\\\\u0965\\\\3",x)
	x<-gsub("&phi;","\\\\u0966\\\\3",x)
	x<-gsub("&chi;","\\\\u0967\\\\3",x)
	x<-gsub("&psi;","\\\\u0968\\\\3",x)
	x<-gsub("&omega;","\\\\u0969\\\\3",x)
	
	x<-gsub("TRUE","Yes",x)
	x<-gsub("FALSE","No",x)
	x
}

.add.newline<-function(n=NULL, font.size=10) {
	# return("\\line ")

	ret<-paste("{\\pard\\fi0\\f2\\fs",(font.size*2),sep="")
	
	if(is.null(n)) {
		if(n>=2) {
			ret<-paste(ret,paste(rep("\\line",n),"\n",collapse="",sep=""),sep="")
		}
	}
	paste(ret,"\\par}",sep="")
}

.chunk.vector<-function(tokens,n=10) {
 	nlines<-as.integer(length(tokens)/n)+1
	ntokens.line <- n #ceiling(length(tokens) / nlines) # tokens per line
	token.list <- split(tokens, rep( 1:ntokens.line, each=ntokens.line, len=length(tokens)))
	munged<-lapply(token.list,paste,collapse="")
	do.call(paste,list(munged,collapse="\n"))
}

# width and height are in inches
.add.png<-function(file,width=3,height=3,verbose=FALSE) {
	# return a hexadecimal version of a file
	max.bytes<-50000000  # maximum file size in bytes (~50MB)
	dat<-readBin(file, what="raw", size=1, signed=TRUE, endian="little",n=max.bytes);
	if(verbose) {
		cat(paste(length(dat),"bytes read\n"))
	}
	paste("{\\rtf1\\ansi\\deff0{\\pict\\pngblip\\picwgoal",round(width*1440),"\\pichgoal",round(height*1440)," ",paste(dat,collapse=""),"}}",sep="")
	# paste("{\\rtf1\\ansi\\deff0{\\pict\\pngblip\\picwgoal",round(width*1440),"\\pichgoal",round(height*1440)," \n",.chunk.vector(dat),"}}",sep="")
}

.rtf.plot<-function(plot.fun,tmp.file="temp.png",width=3.0,height=0.3,res=300, ...) {
	width.px<-round(width*res)
	height.px<-round(height*res)
	#png(tmp.file,width=width.px,height=height.px,units="px",pointsize=8,bg = "white",res=res)
	png(tmp.file,width=width.px,height=height.px,units="px",pointsize=8,bg = "transparent",res=res)
	plot.fun(...)
	dev.off()
	.add.png(tmp.file,width=width,height=height)
}

.rtf.trellis.object<-function(trellis.object,tmp.file="temp.png",width=3.0,height=0.3,res=300,rotate=NULL, ...) {
	if(class(trellis.object) != "trellis") {
		stop("Not a trellis object!")
	}
	
	ret<-""
	
# 	if(is.null(trellis.object$layout)) {
# 		# single page
# 		width.px<-round(width*res)
# 		height.px<-round(height*res)
# 		png(tmp.file,width=width.px,height=height.px,units="px",pointsize=8,bg = "white",res=res)
# 		print(trellis.object)
# 		dev.off()
# 		
# 		if(!is.null(rotate)) {
# 			system(paste("convert -rotate ",rotate," '",tmp.file,"' '",tmp.file,"'",sep=""))
# 			ret<-.add.png(tmp.file,width=height,height=width) # swap width and height
# 		} else {
# 			ret<-.add.png(tmp.file,width=width,height=height)
# 		}
# 		
# 	} else {
# 		plot.cnt<-dim(trellis.object)
# 		per.page<-trellis.object$layout[1]*trellis.object$layout[2]
# 		pages<-floor((plot.cnt-1)/per.page)+1
# 		
# 		for(pg in 1:pages) {
# 			plot.start<-(pg-1)*per.page+1
# 			plot.end<-min(plot.cnt,pg*per.page)
# 			
# 			width.px<-round(width*res)
# 			height.px<-round(height*res)
# 			png(tmp.file,width=width.px,height=height.px,units="px",pointsize=8,bg = "white",res=res)
# 			print(trellis.object[plot.start:plot.end])
# 			dev.off()
# 			
# 			if(!is.null(rotate)) {
# 				system(paste("convert -rotate ",rotate," '",tmp.file,"' '",tmp.file,"'",sep=""))
# 				ret<-paste(ret,.add.png(tmp.file,width=height,height=width),sep="\n") # swap width and height
# 			} else {
# 				ret<-paste(ret,.add.png(tmp.file,width=width,height=height),sep="\n")
# 			}
# 		}
# 	}
	
	if(is.null(trellis.object$layout)) {
		# single page
		width.px<-round(width*res)
		height.px<-round(height*res)
		
		if(!is.null(rotate)) {
			
			png(tmp.file,width=height.px,height=width.px,units="px",pointsize=8,bg = "transparent",res=res)
			#grid::grid.newpage()
			grid::pushViewport(grid::viewport(width=grid::unit(width,"inches"),height=grid::unit(height,"inches"),angle = rotate))
			print(trellis.object, newpage = FALSE)
			#grid::upViewport()
			dev.off()
			
			ret<-.add.png(tmp.file,width=height,height=width) # swap width and height
			
		} else {
			png(tmp.file,width=width.px,height=height.px,units="px",pointsize=8,bg = "transparent",res=res)
			print(trellis.object)
			dev.off()
			
			ret<-.add.png(tmp.file,width=width,height=height)
		}
		
	} else {
		plot.cnt<-dim(trellis.object)
		per.page<-trellis.object$layout[1]*trellis.object$layout[2]
		pages<-floor((plot.cnt-1)/per.page)+1
		
		for(pg in 1:pages) {
			plot.start<-(pg-1)*per.page+1
			plot.end<-min(plot.cnt,pg*per.page)
			
			width.px<-round(width*res)
			height.px<-round(height*res)
			
			if(!is.null(rotate)) {
				png(tmp.file,width=height.px,height=width.px,units="px",pointsize=8,bg = "transparent",res=res)
				#grid::grid.newpage()
				grid::pushViewport(grid::viewport(width=grid::unit(width,"inches"),height=grid::unit(height,"inches"),angle = rotate))
				print(trellis.object[plot.start:plot.end], newpage = FALSE)
				#grid::upViewport()
				dev.off()
				ret<-paste(ret,.add.png(tmp.file,width=height,height=width),sep="\n") # swap width and height
			} else {
				png(tmp.file,width=width.px,height=height.px,units="px",pointsize=8,bg = "transparent",res=res)
				print(trellis.object[plot.start:plot.end])
				dev.off()
				
				ret<-paste(ret,.add.png(tmp.file,width=width,height=height),sep="\n")
			}
		}
	}
	
	ret
}

.max.col.nchar<-function(x,include.row.names=FALSE,wrap.headers=TRUE) {
	# for each column, returns the maximum width (in characters) of either the 
	# largest header word or largest entire column field.  This allows 
	# headers to be longer and have word wrapping while attempting
	# to keep each field of table data on a single line.
	
	contents<-apply(x,2,function(x){max(nchar(x))} )
	col.nchar<-c()
	
	if(wrap.headers==TRUE) {
		# find maximum word width in each header (split on spaces)
		headers<-sapply(names(x), function(x){ max(nchar(strsplit(x," ")[[1]])) } )
		col.nchar<-mapply(max,contents,headers)
	} else {
		headers<-nchar(names(x))
		col.nchar<-mapply(max,contents,headers)
	}
	
	if(include.row.names==TRUE) {
		row.names.nchar<-max(sapply(rownames(x),nchar))
		col.nchar<-c(row.names.nchar,col.nchar)
	}
	
	col.nchar
}

# TODO this code is not so great.  Can we take advantage of R's strwidth function
# and plug in the font faces and size in points to estimate the width. E.g.
# 
# strwidth("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.?,;:'\"!@#$%^&*()-=+_[]{}|\\",units="inches",family="Times",ps=12)
# 
# where face = (1=plain, 2=bold, 3=italic, 4=bold-italic)
# see also : ps.options()

.optimize.col.widths<-function(x,include.row.names=FALSE,max.table.width=6.5,font.size=9,col.padding=0.1) {
	letter.width<-font.size * 1/144    # font point size to width (roughly 1/144 inch)
	letter.width<-letter.width + 0.03  # fine tuning to account for all caps columns or bold
	
	max.nchars <- .max.col.nchar(x,include.row.names,wrap.headers=TRUE)
	col.widths <- max.nchars*letter.width+2*col.padding
	
	# This could also be tweaked with a more formal character width analysis:
	# http://stephensite.net/WordPressSS/2008/02/19/how-to-calculate-the-character-width-accross-fonts-and-points/

	# If table is still too wide, resize each column proportionally to the content to
	# fit the maximum table width allowed
	if(sum(col.widths) > max.table.width) {
		col.widths<-col.widths/sum(col.widths) * max.table.width
	}
	
	col.widths
}

.hex2dec <- function(hexadecimal) {
	hexdigits <- c(0:9, LETTERS[1:6])
	hexadecimal <- toupper(hexadecimal)    # treat upper/lower case the same
	decimal <- rep(0, length(hexadecimal))
	for (i in 1:length(hexadecimal)) {
		digits <- match(strsplit(hexadecimal[i],"")[[1]], hexdigits) - 1
		decimal[i] <- sum(digits * 16^((length(digits)-1):0))
	}
	return(decimal)
}