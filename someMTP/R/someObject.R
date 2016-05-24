#==========================================================
# CLASS DEFINITION *** CLASS DEFINITION *** CLASS DEFINITION
#==========================================================

setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("vectorOrNULL", c("vector", "NULL"))
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("callOrNULL", c("call", "NULL"))


setClass("someMTP.object", 
  representation(
	call = "callOrNULL",
    rej = "vectorOrNULL", 
	p = "vectorOrNULL",
	ord = "vectorOrNULL", 
	idOrd = "vectorOrNULL", 
    MTP = "character",
	GD = "logical",
    q = "numericOrNULL",
	k = "numericOrNULL",
	J = "numericOrNULL",
    alpha = "numericOrNULL",
	alphaprime = "numericOrNULL"
  ),
  prototype = list(
	call = NULL,
    rej = NULL,
	p = NULL,
	ord = NULL,
	idOrd = NULL,
    MTP = NULL,
	GD = FALSE,
    q = NULL,
	k = NULL,
	J = NULL,
    alpha = NULL,
	alphaprime = NULL
  )
)

################################################
setClass("lsd.object", 
  representation(
	globalP = "vectorOrNULL", 
    MTP = "character",
	F = "numericOrNULL",
	df = "vectorOrNULL",
	D = "matrixOrNULL",
	call = "callOrNULL"
  ),
  prototype = list(
	globalP = NULL, 
    MTP = "lsd",
    F = NULL,
	df = NULL,
	D = NULL,
	call = NULL
  )
)



#==========================================================
# PUBLIC METHODS *** PUBLIC METHODS *** PUBLIC METHODS
#==========================================================

#==========================================================
# Function "show" prints a "gt.object" object
#==========================================================

setMethod("show", "someMTP.object", function(object)
{
  cat("someMTP result:\n")
  switch(object@MTP, 
  fdrOrd= cat(paste("Ordered FDR procedure ", ifelse(object@GD," for General Dependence", "" ),"\n ",
        length(object@rej)," tests, q=",round(object@q,digits=5),", individual alpha threshold=",round(object@alphaprime,digits=5),"\n ",sum(object@rej)," rejections\n\n",sep="")),
  kfweOrd= cat(paste("Ordered k-FWER procedure ", ifelse(object@GD," for General Dependence", "" ),"\n ",
        length(object@rej)," tests, alpha=",round(object@alpha,digits=5),", individual alpha threshold=",round(object@alphaprime,digits=5),",\n allowed Jumps=",object@J,"\n "
		,sum(object@rej)," rejections\n\n",sep="")),		
  object )
   cat("\n")
})


setGeneric("summary")
setMethod("summary", "someMTP.object", function(object, ...)
{
  cat("someMTP result:\n")
  switch(object@MTP, 
  fdrOrd = cat(paste("Ordered FDR procedure ", ifelse(object@GD," for General Dependence", "" ),"\n ",
        length(object@rej)," tests, q=",round(object@q,digits=5),", individual alpha threshold=",round(object@alphaprime,digits=5),"\n ",sum(object@rej)," rejections\n\n",sep="")),
  kfweOrd= cat(paste("Ordered k-FWER procedure ", ifelse(object@GD," for General Dependence", "" ),"\n ",
        length(object@rej)," tests, alpha=",round(object@alpha,digits=5),", individual alpha threshold=",round(object@alphaprime,digits=5),",\n allowed Jumps=",object@J,"\n "
		,sum(object@rej)," rejections\n\n",sep="")),		
  none = "method = none")
   cat("\n")
})



#==========================================================
# A sort method for "gt.object"
#==========================================================
setMethod("sort", matchSignature(signature(x = "someMTP.object"), sort),
  function(x) {
  x @ rej = x @ rej[x @idOrd]
  x @ ord = x @ ord[x @idOrd]
  x @ p   = x @ p[x @idOrd]
  x @ idOrd = x @ idOrd[x @idOrd]
  return(x)
  }
)

#==========================================================
# The length method for "gt.object"
#==========================================================
setMethod("length", "someMTP.object", 
            function(x) {
  length(x@rej)
})            


#==========================================================
# The names and alias methods for "gt.object" 
# (applies to pathwaynames)
#==========================================================
setMethod("names", "someMTP.object", 
            function(x) 
{
  names(x@rej)
})      


setMethod("names<-", "someMTP.object", 
            function(x, value) 
{
  names(x@rej) <- value
  x
})            




#==========================================================
# Graph plot for focus level and inheritance procedures
#==========================================================
draw <- function(object, what = c("all","ordVsP", "stepVsR"), pdfName = NULL) {

  
  # find type if missing
  if (missing(what)) 
    what <- "all"
  else
    what <- match.arg(what,c("all","ordVsP", "stepVsR"))  
  
  
  # make the graph object
  if (!is.null(pdfName)) pdf(pdfName,width=ifelse((what == "all")& (object@MTP=="fdrOrd"), 20,10))
  
  if ((what == "all")& (object@MTP=="fdrOrd"))  par(mfrow=c(1,2))
  if (what %in% c("all", "ordVsP")) {
  par(cex=1.5)
	plot(object@ord,-log10(object@p),xlab="Ordering values", ylab="-log10(p-values)",pch=20,axes=TRUE,main="Ordering Criterion vs -log10(p-values)",col=object@rej+1)
	abline(-log10(ifelse(object@MTP=="fdrOrd", object@q,object@alpha)),0,col="gray",lwd=2)
	vline= switch(object@MTP,
	       fdrOrd=(length(object@rej)-which.max(which(cumsum(object@rej[object@idOrd[length(object@idOrd):1]])==0))+1)[1], 
		   kfweOrd = which(cumsum((1-object@rej)[object@idOrd])==(object@J+1))[1] ,NULL) 
	abline(v=ifelse(!is.na(vline),object@ord[object@idOrd[vline]],0) ,col="gray",lwd=2)
	# axis(1)
	# axis(2)
  }
  if ((what %in% c("all", "stepVsR")) & (object@MTP=="fdrOrd")) {
   par(cex=1.5)
	plot(cumsum(object@p[object@idOrd]<= object@alphaprime),xlab="Steps", ylab=paste("# of rejections (q=",object@q,")",sep=""),lwd=2,axes=TRUE,main="Step vs Number of Rejections",type="l")
    legend(x=.05, y=sum(object@p<=object@alphaprime),legend=c("Maximum","Stop if below"),lty=1,col=c("gray","red"),lwd=2,bty="n")
	abline(0,1/(2-object@q),col="red",lwd=2)
	abline(0,1,col="gray",lwd=2)
	par(cex=1)
	
	# axis(1)
	# axis(2)
  }
  if (!is.null(pdfName)) dev.off()
  
}


#############################################################################################################################################################
############ lsd.object

setMethod("show", "lsd.object", 
function(object)
{
  cat("someMTP result:\n")
  cat("Left Spherically Distributed - Test \n Call : ")
  #" ",dim(resp)[2]," variables Y, \n" ,dim(alternative)[2]," dependent variables X, \n  ",ifelse( ((dim(null)[2]>0) | is.null(null)),dim(null)[2],"NO"), " covariates Z \n
  cat(deparse(object@call), "\n df = ",paste(object@df, collapse=", "),"\n F = ",round(object@F,digits=5),"\n p-value = ",round(object@globalP,digits=5),"\n ") 
})

setGeneric("summary")
setMethod("summary", "lsd.object", function(object, showD=TRUE, ...)
{
  cat("someMTP result:\n")
  cat(deparse(object@call))
  cat("\n df = ",paste(object@df, collapse=", "),
       "\n F = ",round(object@F,digits=5),
	   "\n p-value = ",round(object@globalP,digits=5),"\n")
  if(showD) {cat("\n Variables Contribution:\n");
          print(object@D, digits = 3)
		  }
})



# #==========================================================
setGeneric("p.value", function(object, ...) standardGeneric("p.value"))
setMethod("p.value", "lsd.object",
  function(object) {
    object@globalP
  }
)
# #==========================================================
setGeneric("weights", function(object, ...) standardGeneric("weights"))
setMethod("weights", "lsd.object",
  function(object) {
    object@D
  }
)
