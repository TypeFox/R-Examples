#==========================================================
# CLASS DEFINITION *** CLASS DEFINITION *** CLASS DEFINITION
#==========================================================

#setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("numericOrmatrixOrNULL", c("numeric","matrix", "NULL"))
setClassUnion("arrayOrNULL", c("array", "NULL"))
setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))
setClassUnion("numericOrmatrixOrcharacterOrNULL", c("numeric","matrix", "NULL","character"))
setClassUnion("envOrNULL", c("environment", "NULL"))

# #############da togliere per compilazione (esistno gia in someMTP)
#  setClassUnion("numericOrNULL", c("numeric", "NULL"))
#  setClassUnion("listOrNULL", c("list", "NULL"))
#  require(e1071)


setClass("flip.object", 
  representation(
    res = "data.frameOrNULL",
    call = "call", 
	permP="arrayOrNULL",
	permT="arrayOrNULL",
	permSpace="listOrNULL",
	permY="arrayOrNULL",
    #functions = "environment",#"list",
    #subsets = "listOrNULL",
    #structure = "listOrNULL",
    #weights = "listOrNULL",
    tail = "numericOrmatrixOrcharacterOrNULL",
    #Z = "matrixOrNULL",
    #directional = "logical",
    data = "listOrNULL",
    call.env = "envOrNULL"
    #model = "character"
  ),
  prototype = list(
  res = NULL,
	permP=NULL,
	permT=NULL,
	permSpace=NULL,
	permY=NULL,
	data=NULL,
  call.env=NULL
  )
)

#==========================================================
# Function "show" prints a "flip.object" object
#==========================================================
setMethod("show", "flip.object", function(object)
{
  result(object)
})

setGeneric("summary")
setMethod("summary", "flip.object", function(object,star.signif=TRUE,only.p.leq=NULL,...)
{
  nperms= as.list(object@call$perms)
  #cat(" \"flip.object\" object of package flip\n")
  cat(" Call:\n ")
  cat(deparse(object@call), "\n")
  # cat(ifelse(is.null(nperms$seed),"all",""), nperms$B, ifelse(is.finite(nperms$seed),"random",""), "permutations.\n",   
	# ifelse(is.finite(nperms$seed),paste("(seed: ",is.finite(nperms$seed)," )",sep=""),"")) 
  cat(object@permSpace$B-1, "permutations.",sep=" ")
  cat("\n")
  if(!is.null(only.p.leq))  object@res=object@res[object@res[,ncol(object@res)]<=only.p.leq,,drop=FALSE]
  
  if(is.null(star.signif)) star.signif=TRUE 
  if(is.logical(star.signif)) {
    if(star.signif) {
      #takes the last column among raw and Ajusted p-value
      column=c("p-value",names(object@res)[grep("Adjust",names(object@res))])
      column=column[length(column)]
    }
  } else { column=star.signif; star.signif=TRUE}
  if(star.signif){object@res$"sig."=sapply((object@res[,column]<=.05)+
                                  (object@res[,column]<=.01)+
                                  (object@res[,column]<=.001),function(x) paste(rep("*",x),collapse=""))
  }
  result(object,...)
})


#==========================================================
# Functions to extract relevant information from 
# a flip.object object
#==========================================================
setGeneric("result", function(object,...) standardGeneric("result"))
setMethod("result", "flip.object",
  function(object,...) { 
    for(i in c("p-value",names(object@res)[grep("Adjust",names(object@res))])){
      lower=(object@res[,i]<.0001)
      object@res[,i]=formatC(object@res[,i],digits=4,drop0trailing=FALSE,format="f")
      object@res[lower,i]="<.0001"
    }
      print(object@res, digits = 4)  
})

# #==========================================================
# setGeneric("subsets", function(object, ...) standardGeneric("subsets"))
# setMethod("subsets", "flip.object", function(object, ...) {
  # object@subsets
# })


#==========================================================
setGeneric("p.value", function(object, ...) standardGeneric("p.value"))
setMethod("p.value", "flip.object",
  function(object) {
    x=object@res[,"p-value"] 
	names(x)=names(object)
	x
	
  }
)


#==========================================================
#it concatenates flip objects
#warning("More than one test statistic is imputed, only the first perms space will be stored.")
cFlip <- function(...) {
  res=list(...)[[1]]
    if(length(list(...))>1){
		nperms=sapply(list(...),function(xx) nrow(xx@permT))
		if(length(unique(nperms))>1) {
      warning("The flip-objects have different number of permutations, the minimum number will be retained for each test.")
      nperms=min(nperms)
		} else nperms=nperms[1]

    for(i in 2:length(list(...)))  res@permT=cbind(res@permT[1:nperms,,drop=FALSE],list(...)[[i]]@permT[1:nperms,,drop=FALSE])
		
		res@tail = as.vector(unlist(sapply(1:length(list(...)), function(i)  rep(if(is.null(list(...)[[i]]@tail)) 0 else list(...)[[i]]@tail,length.out=ncol(list(...)[[i]]@permT))
                      )))
		# migliore questo output, ammettere la presenza di altri elementi in extraInfoPre
		resNames=unique(as.vector(sapply(list(...),function(xx) colnames(xx@res))))
		resNames=c(setdiff(resNames,c("Stat","p-value")),c("Stat","p-value"))
		res@res[,setdiff(resNames,colnames(res@res))]=NA
		res@res=res@res[,resNames]
		for(i in 2:length(list(...)))  {
      res@res[nrow(res@res)+(1:nrow(list(...)[[i]]@res)),colnames(list(...)[[i]]@res)]=list(...)[[i]]@res
		}
	}
	res
}

#==========================================================
setGeneric("size", function(object, ...) standardGeneric("size"))
setMethod("size", "flip.object",
  function(object) {
    dim(object@permT)
  }
)
# ==========================================================
# setGeneric("dim", function(object, ...) standardGeneric("dim"))
# setMethod("dim", "flip.object",
  # function(object) {
    # c(object@nperms$B , dim(object@res)[1])
  # }
# )


#==========================================================
# The subsetting methods for "flip.object"
#==========================================================
setMethod("[", "flip.object", 
            function(x, i, j,...,drop) 
{
  iii=i
	if(is.character(i) && !all(i %in% names(x))){ 
		search=which(!(i %in% names(x)))
		extended= lapply(i, function(ii) names(x)[if(ii %in% names(x)) ii else grep(ii, names(x))] )
		i=unlist(extended)	
	}


  if (all(i %in% names(x)) || 
          all(i %in% 1:length(x)) ||
          all(i %in% -1:-length(x)) ||
          (is.logical(i) && (length(i)== length(x)))) {
    x@res <- x@res[i, ,drop=FALSE]
    if (!is.null(x@permP)) #if(i <= ncol(x@permP)) 
		x@permP <- x@permP[,i,drop=FALSE]
    if (!is.null(x@permT)) #if(i <= ncol(x@permT)) 
		x@permT <- x@permT[,i,drop=FALSE]
    if (!is.null(x@tail)) # if((i <= length(as.vector(x@tail))) || (length(as.vector(x@tail))==1))
								x@tail <- x@tail[min(length(as.vector(x@tail)),1)]
     if("data"%in%slotNames(x)) if (!is.null(x@data)) {
       #le colonne di Y non sono le stesse delle colonne di permT (a meno che non ci sia un'unica colonna X)
       search=which(!(iii %in% colnames(x@data$Y)))
       extended= lapply(iii, function(ii) colnames(x@data$Y)[if(ii %in% colnames(x@data$Y)) ii else grep(ii, colnames(x@data$Y))] )
       i=unlist(extended)        
      if(!is.null(x@data$Y)) x@data$Y <- x@data$Y[,i,drop=FALSE]
      if(!is.null(x@data$se)) x@data$se <- x@data$se[,i,drop=FALSE]
      if(!is.null(x@data$df.mod)) x@data$df.mod <- x@data$df.mod[,i,drop=FALSE]
      if(!is.null(x@data$df.res)) if(ncol(x@data$df.res)>1) x@data$df.res <- x@data$df.res[,i,drop=FALSE]
      if(!is.null(x@data$covs)) x@data$covs <- x@data$covs[,i,i,drop=FALSE]
      if(!is.null(x@data$Su)) x@data$Su <- x@data$Su[i,i,drop=FALSE]
      if(!is.null(x@data$dispersion)) if(ncol(x@data$dispersion)>1) x@data$dispersion <- x@data$dispersion[,i,drop=FALSE]
      if(!is.null(x@data$coeffWithin)) x@data$coeffWithin <- x@data$coeffWithin[,i,drop=FALSE]
    }
    #if (!is.null(x@weights)) x@weights <- x@weights[i]
    x
  } else {
    stop("invalid index set", call. = FALSE)
  }
})            

setMethod("[[", "flip.object", 
            function(x, i, j,...,exact) 
{
   x[i]
})

#==========================================================
# The length method for "flip.object"
#==========================================================
setMethod("length", "flip.object", 
            function(x) 
{
  nrow(x@res)
})



#==========================================================
# The names and alias methods for "flip.object" 
# (applies to pathwaynames)
#==========================================================
setMethod("names", "flip.object", 
            function(x) 
{
  rownames(x@res)
})      


setMethod("names<-", "flip.object", 
            function(x, value) 
{
  rownames(x@res) <- value
  if (!is.null(x@permP)) colnames(x@permP) <- value
  if (!is.null(x@permT)) colnames(x@permT) <- value
  x
})            



#==========================================================
# A sort method for "flip.object"
#==========================================================
setGeneric("sort") 
setMethod("sort", "flip.object",
  function(x, decreasing = FALSE ) {
      ix <- order(p.value(x), decreasing=decreasing)
    x[ix]
  }
)

# #==========================================================
# # Model.matrix extracts the model matrix (only if x=TRUE)
# #==========================================================
# setMethod("model.matrix", matchSignature(signature(object = "flip.object"), model.matrix),
  # function(object, ... ) {
    # list(tail = object@tail, Z = object@Z)
  # }
# )


#==========================================================
# Multiple testing correction for "flip.object" object
#==========================================================

setGeneric("p.adjust", function(p, method = p.adjust.methods, n = length(p)) standardGeneric("p.adjust"))
setMethod("p.adjust", matchSignature(signature(p = "flip.object"), p.adjust),
  function(p, method = p.adjust.methods, n = length(p)) {
    method <- method[1]
    method <- p.adjust.methods[grep(method, p.adjust.methods, ignore.case=T)]
    if(length(method)==(0))   # this is just to get a good error message
      method <- match.arg(method)
	if (missing(n))
      p@res <- cbind(p@res, p.adjust(p.value(p), method=method))
    else
      p@res <- cbind(p@res, p.adjust(p.value(p), method=method, n=n))
	
	colnames(p@res)[length(colnames(p@res))]=paste("Adjust:",method,sep="")	
	
    p
  }
)



#==========================================================
# Histogram method to visualize permutations
#==========================================================
setGeneric("hist", function(x,...) standardGeneric("hist"))
#setMethod("hist", matchSignature(signature(x = "flip.object"), hist), 
setMethod("hist", "flip.object", function(x, ...)  {
  
  flip.hist <- function(x, breaks=100, main="", xlab = "Permutation test statistics", ...) {

     if (length(x) > 1)
     stop("length(object) > 1. Please reduce to a single test result")
  
    # if (is.null(x@weights)) 
      # weights <- rep(1, size(x))
    # else
      # weights <- x@weights[[1]]
    # if (is.null(x@subsets))
      # subset <- seq_len(size(x))
    # else
      # subset <- x@subsets[[1]]
  
    #recalculate <- x@functions$permutations(subset, weights)
    
    Q <- x@permT[1,]
    nperm <- length(x@permT-1)
    hst <- hist(x@permT, xlim = c(1.1 * min(0, x@permT), 1.1 * max(x@permT)), breaks = breaks, 
      main = main, xlab = xlab, ...)
     if(is.null(list(...)$freq) & is.null(list(...)$probability)) 
       h <- max(hst$counts) else {
         if(c(1-list(...)$freq,list(...)$probability)>0 ) 
           h <- max(hst$density) else  
             h <- max(hst$counts)
       }
     
    arrows( Q, h/2, Q, 0 , lwd=2)
    text( Q, h/2, 'Observed\ntest\nstatistic' , pos=3)
  
    # No output
    invisible(list(statistic = Q, histogram = hst))
  }
  
  flip.hist(x,...)
})        


# #==========================================================
# # Graph plot for flip-oject
# #==========================================================


setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", "flip.object", 
 function(x, y, ...) {
#setMethod("plot", "flip.object", function(x, y, ...) {
  if(!exists("main")) main=NULL 
  if(!exists("xlab")) xlab = NULL
  if(!exists("ylab")) ylab=NULL 
  
  plot.flip <- function(x, y=NULL, main, xlab, ylab,...){
    #draw <- function(x, main, xlab, ylab,...){
    if (length(x)==1 ){
      hist(x, ...)
    } else if (length(x)==2 ){
      plot(x@permT, xlim=range(x@permT[,1]), ylim=range(x@permT[,2]),
           xlab=colnames(x@permT)[1],
           ylab=colnames(x@permT)[2],
           main= "Permutation Space" ,
           col="darkgrey",
           bg="orange",pch=21,lwd=1,cex=1,asp=1)
      
      points(x@permT[1,1],x@permT[1,2],col="darkgrey",bg="blue",cex=2,lwd=2,pch=21)
      text(x@permT[1,1],x@permT[1,2],labels="ObsStat")
    } else { 
      pc=prcomp(x@permT,scale. =FALSE,center=FALSE)
      #obs is always on top-right quadrant:
      pc$rotation[,1]=pc$rotation[,1]*sign(pc$x[1,1]) 
      pc$rotation[,2]=pc$rotation[,2]*sign(pc$x[1,2]) 
      pc$x[,1]=pc$x[,1]*sign(pc$x[1,1])
      pc$x[,2]=pc$x[,2]*sign(pc$x[1,2]) 
      
      pc$x=pc$x[,1:2]/pc$sdev[1:2]      
      datapc=pc$rotation[,1:2]*sqrt(nrow(pc$rotation))*1.3
      
      plot(pc$x, xlim=range(c(pc$x[,1],datapc[,1])), ylim=range(c(pc$x[,2],datapc[,2])),
           xlab=paste("PC1 (",round(pc$ sdev [1]^2 /sum(pc$ sdev ^2) *100,2)," %)",sep=""),
           ylab=paste("PC2 (",round(pc$ sdev [2]^2 /sum(pc$ sdev ^2) *100,2)," %)",sep=""),
           main= "PCA of Permutation Space" ,
           col="darkgrey",
           bg="orange",pch=21,lwd=1,cex=1,asp=1)
      
      points(pc$x[1,1],pc$x[1,2],col="darkgrey",bg="blue",cex=2,lwd=2,pch=21)

      arrows(0,0,datapc[,1],datapc[,2],col=ifelse( p.value(x)<.05,"red","blue"),lwd=2,angle=15,length=.1)
      text(datapc[,1],datapc[,2],labels=rownames(datapc))
      text(pc$x[1,1],pc$x[1,2],labels="ObsStat")
      
      # lam <- pc$sdev[1:2] #* sqrt(dim(pc$x)[1])
      # #plot(pc$x[, 1:2]/lam)
      # pc$x[,1:2]=pc$x[,1:2] / lam
      # pc$rotation[,1:2]=pc$rotation[,1:2]*lam
      # plot(pc$x[,1],pc$x[,2],lwd=1,pty="o",xlim=range(pc$x[,1])*1.2,ylim=range(pc$x[,2])*1.2,
      # xlab=paste("PC1 (",round(pc$ sdev [1]^2 /sum(pc$ sdev ^2) *100,2)," %)",sep=""),
      # ylab=paste("PC2 (",round(pc$ sdev [2]^2 /sum(pc$ sdev ^2) *100,2)," %)",sep=""),col="gray",pch=21,bg="gray")
      # points(pc$x[1,1],pc$x[1,2],col="red",lwd=3,pch=21,bg="red")
      # text(pc$x[1,1]*1.1,pc$x[1,2]*1.1,col="red","Obs")
      # arrows( 0, 0, 2*pc$rotation[,1], 2*pc$rotation[,2], lwd=1,col="gray")
      # text(2.1*pc$rotation[,1], 2.1*pc$rotation[,2], rownames(pc$rotation), cex=1.5,col="black")
      # title("PCA of Permutation Space") 
    }
  }
  plot.flip(x,y=NULL, main=main, xlab=xlab, ylab=ylab,...)
})



# #==========================================================
# # get elements of flip-object
# #==========================================================

getFlip <- function(obj,element){
  if(element%in%slotNames(obj))
    return(slot(obj,element))
  
  if(element%in%names(obj@res))
    return(obj@res[element])
  
  if(tolower(element)%in%c(tolower("Adjust:"),tolower("Adjust")))
    return(obj@res[grep("Adjust",colnames(obj@res))])
    
  if(!is.null(obj@data)){
    if(element%in%names(obj@data))
      return(obj@data[element])
    
    if(substr(element,1,5)=="data$")
      return(obj@data[substr(element,6,nchar(element))])
    
  } else if(!is.null(obj@call.env$data)){
    
    if(element%in%names(obj@data))
      return(obj@data[element])
    
    if(substr(element,1,5)=="data$")
      return(obj@data[substr(element,6,nchar(element))])    
    
  }
  
  if(substr(element,1,4)=="res$")
    return(obj@res[substr(element,5,nchar(element))])
  
  if(element%in%c("nperms","perms","B"))
    return(obj@permSpace$B)
}