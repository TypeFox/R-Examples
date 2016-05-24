###----Data--------
setClass("nExpressionSet", representation(exprs = "matrix",phenoData = "data.frame", featureData="data.frame",annotation="character"))
setClass("nxprnSet", representation(es = "nExpressionSet"))
setClass("nXprnSet", representation("nxprnSet"))
setClass("nxprnSetPair", representation(x = "nxprnSet", y = "nxprnSet",annotation="character"))
#
setGeneric("nexprs", function(object,...) standardGeneric("nexprs"))
setGeneric("nfeatureNames", function(object,...) standardGeneric("nfeatureNames"))
setGeneric("npData", function(object,...) standardGeneric("npData"))
setGeneric("nfData", function(object,...) standardGeneric("nfData"))
setGeneric("nannotation", function(object,...) standardGeneric("nannotation"))
setGeneric("nexprs<-", function(x, value) standardGeneric("nexprs<-"))
setReplaceMethod("nexprs", signature(x="nExpressionSet",value="matrix"),function(x, value) {x@exprs<-value})
setReplaceMethod("nexprs", signature(x="nxprnSet",value="matrix"),function(x, value) {x@es@exprs<-value})
setReplaceMethod("nexprs", signature(x="nxprnSet",value="numeric"),function(x, value) {x@es@exprs<-as_rowmatrix(value)})
setReplaceMethod("nexprs", signature(x="nExpressionSet",value="numeric"),function(x, value) {x@exprs<-as_rowmatrix(value)})

setGeneric("nannotation<-", function(x, value) standardGeneric("nannotation<-"))
setReplaceMethod("nannotation", signature(x="nExpressionSet",value="character"),function(x, value) {x@annotation<-value})
setReplaceMethod("nannotation", signature(x="nxprnSet",value="character"),function(x, value) {x@es@nannotation<-value})

#
setMethod("nannotation", signature(object = "nExpressionSet"), function(object){object@annotation})
setMethod("nexprs", signature(object = "nExpressionSet"), function(object){object@exprs})
setMethod("npData", signature(object = "nExpressionSet"), function(object){object@phenoData})
setMethod("nfData", signature(object = "nExpressionSet"), function(object){object@featureData})
setMethod("nfeatureNames", signature(object = "nExpressionSet"), function(object){rownames(nexprs(object))})
#
setMethod("nannotation", signature(object = "nxprnSet"), function(object){object@es@annotation})#{nannotation(as(object, "nExpressionSet"))})
setMethod("nexprs", signature(object = "nxprnSet"), function(object){object@es@exprs})#{nexprs(as(object, "nExpressionSet"))})#
setMethod("npData", signature(object = "nxprnSet"), function(object){object@es@phenoData})#{npData(as(object, "nExpressionSet"))})#
setMethod("nfData", signature(object = "nxprnSet"), function(object){object@es@featureData})#{nfData(as(object, "nExpressionSet"))})#
setMethod("nfeatureNames", signature(object = "nxprnSet"), function(object){rownames(nexprs(object))})

setMethod("nannotation", signature(object = "nxprnSetPair"), function(object){
  paste(nannotation(object@x),nannotation(object@y),sep=" - ")})
setMethod("nexprs", signature(object = "nxprnSetPair"), function(object){
	  list(x=nexprs(object@x),y=nexprs(object@x))})
setMethod("npData", signature(object = "nxprnSetPair"), function(object){
  list(x=npData(object@x),y=npData(object@y))})
setMethod("nfData", signature(object = "nxprnSetPair"), function(object){ncbind(nfData(object@x),nfData(object@y))})#{object@es@featureData})
setMethod("nfeatureNames", signature(object = "nxprnSetPair"), function(object){rownames(nexprs(object@x))})

#
setMethod("names", signature(x = "nExpressionSet"), function(x){nfeatureNames(x)})
setMethod("names", signature(x = "nxprnSet"), function(x){nfeatureNames(x)})
setMethod("names", signature(x = "nxprnSetPair"), function(x){nfeatureNames(x@x)})
setMethod("print", signature(x = "nExpressionSet"), function(x){print(str(x))})
setMethod("print", signature(x = "nxprnSet"), function(x){message(class(x), "; data:\n");print(as(x, "nExpressionSet"))})    
setMethod("print", signature(x = "nxprnSetPair"), function(x){
  message(class(x), "; data:\n");message("x:\n");print(as(x@x, "nExpressionSet"))
  message("y:\n");print(as(x@y, "nExpressionSet"))})    

setMethod("length", signature(x = "nExpressionSet"), function(x){ # like Size
	size <- sapply(1:length(names(x)), function(i){sum(is.finite(nexprs(x)[i, ]))})
	names(size) <- names(x);size
})
setMethod("length", signature(x = "nxprnSet"), function(x){length(as(x, "nExpressionSet"))})
setMethod("length", signature(x = "nxprnSetPair"), function(x){length(x@x)})

setMethod("[", signature(x = "nExpressionSet", i = "ANY", j = "missing"), function(x, i, j, drop){
  x@exprs <- x@exprs[i, ];x
})
setMethod("[", signature(x = "nExpressionSet", i = "missing", j = "ANY"), function(x, i, j, drop){
  x@exprs <- x@exprs[, j];x
})
setMethod("[", signature(x = "nExpressionSet", i = "ANY", j = "ANY"), function(x, i, j, drop){
  x@exprs <- x@exprs[i, j];x
})

setMethod("[", signature(x = "nxprnSet", i = "ANY", j = "missing"), function(x, i, j, drop){
  es <- as(x, "nExpressionSet");x@es <- es[i, ];x
})
setMethod("[", signature(x = "nxprnSet", i = "missing", j = "ANY"), function(x, i, j, drop){
  es <- as(x, "nExpressionSet");x@es <- es[, j];x
})
setMethod("[", signature(x = "nxprnSet", i = "ANY", j = "ANY"), function(x, i, j, drop){
  es <- as(x, "nExpressionSet");x@es <- es[i, j];x
})

setMethod("[", signature(x = "nxprnSetPair", i = "ANY", j = "missing"), function(x, i, j, drop){
  nexprs(x@x) <- nexprs(x@x) [i, ];nexprs(x@y) <- nexprs(x@y) [i, ];x
})

setMethod("logb", signature(x = "nExpressionSet", base = "missing"), function(x, base)
{
  nexprs(x) <- logb(nexprs(x))
  x
})
setMethod("logb", signature(x = "nxprnSet", base = "missing"), function(x, base)
{
  x@es <- logb(as(x, "nExpressionSet"))
  x
})
setMethod("logb", signature(x = "nxprnSetPair", base = "missing"), function(x, base)
{
   x@x <- logb(x@x);x@y <- logb(x@y)
  x
})
setMethod("exp", signature(x = "nExpressionSet"), function(x)
{
  nexprs(x) <- exp(nexprs(x))
  x
})
setMethod("exp", signature(x = "nxprnSet"), function(x)
{
  x@es <- exp(as(x, "nExpressionSet"))
  x
})
setMethod("exp", signature(x = "nxprnSetPair"), function(x)
{
  x@x <- exp(x@x);x@y <- exp(x@y)
  x
})

setMethod("dim", signature(x = "nxprnSet"), function(x)
{
	dim(nexprs(x))
})
setMethod("dimnames", signature(x = "nxprnSet"), function(x)
{
	dimnames(nexprs(x))
})
setMethod("dim", signature(x = "nxprnSetPair"), function(x)
{
	list(x=dim(nexprs(x@x)),y=dim(nexprs(x@y)))
})
setMethod("dimnames", signature(x = "nxprnSetPair"), function(x)
{
	list(x=dimnames(nexprs(x@x)),y=dimnames(nexprs(x@y)))
})

##
setMethod("nfeatureNames", signature(object = "nxprnSetPair"), function(object){nfeatureNames(object@x)}) # new 24 April 2008; for bias.r
setMethod("names", signature(x = "nxprnSetPair"), function(x){nfeatureNames(x)})# new 24 April 2008; for bias.r
setMethod("logb", signature(x = "nxprnSetPair", base = "missing"), function(x, base)
{
	if(is(x@x, "nXprnSet"))
	{
		x@x <- logb(x@x)
		x@y <- logb(x@y)
	}
	else
		warning("no log taken")
  x
})
setMethod("[", signature(x = "nxprnSetPair", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@x <- x@x[i, ]
	x@y <- x@y[i, ]
	stopifnot(validObject(x))
	x
})


#
setMethod("nexprs", signature(object = "nXprnSet"), function(object){nMatrix(nexprs(as(object, "nExpressionSet")))})
setMethod("print", signature(x = "nXprnSet"), function(x)
{
  message(class(x), "; ratio or nonnegative intensity data:\n")
  print(as(x, "nExpressionSet"))
})
setMethod("logb", signature(x = "nXprnSet", base = "missing"), function(x, base){logb(as(x, "nxprnSet"))})
#
####setIs("nxprnSet", "nExpressionSet")
setValidity("nExpressionSet", function(object){
        oks<-c(NA,NA,TRUE,TRUE)
        fd<-nfData(object);nfd<-nrow(fd)
        pd<-npData(object);npd<-nrow(pd)
        dd<-nexprs(object);ndr<-nrow(dd);ndc<-ncol(dd)
	oks[1]<-nfd%in%c(0,ndr)
        oks[2]<-npd%in%c(0,ndc)
        if (nfd==ndr){oks[3]<-identical(rownames(dd) == rownames(fd))}
        if (npd==ndc){oks[4]<-identical(colnames(dd) == rownames(pd))}
        if(!all(oks)){stop("Invalid nExpressionSet")}
        
})

setValidity("nxprnSet", function(object){validObject(object@es)})#setValidity("nxprnSet", function(object){all(colnames(exprs(object)) == rownames(pData(object)))})
setValidity("nXprnSet", function(object){validObject(object@es);all(as.numeric(object@es@exprs) >= 0, na.rm = TRUE)})
setValidity("nxprnSetPair", function(object)
{
	cla.ok <- !is(object@x, "nXprnSet") && !is(object@y, "nXprnSet")
	fea.x <- nfeatureNames(object@x)
	fea.y <- nfeatureNames(object@y)
	fea.ok <- length(fea.x) == length(fea.x) && is.character(fea.x) && all(fea.x == fea.y)
	ok <- cla.ok && fea.ok
	if(!ok){printInvalid(object); browser()}
	ok
})
#
setAs(from = "nExpressionSet", to = "nxprnSet", function(from)
{
	coercenExpressionSet(from = from, to.fun = new_nxprnSet)
})
setAs(from = "nExpressionSet", to = "nXprnSet", function(from)
{
	coercenExpressionSet(from = from, to.fun = new_nXprnSet) # as(as(from, "nxprnSet"), "NxprnSet")
})

setAs(from = "nxprnSetPair", to = "nxprnSet", function(from)
{
        featureData <- rbind(nfData(from@x), nfData(from@y))
	phenoData <- rbind(npData(from@x), npData(from@y))
	exprs <- cbind(nexprs(from@x), nexprs(from@y))
	new_nxprnSet(phenoData = phenoData, exprs = exprs, featureData=featureData)
})

setAs(from = "nxprnSet", to = "nExpressionSet", function(from){from@es})
setAs(from = "nxprnSet", to = "numeric", function(from){as(nexprs(from), "numeric")})
setAs(from = "nExpressionSet", to = "nxprnSet", function(from){new("nxprnSet",es=from)})
setAs(from = "nExpressionSet", to = "numeric", function(from){as(nexprs(from), "numeric")})
setAs(from = "nxprnSetPair", to = "nxprnSet", function(from){
  esx<-cbind(nexprs(from)$x,nexprs(from)$y)
  new_nxprnSet(phenoData = as.data.frame(NULL), exprs = esx,featureData=nfData(from),
	       annotation=nannotation(from))
  })

##
new_nExpressionSet<-function(x = matrix(0), phenoData = as.data.frame(NULL), featureData=as.data.frame(NULL), annotation=character(0)){
    if(is(x,"numeric")){nx<-matrix(x, length(x));colnames(nx)<-names(x);x<-nx}
    k.2dfr<-function(y){
        y<-MakeNames(y)
        if(Is(y,c("numeric","logical","character"))&&length(y)>0){
            y<-data.frame(y=y,row.names=names(y),stringsAsFactors=FALSE)}
        else if (is(y,"matrix")&&length(y)>0){y<-as.data.frame(y)}
    }
    
    assert.is(x,"matrix")
    assert.is(phenoData,"data.frame")
    assert.is(featureData,"data.frame")
    
    x<-MakeNames(x)
    if(length(phenoData)>0){
        phenoData<-MakeNames(phenoData)
        phenoData<-phenoData[,colnames(x)]
        if(nrow(phenoData)!=ncol(x)){stop("problem with names of phenoData")}
        }
    if(length(featureData)>0){
        featureData<-MakeNames(featureData);
        featureData<-featureData[rownames(x),]
        if(nrow(featureData)!=nrow(x)){stop("problem with names of featureData")}}
        
    new("nExpressionSet",exprs = x,phenoData = phenoData, featureData=featureData,annotation=annotation)
}
new_nxprnSet<-function(phenoData = as.data.frame(NULL), exprs = matrix(0),featureData=as.data.frame(NULL), annotation=character(0)){
  
  z<-new_nExpressionSet(x = exprs, phenoData = phenoData, featureData=featureData, annotation=annotation)
  new("nxprnSet",es=z)
}
new_nXprnSet<-function(phenoData = as.data.frame(NULL), exprs = matrix(0),featureData=as.data.frame(NULL), annotation=character(0)){
    #exprs<-nMatrix(exprs)
    z<-new_nExpressionSet(x = exprs, phenoData = phenoData, featureData=featureData, annotation=annotation)
    nz<-new("nXprnSet",new("nxprnSet",es=z))
}
coercenExpressionSet <- function(from, to.fun,...)
{
	assert.is(from, "nExpressionSet")
	assert.is(to.fun, "function")
	to.fun(phenoData = npData(from), exprs = nexprs(from), featureData=nfData(from), annotation=nannotation(from), ...)#April 2013: marta added featureData=fData(from)
}
#
setGeneric("nxprnSubset", function(object,...) standardGeneric("nxprnSubset"))
setMethod("nxprnSubset", signature(object = "nExpressionSet"), function(object, level, factor.name, ...){
	if(missing(factor.name))
	{
		fac.boo <- sapply(npData(object), is.factor)
		if(!any(fac.boo))
			stop("there are no factors in pData(object)")
		factor.name <- names(npData(object))[fac.boo][1]
		message("factor.name = ", factor.name)
	}
	fac <- npData(object)[, factor.name] # [,  == "MyoT"]
	stopifnot(is.factor(fac))
	if(!all(level %in% as.character(fac)))
	{ message("cannot complete xprnSubset"); browser()}
	boo <- fac %in% level
	stopifnot(sum(boo) >= 1 && length(boo) == ncol(nexprs(object)))
	es <- object[, boo]
	ann <- paste(level, collapse = "&")
	nannotation(es) <- if(length(nannotation(es)) == 1)
		paste(ann, nannotation(es), sep = " of ")
	else
		ann
	es
})
setMethod("nxprnSubset", signature(object = "nxprnSet"), function(object, ...){object@es <- nxprnSubset(object = object@es, ...);object})
#

setClassUnion(name = "nxprnSetObject", members = c("nxprnSet", "nxprnSetPair")) # used in estimate.s, biasEstimate

setGeneric("new_nxprnSetPair", function(x, y, factor.name,...) standardGeneric("new_nxprnSetPair"))
setMethod("new_nxprnSetPair", signature(x = "nxprnSet", y = "missing", factor.name = "character"), function(x, y, factor.name, level, verbose)
{
	if(missing(verbose))
		verbose <- FALSE
	fac <- npData(x)[, factor.name]
	if(!is.factor(fac)) stop("!is.factor(fac)")
	if(missing(level))
		level <- levels(unique((fac)))#as.character
	if(length(level) != 2) stop("length(level) != 2; try specifying different level argument.")
	get.subset <- function(level)
	{
		if(!all(level %in% as.character(fac)))
		{ message("cannot complete nxprnSetPair"); browser()}
		if(verbose)
			message("calling nxprnSubset with level = ", level, ", factor.name = ", factor.name)
		su <- nxprnSubset(object = x, level = level, factor.name = factor.name)
		if(verbose)
			message("finished nxprnSubset with level = ", level, ", factor.name = ", factor.name)
		su
	}
	x.sub <- get.subset(level = level[1])
	y.sub <- get.subset(level = level[2])
	new_nxprnSetPair(x = x.sub, y = y.sub)
})
setMethod("new_nxprnSetPair", signature(x = "nxprnSet", y = "nxprnSet", factor.name = "missing"), function(x, y, factor.name,...)
{
	if(is(x, "nXprnSet") || is(y, "nXprnSet"))
		stop("x, y inconsistency")
	else
		new("nxprnSetPair", x = x, y = y)
})
setMethod("new_nxprnSetPair", signature(x = "nXprnSet", y = "nXprnSet", factor.name = "missing"), function(x, y, factor.name, ...)
{
	new_nxprnSetPair(x = logb(x), y = logb(y))
})

setMethod("new_nxprnSetPair", signature(x = "matrix", y = "matrix", factor.name = "missing"), function(x, y, factor.name, fdata=as.data.frame(NULL), annotation=character(0), pdatax = as.data.frame(NULL), pdatay = as.data.frame(NULL), paired=F,rm.na=T)
{#new_nxprnSet<-function(phenoData = as.data.frame(NULL), exprs = matrix(0),featureData=as.data.frame(NULL), annotation=character(0)){
    z<-prep.2matrices(x=x,y=y,paired=paired,rm.na=rm.na)
    xx<-new_nxprnSet(phenoData = pdatax, exprs = z$x,featureData=fdata, annotation=character(0))
    yy<-new_nxprnSet(phenoData = pdatay, exprs = z$y,featureData=fdata, annotation=character(0))
    zz<-new_nxprnSetPair(x = xx, y = yy)
    zz@annotation<-annotation
    zz
})

#
setClass("nxprnSetObjects", representation("list")) # new 28 April 2008
setValidity("nxprnSetObjects", function(object){all(sapply(object, is, class2 = "nxprnSetObject"))})
setClass("nxprnSetObjectPair", representation(training = "nxprnSetObject", test = "nxprnSetObject")) # new 24 April 2008; for bias.r
setValidity("nxprnSetObjectPair", function(object)
{
	cla.ok <- TRUE
	fea.training <- nfeatureNames(object@training)
	fea.test <- nfeatureNames(object@test)
	fea.ok <- length(fea.training) == length(fea.training) && is.character(fea.training) && all(fea.training == fea.test)
	ok <- cla.ok && fea.ok
	if(!ok){printInvalid(object); browser()}
	ok
})
setMethod("nfeatureNames", signature(object = "nxprnSetObjectPair"), function(object){nfeatureNames(object@test)}) # new 24 April 2008; for bias.r
nprnSet2matrix<-function(x,y=NULL,paired=FALSE){rm.na<-T
	assert.is(x,c('nxprnSet','nxprnSetPair','matrix'))
	assert.is(y,c('nxprnSet','NULL','matrix'))
	if(is(x,'nprnSet')){x<-nexprs(x)}
	if(is(y,'nprnSet')){y<-nexprs(y)}
	if(is(x,'nprnSetPair')){y<-nexprs(x)$y;x<-nexprs(x)$x}
	y<-MakeNames(y);x<-MakeNames(x)

	if(!is_vide(y)&&paired==T){
		z<-prep.2matrices(x=x,y=y,paired=T,rm.na=rm.na)
		x<-z$x;y<-z$y;paired<-z$info$paired
		pok<-is_paired(x=x,y=y)
		if(!pok&&paired==T){stop('data is not paired')}
		if(paired==T){x<-x-y;y<-NULL;paired<-FALSE}}
	assert.is(x,'matrix')
	assert.is(y,c('matrix','NULL'))
	
	list(x=x,y=y,paired=paired)}
#----------------------
setGeneric("removeMissing", function(object,...) standardGeneric("removeMissing"))
setMethod("removeMissing", signature(object = "nxprnSet"), function(object)
{
  missing.gene <- as.logical(apply(nexprs(object), 1, function(ro){all(is.na(ro))}))
  stopifnot(length(missing.gene) == nrow(nexprs(object)))
  if(all(missing.gene))
  {
  	message("all genes are missing in ", class(object))
  	browser()
  }
  if(any(missing.gene))
  	object[!missing.gene]
  else
	  object
})

setMethod("removeMissing", signature(object = "nxprnSetPair"), function(object)
{
  x <- removeMissing(object@x)
  y <- removeMissing(object@y)
  nam <- intersect(nfeatureNames(object@x), nfeatureNames(y))
  if(!is.character(nam) || length(nam) == 0)
  { message("x and y have no non-missing genes in common"); browser()}
  object@x <- x[nam]
  object@y <- y[nam]
  if(!validObject(object))
  { message(class(object), " is no longer valid"); browser()}
  object
})

 
#
#Rep <- function(object, times, ...){stop(paste("bad class for Rep (", class(object), ") or times (", class(times), ")", sep = ""))}
setGeneric("Rep", function(object,...) standardGeneric("Rep"))
setMethod("Rep", signature(object = "Vector"), function(object, times, ...)
{
	nam <- RepNames(object = object, times = times, ...)
	vec <- rep(object, times)
	stopifnot(length(vec) == length(nam))
	names(vec) <- nam
	vec
})
setMethod("Rep", signature(object = "matrix"), function(object, times, ...)
{
	nam <- RepNames(object = object, times = times, ...)
	mat <- matrix(rep(as.numeric(object), times), ncol = ncol(object), byrow = TRUE)
	stopifnot(nrow(mat) == times * nrow(object))
	stopifnot(nrow(mat) == length(nam))
	rownames(mat) <- nam
	colnames(mat) <- colnames(object)
	mat
})
setMethod("Rep", signature(object = "nxprnSet"), function(object, times, ...)
{
	datf <- npData(object)
	stopifnot(is(datf, "data.frame"))
	mat <- Rep(object = nexprs(object), times = times, ...)
	stopifnot(is.character(rownames(mat)) && !any(duplicated(rownames(mat))))
	es <- new_nxprnSet(phenoData = datf, exprs = mat, annotation = "normal.rxprnSet")
	nam <- names(es)
	stopifnot(length(nam) == times * length(names(object)) && all(nam == names(es)))
	es
})



#

#------------------------------------------------------------------
sorted <- function(object, ...)
{
	Slot(object = object, name = "sorted", ...)
}



##-------