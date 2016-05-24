#' Adds phenotypic variables to gwaa.data object
#' 
#' Adds phenotypic variables to \code{phdata} part 
#' of an \code{\link{gwaa.data-class}} object
#' 
#' If "newphdata" is a data frame, it is simply merged to 
#' the \code{phdata} part of the "data", and is sorted 
#' according to the right order. In this case,  
#' The "newphdata" frame should contain single variable 
#' named "id", preferably of "character" class. It may 
#' contain "sex" variable, but that will be re-named 
#' to avoid duplication with the default sex variable
#' presented in \code{phdata}.
#' 
#' If 'newphdata' is a vector, it should be of the 
#' same length as the number of people in the 'data' 
#' and is assumed to have the same order. In this case,
#' you also need to supply the name of the new phenotype 
#' via the 'name' argument 
#' 
#' @param data an object of \code{\link{gwaa.data-class}}
#' @param newphdata data frame or a vector with new 
#' phenotypic data
#' @param name if 'newphdata' is a vector, the name
#' of new variable should be specified in 'name'
#'
#' @return An (updated) object of \code{\link{gwaa.data-class}}
#' 
#' @author Yurii Aulchenko
#' 
#' @seealso \code{\link{merge.gwaa.data}} \code{\link{merge.snp.data}}
#' 
#' @keywords manip
#' 
#' @examples 
#' data(srdta)
#' # take a small subset for this example
#' srdta <- srdta[1:10,1:5]
#' srdta
#' # add single var
#' rnd <- rnorm(nids(srdta))
#' srdta1 <- add.phdata(srdta,rnd,name="random")
#' srdta1
#' # add > 1 var
#' # generate id names
#' ids <- paste("p",c(2,1,7,3,5,9,11,22,27),sep="")
#' # generate some random trait values
#' newtra1 <- rnorm(9)
#' newtra2 <- rnorm(9)
#' # make data frame
#' x <- data.frame(id=ids,newtra1=newtra1,newtra2=newtra2)
#' x
#' # now add this new trait to the data
#' srdta1 <- add.phdata(srdta,x)
#' srdta1
#' 

"add.phdata" <- 
		function(data,newphdata,name) {
	if (!is(data,"gwaa.data")) stop("data argument must be of gwaa.data-class")
	if (is(newphdata,"vector")) {
		if (missing(name)) stop("'name' should be provided if newphdata is a vector")
		if (!is(name,"character")) stop("'name' should be character")
		if (length(name) != 1) stop("only one name should be specified")
		if (length(newphdata) != dim(data@phdata)[1]) stop("dimensions of newphdta and data do not match")
		tmp <- data@phdata[,c("id","sex")]
		tmp[,name] <- newphdata
		newphdata <- tmp
		newphdata <- newphdata[,c(1,3)]
	}
	if (!is(newphdata,"data.frame")) stop("newphdata argument must be of data.frame-class")
	nidcols <- sum(names(newphdata) %in% "id")
	if (nidcols==0) stop("can not find \"id\" column in newphdata")
	if (nidcols>1) stop("more than one \"id\" column in newphdata")
	if (length(unique(newphdata$id)) != length(newphdata$id)) stop("duplicated id names in newphdata")
	if (class(newphdata$id) != "character") {
		warning("newphdata id variable does not have character class. Converted")
		newphdata$id <- as.character(newphdata$id)
	}
	
	oldph <- data@phdata
	if (length(unique(oldph$id)) != length(oldph$id)) stop("duplicated id names in oldph!!!")
	if (class(oldph$id) != "character") {
		stop("oldph id variable does not have character class!!!")
	}
	nadded <- sum((oldph$id %in% newphdata$id),na.rm=T)
	if (nadded < length(oldph$id)) 
		warning(paste("Of",length(oldph$id),"IDs present in data, only",nadded,"found in newphdata"))
	newph <- merge(oldph,newphdata,by="id",all.x=T,all.y=F)
	if (dim(newph)[1] != dim(oldph)[1]) stop("something terribly wrong in merge...")
	newph <- newph[match(oldph$id,newph$id),]
	rownames(newph) <- newph$id
	
	if (is.na(match("sex",names(newph))))
		if (!is.na(match("sex.x",names(newph)))) {
			nvar <- match("sex.x",names(newph))
			names(newph)[nvar] <- "sex"
		}
	
	data@phdata <- newph
	data
}
