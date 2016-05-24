
 
facetshade <- function( data, mapping, f, geom, geom.mapping, bg.all = TRUE, keep.orig = FALSE, ... ){

body <-  ggplot
if(missing(mapping)){
	mapping <- aes()
}
if(missing(geom)){
	geom <- NULL
}
	p <- body
	
gnames <- switch(class(f)[1],
	formula = rownames(attr(terms(f),"factors")),
	wrap = sapply(f$facets[1:2],toString),
	grid = sapply(f[1:2],toString))
if(inherits(f,'formula')){
	f <- facet_grid(f)
}
gnames <- gnames[which(gnames %in% names(data))]

	ind <- which(names(data) %in% gnames)
	
	ord <- do.call(order,data[,ind,drop=FALSE])
	data <- data[ord,]
	
	mdata <- subtable(data,ind)
	gs <- mdata$Freq
	ng <- nrow(mdata)
	mdata <- mdata[rep(1:ng, each=nrow(data)),]
	
	
	xn <- toString(mapping$x)
	yn <- toString(mapping$y)
	
	mdata[xn] <- data[xn]
	
	if(yn != ""){
		mdata[yn] <- data[yn]
	}	
vars <- sapply(mapping, toString)
	for(v in vars){
		if(!v %in% names(mdata)){
			mdata[v] <- rep(unlist(data[v]),ng)
		}
}
	
	if(keep.orig){
		for(i in gnames){
			mdata[paste("orig",i,sep=".")] <- rep(unlist(data[i]),ng)
		}
	}
	if(!bg.all){
		#mdata <- subset(mdata,gv == 0)
		cgs <- cumsum(gs)
		cgs <- cbind(c(1,cgs[-ng]+1),cgs) + (0:(ng-1))*nrow(data)
		rm.ind <- unlist( apply(cgs,1,function(z) z[1]:z[2]) )
		mdata <- mdata[-rm.ind,]
	}
	if(!is.null(geom)){
		# facetshade as first layer using the specified geom
		if(missing(geom.mapping)){
			shade.layer <- geom(data=mdata, ... )
		}else{
			shade.layer <- geom(data=mdata, mapping = geom.mapping, ... )
		}
		p <- body(data=data, mapping = mapping) + shade.layer
	}else{
		p <- body(data=mdata, mapping = mapping, ...)
	}
	p <- p + f + guides(colour = guide_legend(title=NULL))
	return(p)
}

duplicate.gg<-function(x) {
    #require(proto)
    if (requireNamespace("proto", quietly = TRUE)) {
		 r<-x
    r$layers <- lapply(r$layers, function(x) {
        proto::as.proto(as.list(x), parent=x)
    })
    r
	}else{
		stop("Please install package 'proto' to use this reordering function.")
	}
   return(r)
}

duplicate.layer<-function(x) {
    #require(proto)
    if (requireNamespace("proto", quietly = TRUE)) {
		r <-   proto::as.proto(as.list(x), parent=x)
	}else{
		stop("Please install package 'proto' to use this reordering function.")
	}
	return(r)	 
}

 
facetshade2 <- function( gg, copy.layer = 1, geom, mapping, bg.all = TRUE, keep.orig = FALSE,  ... ){
if(missing(geom)){
	geom <- NULL
}
data <- gg$data
k <- 0
while(length(data)==0){
	data <- gg$layers[[k <- k+1]]$data
}
num.layers <- length(gg$layers)
if(missing(mapping)){
	mapping <- aes()	
}
	mapping <- c(mapping,gg$layers[[copy.layer]]$mapping,gg$mapping)
	mapping <- mapping[!duplicated(names(mapping))]
	class(mapping) <- 'uneval'

	xn <- toString(mapping$x)
	yn <- toString(mapping$y)
	stopifnot(xn != '')

# get facet variables
facets <- gg$facet$facets
if(is.null(facets)){
	# facet grid was used (need better solution here)
	facets <- gg$facet
}
	gnames <- sapply(facets[1:2],toString)
	gnames <- gnames[gnames != '']

# remove all mappings for facet variables
exclude <- sapply(mapping,toString) %in% gnames
mapping <- mapping[!exclude]
class(mapping) <- 'uneval'
	
	ind <- which(names(data) %in% gnames)
	ord <- do.call(order,data[,ind,drop=FALSE])
	
	data <- data[ord,]
	mdata <- subtable(data,ind)

	gs <- mdata$Freq
	ng <- nrow(mdata)

	mdata <- mdata[rep(1:ng, each=nrow(data)),]
	mdata[xn] <- data[xn]
	if(yn != ""){
		mdata[yn] <- data[yn]
	}	
vars <- sapply(mapping, toString)
	for(v in vars){
		if( (!v %in% names(mdata)) & (v %in% names(data))){
			mdata[v] <- rep(unlist(data[v]),ng)
		}
}
	if(keep.orig){
		for(i in gnames){
			mdata[paste("orig",i,sep=".")] <- rep(unlist(data[i]),ng)
		}
	}
	if(!bg.all){
		#mdata <- subset(mdata,gv == 0)
		cgs <- cumsum(gs)
		cgs <- cbind(c(1,cgs[-ng]+1),cgs) + (0:(ng-1))*nrow(data)
		rm.ind <- unlist( apply(cgs,1,function(z) z[1]:z[2]) )
		mdata <- mdata[-rm.ind,]
	}
	
if(is.null(geom)){
	shade.layer <- duplicate.layer(gg$layers[[copy.layer]])
	shade.layer$mapping <- mapping
	shade.layer$data <-  mdata
	params <- c(...)
	for(i in names(params)){
		shade.layer$geom_params[i] <-  params[i]
	}
}else{
	shade.layer <- geom(data = mdata, mapping = mapping, ... )
}
	gg$layers <- c(shade.layer, gg$layers )
	return(gg)
}			


# facetshade3 <- function( gg, geom, bg.all = TRUE, keep.orig = FALSE,  ... ){

# if(is.null(bg.all)){
	# bg.all <- TRUE
# }
# if(is.null(keep.orig)){
	# keep.orig <- FALSE
# }
# data <- gg$data
# k <- 0
# while(length(data)==0){
	# data <- gg$layers[[k <- k+1]]$data
# }
# print(summary(data))
	# num.layers <- length(gg$layers)
	# mapping <- c(geom$mapping,gg$mapping)
	# mapping <- mapping[!duplicated(names(mapping))]
	# class(mapping) <- 'uneval'
# print(mapping)

	# xn <- toString(mapping$x)
	# yn <- toString(mapping$y)
	# stopifnot(xn != '')

# # get facet variables
# facets <- gg$facet$facets
# if(is.null(facets)){
	# # facet grid was used (need better solution here)
	# facets <- gg$facet
# }
	# gnames <- sapply(facets[1:2],toString)
	# gnames <- gnames[gnames != '']

# # remove all mappings for facet variables
# exclude <- sapply(mapping,toString) %in% gnames
# mapping <- mapping[!exclude]
# class(mapping) <- 'uneval'
	
	# ind <- which(names(data) %in% gnames)
	# ord <- do.call(order,data[,ind,drop=FALSE])
	
	# data <- data[ord,]
	# mdata <- subtable(data,ind)

	# gs <- mdata$Freq
	# ng <- nrow(mdata)

	# mdata <- mdata[rep(1:ng, each=nrow(data)),]
	# mdata[xn] <- data[xn]
	# if(yn != ""){
		# mdata[yn] <- data[yn]
	# }	
# vars <- sapply(mapping, toString)
	# for(v in vars){
		# if( (!v %in% names(mdata)) & (v %in% names(data))){
			# mdata[v] <- rep(unlist(data[v]),ng)
		# }
# }
	# if(keep.orig){
		# for(i in gnames){
			# mdata[paste("orig",i,sep=".")] <- rep(unlist(data[i]),ng)
		# }
	# }
	# if(!bg.all){
		# #mdata <- subset(mdata,gv == 0)
		# cgs <- cumsum(gs)
		# cgs <- cbind(c(1,cgs[-ng]+1),cgs) + (0:(ng-1))*nrow(data)
		# rm.ind <- unlist( apply(cgs,1,function(z) z[1]:z[2]) )
		# mdata <- mdata[-rm.ind,]
	# }
	# print(summary(mdata))
	# geom <- geom %+% mdata
	# print(geom)
	# gg <- gg + geom
	# return(gg)
# }

#  "+.gg" <- function(x, fs) {
#    if(exists('shade',envir = fs,inherits = FALSE)){
#      # create a shade layer
#    }
#    bg.all <- fs$stat_params$bg.all
#    keep.orig <- fs$stat_params$keep.orig
#    x <- facetshade3(x, geom = fs, bg.all = bg.all, keep.orig = keep.orig )
#  }
# "+.gg" <- function(e1, e2) {
#   # Get the name of what was passed in as e2, and pass along so that it
#   # can be displayed in error messages
#   e2name <- deparse(substitute(e2))
#   
#   if      (is.theme(e1))  add_theme(e1, e2, e2name)
#   else if (is.ggplot(e1)) add_ggplot(e1, e2, e2name)
# }

