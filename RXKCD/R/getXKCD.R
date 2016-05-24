read.xkcd <- function(file = NULL)
{
  if(!is.null(file) && file.exists(file)) {
    xkcd <- file
  } else {
    path <- system.file("xkcd", package = "RXKCD") # fix requested by Brian Ripley
    datafiles <- list.files(path)
    if(!is.null(file) && file.exists(file.path(path, file))) {
      xkcd <- file.path(path, file)
    } else {
      if(!is.null(file)) stop("sorry, ", sQuote(file), " not found")
      file <- datafiles
      xkcd <- file.path(path, file)
    }
  }
  out <-read.csv(xkcd)
  return(out)
}
#'
#' Update the XKCD database saved in the user directory
#'
#' This function update the local version of the XKCD database used by searchXKCD
#'
#' @references http://xkcd.com/license.html
#'
#' @export
#'
updateConfig <- function(){
	home <- Sys.getenv("HOME") # user's home directory
	if( !file.exists( paste(home, ".Rconfig/rxkcd.rda", sep="/") ) ) {
		stop("Use saveConfig() to save your xkcd database locally!")
		} else load( paste(home, ".Rconfig/rxkcd.rda", sep="/") )
	from <- dim(xkcd.df)[[1]]
	current <- getXKCD("current", display=F)
	if ( current$num == xkcd.df$id[dim(xkcd.df)[[1]]] ) stop("Your local xkcd is already updated!")
	tmp <- NULL
	for( i in c((from+1):(current$num)) ){
		if (is.null(tmp)) tmp <- getXKCD(i, display=F)
		else tmp <- rbind(tmp, getXKCD(i, display=F))
	}
	suppressWarnings(tmp <- data.frame(tmp))
	row.names(tmp) <- tmp$num
	xkcd2add <- cbind( 
	"id"=unlist(tmp[["num"]]),
	"img"=unlist(tmp[["img"]]),
	"title"=unlist(tmp[["title"]]),
	"month"=unlist(tmp[["month"]]),
	"num"=unlist(tmp[["num"]]),
	"link"=unlist(tmp[["link"]]),
	"year"=unlist(tmp[["year"]]),
	"news"=unlist(tmp[["news"]]),
	"safe_title"=unlist(tmp[["safe_title"]]),
	"transcript"=unlist(tmp[["transcript"]]),
	"alt"=unlist(tmp[["alt"]]),
	"day"=unlist(tmp[["day"]]) 
	)
	suppressWarnings(xkcd2add <- data.frame(xkcd2add))
	row.names(xkcd2add) <- xkcd2add$num
	xkcd.updated <- rbind(xkcd.df,xkcd2add)
	xkcd.df <- xkcd.updated
	# write.csv(xkcd.updated,file="xkcd.csv",row.names=F)
	save( xkcd.df, file=paste(home, ".Rconfig/rxkcd.rda", sep="/") , compress=TRUE)
}
#'
#' Save XKCD database info into a file in the user directory
#'
#' This function saves the xkcd database as a file in the user's home directory
#'
#' @references http://xkcd.com/license.html
#'
#' @export
#'
saveConfig <- function(){
	home <- Sys.getenv("HOME") # home dir of the user
	if( file.exists( paste(home, ".Rconfig/rxkcd.rda", sep="/") ) ) stop("Use updateConfig() for updating your local xkcd database")
	else {
		dir.create( paste(home, ".Rconfig", sep="/") )
		xkcd.df <- read.xkcd()
		save( xkcd.df, file=paste(home, ".Rconfig/rxkcd.rda", sep="/") , compress=TRUE)
	}
}
#'
#' Search your favorite XKCD comic strip by title/trascript
#'
#' This function use grep to inspect the title and trascript for all the occurrences of a specified string and return a data.frame with both the number and the title of the XKCD comic strips.
#'
#' @param which string.
#' @param xkcd.df A character string giving a xkcd file in csv format. By default the csv file in the data directory of the xkcd package are used.
#'
#' @return a data.frame containing the following fields: \itemize{
#' \item num The num of the XKCD comic strip
#' \item title The title of the XKCD comic strip
#' }
#'
#' @references http://xkcd.com/license.html
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library("RXKCD")
#' searchXKCD(which="significant") 
#' searchXKCD(which="someone is wrong") }
#'
searchXKCD <- function(which="significant"){
	xkcd.df <- NULL # Thanks to Duncan Murdoch
	home <- Sys.getenv("HOME") # user's home directory
	if( file.exists( paste(home, ".Rconfig/rxkcd.rda", sep="/") ) ) {
		load( paste(home, ".Rconfig/rxkcd.rda", sep="/") )
		xkcd.df <- xkcd.df
	} else	xkcd.df <- read.xkcd()
	if(is.character(which)) {
		if(length(which) > 1) which <- sample(which)
	which.tt <- grep(which, xkcd.df["title"][[1]], ignore.case = TRUE, useBytes = TRUE)
	which.tr <- grep(which, xkcd.df["transcript"][[1]], ignore.case =TRUE, useBytes = TRUE)
	which.all <- unique(c(which.tr, which.tt))
	} 
	out <- data.frame(num=xkcd.df[which.all, "num"], title=xkcd.df[which.all, "title"])
	return(out)	
}
#'
#' Display your favourite XKCD comic in R
#'
#' This function fetches a XKCD comic strip (randomly or by number) and displays it on screen.
#'
#' @param which string: either "current" or "random"; or a number indicating the specific strip.
#' @param display  logical; TRUE (default) if you like to display the strip on the screen
#' @param html logical; TRUE if you like to open the XKCD web page for the selected comic in your browser: 
#' if TRUE it sets display and saveImg arguments to FALSE. Default FALSE
#' @param saveImg logical; TRUE if you want to save image in the current directory. Default FALSE
#'
#' @return a list containing the following fields: \itemize{
#' \item imgURL of the XKCD comic strip image (png)
#' \item title Title of the XKCD comic strip
#' \item month
#' \item numNumber of the XKCD comic strip
#' \item link
#' \item year Year of publication
#' \item safe_title
#' \item transcript
#' \item alt
#' \item day
#' }
#'
#' @references http://xkcd.com/license.html
#'
#' @export
#'
#' @examples
#'
#' library("RXKCD")
#' significant <- getXKCD(882, display=FALSE)
#'
getXKCD <- function(which = "current", display = TRUE, html = FALSE, saveImg = FALSE) {
	if (which=="current") xkcd <- fromJSON("http://xkcd.com/info.0.json")
	else if(which=="random"|which=="") {
		current <- fromJSON("http://xkcd.com/info.0.json")
		num <- sample(1:current["num"][[1]], 1)
		xkcd <- fromJSON(paste("http://xkcd.com/",num,"/info.0.json",sep=""))
	} 
	else xkcd <- fromJSON(paste("http://xkcd.com/",which,"/info.0.json",sep=""))
	class(xkcd) <- "rxkcd"	
	if(html) {
		display= FALSE
		browseURL( paste("http://xkcd.com/", as.numeric(xkcd["num"][[1]]),sep="") ) 
	}
	if (display|saveImg) {
		if(grepl(".png",xkcd["img"][[1]])){
			download.file(url=xkcd["img"][[1]], quiet=TRUE, mode="wb", destfile=paste(tempdir(),"xkcd.png",sep="/"))
			xkcd.img <- readPNG( paste(tempdir(),"xkcd.png",sep="/") )
		}
		else if(grepl(".jpg",xkcd["img"][[1]])){
			download.file(url=xkcd["img"][[1]], quiet=TRUE, mode="wb", destfile=paste(tempdir(),"xkcd.jpg",sep="/"))
			xkcd.img <- readJPEG( paste(tempdir(),"xkcd.jpg",sep="/") )
		} else stop("Unsupported image format! Try html = TRUE")
		# show the image if the format is supported
		if(display){
			max.dim = max(dim(xkcd.img))
			plot(1:max.dim, type="n", axes=F, xaxt="n",yaxt="n",xlab="",ylab="")
			rasterImage(xkcd.img, xleft=0, ybottom=0, xright=dim(xkcd.img)[[2]], ytop=dim(xkcd.img)[[1]])
		}
		# save the image
		if(saveImg) writePNG( image=xkcd.img, target=paste(xkcd$title,".png",sep="") )
	}
	return(xkcd)
}

print.rxkcd <- function(x, ...){
	cat("image.url = ", x$img, "\n", sep="")
	cat("title =  ", x$title, "\n", sep="")
	cat("num = ", x$num, "\n", sep="")
	cat("year = ", x$year, "\n", sep="")
	cat("transcript = ", x$transcript,"\n", sep="")
	cat("alt = ", x$alt, "\n", sep="")
}
