devSVGTips <- function (file = "Rplots.svg", width = 7, height = 7,
                        bg = "white", fg = "black", onefile = TRUE,
                        xmlHeader = FALSE, useStyleAttributes=FALSE,
                        toolTipMode=1, toolTipFontSize=10, toolTipOpacity=1.0,
                        title="R SVG Plot", sub.special=TRUE)
{
    dev <- .C("do_SVG", as.character(file),
              as.character(bg), as.character(fg),
              as.double(width), as.double(height),
              as.logical(FALSE), as.logical(xmlHeader),
              as.character(encodeSVGSpecialChars(as.character(title), sub.special)),
              as.integer(toolTipMode),
              as.integer(toolTipFontSize), as.double(toolTipOpacity),
              as.logical(onefile), as.logical(useStyleAttributes),
              PACKAGE="RSVGTipsDevice")

      invisible(dev)
}

setSVGShapeContents <- function(contents)
{
    if (names(dev.cur()) == "devSVG")
        .C("SetSvgShapeContents", as.character(contents), PACKAGE="RSVGTipsDevice")
    invisible(NULL)
}

setSVGShapeURL <- function(url, target=NULL, sub.special=TRUE)
{
    if (names(dev.cur()) == "devSVG") {
        .C("SetSvgShapeURL", encodeSVGSpecialChars(as.character(url), sub.special=sub.special, xent=TRUE), PACKAGE="RSVGTipsDevice")
        if (!is.null(target))
            .C("SetSvgShapeURLTarget", encodeSVGSpecialChars(as.character(target), sub.special=sub.special, xent=TRUE), PACKAGE="RSVGTipsDevice")
    }
    invisible(NULL)
}

setSVGShapeToolTip <- function(title=NULL, desc=NULL, desc1=desc, desc2=NULL, sub.special=TRUE) {
    if (names(dev.cur()) != "devSVG")
        return(invisible(NULL))
    contents <- character(0)
    toolTipMode <- getSVGToolTipMode()
    if (toolTipMode>0) {
        if (!is.null(title))
            contents <- c(contents, paste("<title>", encodeSVGSpecialChars(title, sub.special=sub.special), "</title>", sep=""))
        if (toolTipMode==1) {
            if (!is.null(desc1))
                contents <- c(contents, paste("<desc>", encodeSVGSpecialChars(desc1, sub.special=sub.special), "</desc>", sep=""))
        } else {
            if (!is.null(desc1))
                contents <- c(contents, paste("<desc1>", encodeSVGSpecialChars(desc1, sub.special=sub.special), "</desc1>", sep=""))
            if (!is.null(desc2))
                contents <- c(contents, paste("<desc2>", encodeSVGSpecialChars(desc2, sub.special=sub.special), "</desc2>", sep=""))
        }
        if (length(contents)) {
            .C("SetSvgShapeContents", as.character(paste(contents, collapse="\n")), PACKAGE="RSVGTipsDevice")
        }
    }
    invisible(NULL)
}

encodeSVGSpecialChars <- function(x, sub.special=TRUE, xent=FALSE) {
    if (!sub.special)
        return(x)
    if (regexpr("[<>&'\"]", x) > 0) {
        # use a negative look-ahead assertion to avoid replacing the "&" in an XML entity
        # XML entities are things like "&amp;" (ampersand) and "&#x3b1;" ('alpha' in Greek font)
        x <- gsub("&(?!#?[A-Za-z0-9]+;)", if (xent) "&#38;" else "&amp;", x, perl=TRUE)
        x <- gsub("<", if (xent) "&#60;" else "&lt;", x, fixed=TRUE)
        x <- gsub(">", if (xent) "&#62;" else "&gt;", x, fixed=TRUE)
        x <- gsub("'", if (xent) "&#39;" else "&apos;", x, fixed=TRUE)
        x <- gsub("\"", if (xent) "&#34;" else "&quot;", x, fixed=TRUE)
    }
    x
}

getSVGToolTipMode <- function()
{
    if (names(dev.cur()) != "devSVG")
        return(-1)
    return(.C("GetSvgToolTipMode", mode=integer(1), PACKAGE="RSVGTipsDevice")$mode)
}

#SvgDevicePoint <- function(x,y)
#{
#       tmp <- .C('GetSvgDevicePoint',as.double(x),as.double(y))
#
#       c(tmp[[1]],tmp[[2]])
#}
#
#SvgUserPoint <- function(x,y)
#{
#       tmp <- .C('GetSvgUserPoint',as.double(x),as.double(y))
#
#       c(tmp[[1]],tmp[[2]])
#}
#
#SvgDeviceBoundry <- function()
#{
#       w <- 0.0
#       h <- 0.0
#       tmp <- .C('GetSvgDeviceBoundry',w,h)
#
#       c(tmp[[1]],tmp[[2]])
#}
#
#SvgDevicePoints <- function(x,y)
#{
#       tmp <- .C('GetSvgDevicePoints',as.double(x),
#                               as.double(y),n=as.integer(length(x)))
#
#       list(x=tmp[[1]],y=tmp[[2]])
#}
#
#MetaSvg <- function(size, box.size, points, values)
#{
#       str <- paste("<MetaSvg width=\"",size[1],"\" height=\"",
#              size[2],"\" rect.width=\"",box.size[1],"\" rect.height=\"",
#              box.size[2],"\">",sep="")
#
#       if(length(points) > 0){
#              buf <- paste("<point x=\"",points$x,"\" y=\"",
#              points$y,"\" value=\"",values,"\"/>",sep="",collapse="\n")
#
#       }
#
#       str <- paste(str,buf,"</MetaSvg>",sep="")
#}










