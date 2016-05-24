#'
#' plots result of the voter transition table
#' 
#' @description
#' plots results of the voter transition table. Supports absolute values, relative values
#' and negative values 
#' 
#' @param x \code{matrix} with results
#' @param abs \code{TRUE} if values are not between [0,1] (Default is \code{FALSE})
#' @param bgColors vector with 3 elements: starting colour of new colour palette,
#'          ending colour of new colour palette, length of colour palette
#' @param rectOpts named \code{list} of options for rectangle
#' @param cellOpts named \code{list} of options for inner cell text
#' @param dimNameOpts named \code{list} of options for dimnames
#' 
#' @seealso
#' \code{\link[eiwild]{getBalance}}
#' 
#' 
#' @examples
#' \dontrun{
#' # loading some fake election data
#' data(topleveldat)
#' form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
#' set.seed(1234)
#' res <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
#'                   sample=1000, thinning=2, burnin=100,verbose=100)
#' res2 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
#'                  sample=1000, thinning=2, burnin=100,verbose=100)
#' 
#' 
#' tabs <- summary(res)
#' tabs2 <- summary(res2)
#' plotResult(round(tabs$relative,3))
#' plotResult(tabs$absolut, abs=TRUE)
#' bal <- getBalance(tabs$absolut, which=c("c","GRUN_2"))
#' plotResult(bal, abs=TRUE)
#' 
#' plotResult(round(tabs$relative,3), bgColors=c("white", "darkorange", 9))
#' plotResult(round(tabs$relative,3), bgColors=c("white", "darkorange", 5))
#' 
#' plotResult(round(tabs$relative,3) - round(tabs2$relative, 3), abs=TRUE,
#'            bgColors=c("white", "darkorange", 9))
#' 
#' # ugly ;)
#' plotResult(round(tabs$relative,3), bgColors=c("blue", "red", 5)) 
#' #' }
#' 
#' 
#' @export
#'

plotResult <- function(x, abs=FALSE, bgColors=c("white", "steelblue", 10),
                       rectOpts=list(border=grey(0.7)),
                       cellOpts=list(cex=1),
                       dimNameOpts=list(col=grey(0.2), cex=.8)){
    
        # making colours
    colLength <- as.numeric(bgColors[3])
    if(colLength > 256) 
      warning("Cannot generate a color pallette with more than 256 colors", .call=FALSE)
    if (abs == TRUE) {
      if(min(x) < 0){ 
        # if there are negative values to plot x is replaced with abs(x) to have colors, too
        xSafe <- x
        x <- abs(x)
      }
      steps <- (max(x) - min(x))/colLength
      minx <- min(x) - 0.01
      brkPoints <- c(minx, steps * 1:colLength)
      brkPoints[colLength + 1] <- max(x) + 0.01
    }
    else {
      steps <- 1/colLength
      minx <- -0.01
      brkPoints <- c(minx, steps * 1:colLength)
      brkPoints[colLength + 1] <- 1.01
    }
    colCatIndiz <- as.integer(cut(x, brkPoints))
    colorF <- colorRampPalette(bgColors[1:2])
    cols <- colorF(colLength)
    
    bgCols <- rep("#00000000", length(x))
    ### right colour for every point
    for (j in 1:length(x))
        bgCols[j] <- cols[colCatIndiz[j]]
    bgCols <- matrix(bgCols, nrow(x), ncol(x))
    
    if(exists("xSafe"))
      x <- xSafe
    
        # making plot
    plot.new()
    plot.window(xlim=c(0,ncol(x)+1), ylim=c(nrow(x)+1,0))
        # making background colours
    for(rr in 1:nrow(x)) 
        for(cc in 1:ncol(x)){
            do.call("rect", c(xleft=cc, xright=cc+1, ytop=rr, ybottom=rr+1,
                              col=bgCols[rr,cc], rectOpts))
            do.call("text", c(cc+0.5, rr+.5, x[rr,cc], adj=.5,
                              cellOpts))
        }
    for(rr in 1:nrow(x))
        do.call("text", c(0.5, rr+0.5, rownames(x)[rr], dimNameOpts))
    for(cc in 1:ncol(x))
        do.call("text", c(cc+0.5, 0.5, colnames(x)[cc], dimNameOpts))    
}
