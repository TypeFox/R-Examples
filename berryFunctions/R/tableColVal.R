#' Table with values with value-dependent colored backgrounds in pdf
#' 
#' Table with numbers and corresponding color in the background of each cell.
#' (heatmap)
#' 
#' @details I saw a presentation today with a table of values of differences between several
#' models and datasets of global precipitation. I decided I don't like reading
#' 20+ values, and would like to see a corresponding color in the background of each cell. (heatmap)
#' Writing this function took me about 1 hour and 30 minutes and was a nice coding excercise.
#' Feedback welcome at berry-b@gmx.de!
#' 
#' @return None. PDF or plot produced.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Nov. 2012
#' @seealso \code{\link{pdf}}, \code{\link{heatmap}}
#' @keywords hplot
#' @export
#' @examples
#' 
#' Bsp <- matrix(c(21,23,26,27, 18,24,25,28, 14,17,23,23, 16,19,21,25), ncol=4, byrow=TRUE)
#' colnames(Bsp) <- paste0("Measure", LETTERS[1:4])
#' rownames(Bsp) <- paste("prod", 8:11, sep="_")
#' Bsp
#' 
#' tableColVal(Bsp)
#' tableColVal(Bsp, nameswidth=0.1) # relative to plot width
#' tableColVal(Bsp, namesheight=0.5, argcol=list(srt=90))
#' 
#' tableColVal(Bsp, argrow=list(col="red", cex=2) )
#' tableColVal(Bsp, Range=c(10,40))
#' tableColVal(Bsp, Range=c(20,40))
#' tableColVal(Bsp, palette=heat.colors(12))
#' tableColVal(Bsp, palette=c(2,4,7), argmain=list(labels="last\ncomparison"))
#' 
#' \dontrun{
#' ## Rcmd check --as-cran doesn't like to open external devices such as pdf,
#' ## so this example is excluded from running in the checks.
#' tableColVal(Bsp, pdf=TRUE, width=12)  # further arguments to pdf possible.
#' 
#' Bsp2 <- matrix(sample(1:100, 30), ncol=6, byrow=TRUE)
#' graphics.off(); X11(height=4)
#' tableColVal(Bsp2)
#' }
#' 
#' @param mat Matrix with values
#' @param pdffile Name of File to write to. DEFAULT: "table_col_val.pdf"
#' @param pdf Should table be written to pdffile? (Else it will plot in x11). DEFAULT: FALSE
#' @param nameswidth Relative width of row names at the left, as a percentage of plot. DEFAULT: 0.3
#' @param namesheight Relative height of column names at the top. DEFAULT: 0.1
#' @param palette Color palette for the heatmap. DEFAULT: seqPal(nrow(mat)*ncol(mat))
#' @param Range Range of values mapped linearly to color palette. DEFAULT: range(mat,finite
#' @param argclass List of arguments specifying how to call \code{\link{classify}}, eg. method. DEFAULT: NULL
#' @param argrow,argcol,argcell,argmain List of arguments passed to \code{\link{text}} 
#'        in row and column names, cell content and topleft cell, respectively. Could be cex, col, srt, etc. DEFAULTS: NULL
#' @param \dots Further arguments passed to \code{\link{pdf}}.
#' 
tableColVal <- function(
mat,
pdffile="table_col_val.pdf",
pdf=!missing(pdffile),
nameswidth=0.3,
namesheight=0.1,
palette=seqPal(nrow(mat)*ncol(mat)),
Range=range(mat,finite=TRUE),
argclass=NULL,
argrow=NULL,
argcol=NULL,
argcell=NULL,
argmain=NULL,
...)
{
# expand pdf-path to working directory if only file name (without path) is given:
if(pdf){   if(!grepl("/", pdffile) | grepl("\\", pdffile, fixed=TRUE) )
               pdffile <- paste(getwd(), pdffile, sep="/")  #"
        pdf(pdffile, ...) }# open pdf device
mat <- as.matrix(mat)
nc <- ncol(mat) ; nr <- nrow(mat)
# set plot
op <- par(mai=c(0, 0, namesheight*par()$pin[2], 0), xpd=TRUE )
plot(1, ylim=c(nr+1, 1), xlim=c(0,1), type="n", xaxs="i", yaxs="i", axes=FALSE, ann=FALSE)
# set positions for text and lines
rights <- seq(nameswidth, 1, len=nc+1)
lefts <- c(0, rights[1:nc] )
middles <- nameswidth + (1:nc*2-1) * (1-nameswidth)/nc/2
# define color for each value of mat
cl <- do.call(classify, args=owa(list(x=mat, breaks=length(palette), Range=Range), argclass))
##mod <- lm(c(1, length(palette)) ~ Range)[[1]]
##lincol <- round(as.vector(mat) * mod[2] + mod[1])
# plot rectancles with colors corresponding to values of mat
rect(xleft=rep(lefts[-1], each=nr), xright=rep(rights[-1], each=nr),
     ybottom=rep(2:(nr+1), nc), ytop=rep(1:nr, nc), col=palette[cl$index] , border=NA)
abline(v=rights, h=1:nr)
# add "titles"
ytitles <- 1-(namesheight*nr/2)
do.call(text, args=owa(d=list(x=middles,      y=ytitles,  labels=colnames(mat)), argcol))
do.call(text, args=owa(d=list(x=nameswidth/2, y=ytitles,  labels="tableColVal"), argmain))
do.call(text, args=owa(d=list(x=nameswidth/2, y=1:nr+0.5, labels=rownames(mat)), argrow))
# add text to each cell
do.call(text, args=owa(d=list(x=rep(middles, each=nr), y=rep(1:nr, nc)+0.5, 
                              labels=as.vector(mat)), argcell))
# Set old paramaters again:
par(op)
if(pdf) { dev.off() ; message("PDF-File is located here:", pdffile, "\n") }# close pdf device
} # end of function
