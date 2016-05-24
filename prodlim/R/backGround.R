#' Background and grid color control.
#' 
#' Some users like background colors, and it may be helpful to have grid lines
#' to read off e.g. probabilities from a Kaplan-Meier graph. Both things can be
#' controlled with this function. However, it mainly serves
#' \code{\link{plot.prodlim}}.
#' 
#' 
#' @param xlim Limits for the xaxis, defaults to par("usr")[1:2].
#' @param ylim Limits for the yaxis, defaults to par("usr")[3:4].
#' @param bg Background color. Can be multiple colors which are then switched
#' at each horizontal line.
#' @param fg Grid line color.
#' @param horizontal Numerical values at which horizontal grid lines are
#' plotted.
#' @param vertical Numerical values at which vertical grid lines are plotted.
#' @param border The color of the border around the background.
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @keywords survival
#' @examples
#' 
#' 
#' plot(0,0)
#' backGround(bg="beige",fg="red",vertical=0,horizontal=0)
#' 
#' plot(0,0)
#' backGround(bg=c("yellow","green"),fg="red",xlim=c(-1,1),ylim=c(-1,1),horizontal=seq(0,1,.1))
#' backGround(bg=c("yellow","green"),fg="red",horizontal=seq(0,1,.1))
#' 
#' @export
backGround <- function(xlim,
                       ylim,
                       bg="white",
                       fg="gray77",
                       horizontal=NULL,
                       vertical=NULL,
                       border="black"){
    U <- par("usr")
    if (missing(xlim))
        xlim <- c(U[1],U[2])
    if (missing(ylim))
        ylim <- c(U[3],U[4])
    # background
    if (!is.null(bg)){
        if (length(bg)==1){
            rect(U[1],U[3],U[2],U[4],col=bg[1], border=border)
        }else{
             if (length(bg)>1){
                 ybot <- sort(unique(c(ylim[1],horizontal,ylim[2])))
                 NR <- length(ybot)
                 bcol <- rep(bg,length.out=NR)
                 nix <- sapply(1:(NR-1),function(r){
                                     ## for (r in 1:(NR-1)){
                                     ## rect(xleft=xlim[1],xright=xlim[2],ybottom=ybot[r],ytop=ybot[r+1],col=bcol[r],border=FALSE)
                                     ## polygon(x=c(xlim[1],xlim[1],xlim[2],xlim[2],xlim[1]),
                                     polygon(x=c(U[1],U[1],U[2],U[2],U[1]),
                                             y=c(ybot[r],ybot[r+1],ybot[r+1],ybot[r],ybot[r]),
                                             col=bcol[r],
                                             border=FALSE)
                                     ## do NOT specify: density=100 as this slows this down!
                                 })
             }
         }
    }
    # grid 
    if (length(fg)>0){
        if (length(vertical)>0)
            abline(v=vertical,col=fg)
        if (length(horizontal)>0)
            abline(h=horizontal,col=fg)
    }
}
