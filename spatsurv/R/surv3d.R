##' Do a 3d plot of spatial survival data
##'
##' Uses rgl graphics to make a spinny zoomy plot
##' @title Spatial Survival Plot in 3D
##' @param spp A spatial points data frame
##' @param ss A Surv object (with right-censoring)
##' @param lwd Line width for stems
##' @param lcol Line colour for stems
##' @param lalpha Opacity for stems
##' @param pstyle Point style "point" or "text"
##' @param psize Vector of length 2 for uncensored/censored points size
##' @param pcol Vector of length 2 for uncensored/censored points colours
##' @param ptext Vector of length 2 for uncensored/censored text characters
##' @param palpha Opacity for points/text
##' @param title Main title for plot
##' @param basegrid add a grid at t=0
##' @param baseplane add a plane at t=0
##' @return nothing
##' @examples
##' \dontrun{
##' require(sp)
##' require(survival)
##' d = data.frame(
##'   x=runif(40)*1.5,
##'   y = runif(40),
##'   age=as.integer(20+30*runif(40)),
##'   sex = sample(c("M","F"),40,TRUE)
##' )
##' coordinates(d)=~x+y
##' d$surv = Surv(as.integer(5+20*runif(40)),runif(40)>.9)
##' clear3d();surv3d(d,d$surv,baseplane=TRUE,basegrid=TRUE)
##' clear3d();surv3d(d,d$surv,baseplane=TRUE,basegrid=TRUE,pstyle="t",lalpha=0.5,lwd=3,palpha=1)
##' }
##' @author Barry S Rowlingson
##' @export
surv3d <- function(spp, ss,
                   lwd=2, lcol="black",lalpha=1.0,
                   pstyle=c("point","text"),
                   psize=c(20,10),
                   pcol=c("red","black"),
                   ptext = c("X",""),
                   palpha=1.0,
                   title="Spatial Survival",
                   basegrid=TRUE, baseplane=TRUE){

    pstyle=match.arg(pstyle)
    
    nr = nrow(spp)
    xy = rbind(coordinates(spp),coordinates(spp))
    xyz = cbind(xy,c(rep(0,nr), ss[,"time"]))

    # weave the lines
    xyz = xyz[rep(1:nr,rep(2,nr))+rep(c(0,nr),nr),]

    # segments3d takes pairs for line segments
    rgl::segments3d(xyz,lwd=lwd, col=lcol, alpha=lalpha)

    # add points for uncensored obs
    unc = ss[,"status"] == 1
    xyp = cbind(coordinates(spp),ss[,"time"])

    if(pstyle=="text"){
        psize=psize/12
        rgl::text3d(xyp[unc,,drop=FALSE],texts=ptext[1],col=pcol[1], cex=psize[1], alpha=palpha)
        rgl::text3d(xyp[!unc,,drop=FALSE],texts=ptext[2],col=pcol[2], cex=psize[2], alpha=palpha)
    }
    if(pstyle=="point"){
        rgl::points3d(xyp[unc,,drop=FALSE], size=psize[1], col=pcol[1], alpha=palpha)
        rgl::points3d(xyp[!unc,,drop=FALSE], size=psize[2], col=pcol[2], alpha=palpha)
    }
    rgl::aspect3d(c(1,1,1))
    rgl::title3d(main=title,xlab='x',ylab='y',zlab='Time')
    rgl::axes3d()

    ## minimum visible z coord (zero gets clipped)
    zminvis = min(ss[,"time"])/100    

    if (basegrid){
        x=seq(min(xyp[,1]),max(xyp[,1]),len=20)
        y=seq(min(xyp[,2]),max(xyp[,2]),len=20)
        rgl::abclines3d(x,min(y),zminvis, a=0, b= 1, c=0, col="gray")
        rgl::abclines3d(min(x),y,zminvis, a=1, b= 0, c=0, col="gray")
    }
    if(baseplane){
        rgl::planes3d(0,0,1,-zminvis,col="gray",alpha=0.5)
    }
    invisible(0)
}
