##
##  PURPOSE:   Plot individual longitudinal profiles
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   23/03/2010 (as a stand alone function)
##             26/11/2010:  added to the mixAK package
##             21/06/2013:  highlighting of selected profiles implemented
##
##  FUNCTIONS: plotProfiles
##
## ==========================================================================

## *************************************************************
## plotProfiles
## *************************************************************
##
plotProfiles <- function(ip, data, var, trans, tvar, gvar,
                         auto.layout=TRUE, lines=TRUE, points=FALSE, add=FALSE,
                         xlab="Time", ylab, xaxt="s", yaxt="s", xlim, ylim, main,
                         lcol, col, bg, lty=1, lwd=1, pch=21,
                         cex.points=1,
                         highlight, lines.highlight=TRUE, points.highlight=TRUE,
                         lcol.highlight="red3", col.highlight="red3", bg.highlight="orange",
                         lty.highlight=1, lwd.highlight=2, pch.highlight=23,
                         cex.highlight=1)
{   
  if (lines){
    if (missing(col)) col <- rainbow_hcl(1, start = 230, c = 40)
  }

  if (points){
    if (missing(col)) col <- "darkblue"    
  }
  if (missing(lcol)) lcol <- col  
  if (missing(bg)) bg <- rainbow_hcl(1, start = 230, c = 40)
  
  if (missing(xlim)){
    xlim <- range(data[, tvar], na.rm=TRUE)
    KEEP <- rep(TRUE, nrow(data))
  }else{
    KEEP <- data[, tvar] >= xlim[1] & data[, tvar] <= xlim[2]
    xlim <- range(data[KEEP, tvar], na.rm=TRUE)
  }  
  if (missing(ylim)){
    if (missing(trans)) ylim <- range(data[KEEP, var], na.rm=TRUE)
    else                ylim <- range(trans(data[KEEP, var]), na.rm=TRUE) 
  }  
  if (missing(ylab)){
    if (missing(trans)) ylab <- substitute(var)
    else                ylab <- paste(substitute(trans), "(", var, ")", sep="")    
  }  

  if (!missing(gvar)){
    GROUP <- levels(data[, gvar])
    if (length(col) == 1){
      if (length(bg) == 1) col <- rainbow_hcl(length(GROUP)) else col <- rep(col, length(GROUP))     
    }
    if (length(lty) == 1) lty <- rep(lty, length(GROUP))
    if (length(lcol) == 1) lcol <- rainbow_hcl(length(GROUP))
    if (length(col) != length(GROUP)) stop("incorrect col supplied")
    if (length(lcol) != length(GROUP)) stop("incorrect lcol supplied")    
    if (length(bg) == 1) bg <- rainbow_hcl(length(GROUP))
    if (length(bg) != length(GROUP)) stop("incorrect bg supplied")    
    names(col) <- names(bg) <- names(lcol) <- names(lty) <- GROUP
  }else{
    col <- col[1]
    lcol <- lcol[1]
    lty <- lty[1]
  }  
  
  if (auto.layout & !add){
    oldPar <- par(mfrow=c(1, 1), bty="n")
    on.exit(oldPar)
  }

  if (!add) plot(xlim, ylim, type="n", xaxt=xaxt, yaxt=yaxt, xlab=xlab, ylab=ylab)
  for (i in 1:length(ip)){
    if (!missing(gvar)){
      COL <- col[ip[[i]][1, gvar]]
      BG <- bg[ip[[i]][1, gvar]]
      LCOL <- lcol[ip[[i]][1, gvar]]
      LTY <- lty[ip[[i]][1, gvar]]      
    }else{
      COL <- col[1]
      BG <- bg[1]
      LCOL <- lcol[1]
      LTY <- lty[1]
    }  
    if (missing(trans)){
      if (lines) lines(ip[[i]][, tvar], ip[[i]][, var], col=LCOL, lty=LTY, lwd=lwd)
      if (points) points(ip[[i]][, tvar], ip[[i]][, var], col=COL, bg=BG, pch=pch, cex=cex.points)
    }else{  
      if (lines) lines(ip[[i]][, tvar], trans(ip[[i]][, var]), col=LCOL, lty=LTY, lwd=lwd)
      if (points) points(ip[[i]][, tvar], trans(ip[[i]][, var]), col=COL, bg=BG, pch=pch, cex=cex.points)      
    }  
  }
  if (!missing(main)) title(main=main)

  if (!missing(highlight)){
    for (hh in highlight){
      if (hh < 1 | hh > length(ip)){
        warning("profile number ", hh, " not highlighted.", sep = "")
        next
      }
      if (missing(trans)){
        if (lines.highlight)  lines(ip[[hh]][, tvar], ip[[hh]][, var],  col=lcol.highlight, lty=lty.highlight, lwd=lwd.highlight)
        if (points.highlight) points(ip[[hh]][, tvar], ip[[hh]][, var], col=col.highlight, bg=bg.highlight, pch=pch.highlight, cex=cex.highlight)
      }else{
        if (lines.highlight)  lines(ip[[hh]][, tvar], trans(ip[[hh]][, var]),  col=lcol.highlight, lty=lty.highlight, lwd=lwd.highlight)
        if (points.highlight) points(ip[[hh]][, tvar], trans(ip[[hh]][, var]), col=col.highlight, bg=bg.highlight, pch=pch.highlight, cex=cex.highlight)
      }  
    }  
  }  
  
  return(invisible(ip))
}  
