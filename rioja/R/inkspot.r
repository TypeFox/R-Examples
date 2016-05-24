inkspot <- function(data, gradient=1:nrow(data), use.rank=FALSE, reorder.species = TRUE, x.axis=c("sites", "gradient", "none"), x.axis.top=FALSE, site.names=NULL, spec.names=NULL, pch=21, cex.max=3, col="black", bg="darkgrey", legend.values=c(2, 5, 10, 20, 50), ...) {
   x.axis = match.arg(x.axis)
   ord <- order(gradient)
   grad.srt <- sort(gradient)
   if (reorder.species) {
#     wa.sc <- apply(data[ord, ], 2, function(x, env) { sum(x*env, na.rm=TRUE) / sum(x, na.rm=TRUE) }, env=1:nrow(data))
     wa.sc <- apply(data, 2, function(x, env) { sum(x*env, na.rm=TRUE) / sum(x, na.rm=TRUE) }, env=gradient)
     spec.ord <- order(wa.sc)
   }
   else {
      spec.ord <- 1:ncol(data)
   }
   if (use.rank)
      grad.srt <- 1:nrow(data)
   nR <- nrow(data)
   nC <- ncol(data)
   ss <- list(sites=ord, spec=spec.ord)
   ddd <- as.vector(as.matrix(sqrt(data[ss$sites, ss$spec])))
 	 ra <- max(ddd, na.rm=TRUE)
   ddd <- ddd / ra * cex.max
   r <- rep((1:nC), each=nR)
   c <- rep(grad.srt, times=nC)
   if (!is.null(site.names))
      sn <- site.names[ss$sites]
   else
      sn <- rownames(data)[ss$sites]
   if (!is.null(spec.names))
      spn <- spec.names[ss$spec]
   else
      spn <- colnames(data)[ss$spec]
   plot(c, r, cex=ddd, pch=pch, yaxt="n", xaxt="n", ylab="", xlab="", col=col, bg=bg, ...)
#   axis(side=1, at=1:nR, labels=sn, las=2, ...)
   if (x.axis=="sites")
     axis(side=1, at=grad.srt, labels=sn, las=2, ...)
   else {
     if (x.axis=="gradient")
       axis(side=1, col="black", ...)
     else
       axis(side=1, at=grad.srt, labels=rep("", length(sn), ...))
   }
   axis(side=2, at=1:nC, labels=spn, las=1, ...)
   if (x.axis.top) {
      if (use.rank) {
         x <- pretty(gradient)
         breaks <- apply(data.frame(x), 1, function(x, y) { z <- which(y < x); if(length(z)==0) z=0; max(z, na.rm=TRUE); }, y=sort(gradient)) 
         breaks[breaks==0] <- NA
         axis(side=3, at = breaks, labels=x, ...)
      }
      else {
         axis(side=3, ...)
      }
   }
   if (!is.null(legend.values))
     legend("topleft", as.character(legend.values), pch=pch, bg="white", cex=.8, pt.cex=sqrt(legend.values) / ra * cex.max, col=col, pt.bg=bg)
   invisible(ss)
}

