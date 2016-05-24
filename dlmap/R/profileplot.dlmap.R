`profileplot.dlmap` <- 
function(object, chr, marker.names=TRUE, QTLpos=TRUE, pch=20, ...)
{
   if (missing(chr)) chr <- c(1:length(object$profile))

   if (is.null(object$profile)) {
	cat("No QTL detected, cannot plot profile")
	return(invisible())
   }

   if (attr(object$input, "type")=="other") {
	cat("Cannot plot profiles for association analysis")
	return(invisible())
   }

   if (class(chr)=="character") chr <- match(chr, names(object$profile))

   nplots <- length(chr)

   if (nplots>1)
     op <- par(mfrow=c(ceiling(nplots/2), 2), font.lab=2)

   for (ii in chr)
   {
     prof <- object$profile[[ii]]
     ichr <- names(object$profile)[ii]
     plot(prof, xlab=paste("Chr ", ichr, " Position (cM)", sep=""),
	 ylab="Wald", pch=pch, type="o")

     mrk <- setdiff(1:ncol(prof), grep("loc", names(object$mapp[[ichr]])))
     labels <- names(object$input$mapp[[ichr]])[mrk]
     pos <- prof[1, mrk]

     if (marker.names)
        axis(side=3, at=pos, labels=labels)

     qpos <- as.numeric(as.character(object$Summary$Pos))
     qchr <- which(as.character(object$Summary$Chr)==ichr)

     if (QTLpos)
	abline(v=qpos[qchr])
   }

}
 
