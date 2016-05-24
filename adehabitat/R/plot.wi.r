"plot.wi" <- function(x, caxis=0.7, clab=1, ylog=FALSE, errbar=c("CI", "SE"),
                      main="Manly selectivity measure", noorder = TRUE, ...)
{
    ## Verifications
    errbar<-match.arg(errbar)
    opar<-par(ask=TRUE)
    on.exit(par(opar))
    if (!inherits(x, "wi"))
      stop("x should be of class wi")
    eb<-ifelse(errbar=="SE", 1, abs(qnorm(x$alpha/length(x$wi))) )

    ## Should the selection ratios be ordered?
    if (noorder)
        wi<-sort(x$wi, decreasing=TRUE)
    else
        wi<-x$wi

    ## A log scale?
    if ((any(wi==0))&(ylog)) {
      warning("zero values in x, ylog has been set to FALSE")
      ylog<-FALSE
    }
    logy<-ifelse(ylog, "y", "")


    ## Also order the standard error
    if (noorder)
      sewi<-x$se.wi[order(x$wi, decreasing=TRUE)]
    else
      sewi<-x$se.wi

    ## 0 standard errors in case of NA
    sewi[is.na(sewi)]<-0
    nwi<-names(wi)

    ## Computes the range
    rgy<-range(c(wi, wi+eb*sewi, wi-eb*sewi))

    ## The legend
    textleg<-paste("Selection ratios (+/-", errbar,")")
    if (inherits(x, "wiII")|inherits(x,"wiIII"))
      textleg<-paste("Global Selection ratios (+/-", errbar,")")

    ## Draws the plot
    if (!ylog) rgy[1]<-0
    plot(wi, axes=FALSE, ylim=rgy, ty="n", xlab="",
         ylab=textleg,
         cex.lab=clab, log=logy, main=main, ...)
    axis(side=1, at=c(1:length(wi)), labels=names(wi), cex.axis=caxis, las=2)
    axis(side=2, cex.axis=caxis)
    box()
    points(c(1:length(wi)), wi, pch=16)
    lines(1:length(wi), wi)
    abline(h=1, lwd=2)
    for (i in 1:length(wi)) {
        lines(c(i,i), c(wi[i]-eb*sewi[i], wi[i]+eb*sewi[i]))
        lines(c(i-0.1,i+0.1), c(wi[i]-eb*sewi[i], wi[i]-eb*sewi[i]))
        lines(c(i-0.1,i+0.1), c(wi[i]+eb*sewi[i], wi[i]+eb*sewi[i]))
    }


    ## plots the standardized selection ratios on designs I
    if (inherits(x, "wiI")) {
        if (noorder)
            Bi<-x$Bi[order(x$wi, decreasing=TRUE)]
        else
            Bi<-x$Bi

        ## Draws the plot
        plot(Bi, axes=FALSE, ty="n", xlab="",
             cex.lab=clab, main="Scaled selection ratios", ...)
        axis(side=1, at=c(1:length(wi)), labels=names(wi),
             cex.axis=caxis, las=2)
        axis(side=2, cex.axis=caxis)
        lines(1:length(wi), Bi)
        points(c(1:length(wi)), Bi, pch=16)
        box()

        ## plots the use and availability
        if (noorder) {
            ut<-x$used.prop[order(x$wi, decreasing=TRUE)]
            seu<-x$se.used[order(x$wi, decreasing=TRUE)]
            sea<-x$se.avail[order(x$wi, decreasing=TRUE)]
            av<-x$avail.prop[order(x$wi, decreasing=TRUE)]
        } else {
            ut<-x$used.prop
            seu<-x$se.used
            sea<-x$se.avail
            av<-x$avail.prop
        }

        ## Plot them
        rgy<-range(c(av, ut-eb*seu, ut+eb*seu, av-eb*sea, av+eb*sea))
        rgy<-c(rgy[1], rgy[2]+(rgy[2]-rgy[1])/4)
        plot(ut, axes=FALSE, ty="n", xlab="", cex.lab=clab, ylim=rgy,
             main="Used and available proportions",
             ylab=paste("Porportion (+/-", errbar,")"),...)
        points(1:length(wi)-0.05, av, pch=16)
        points(1:length(wi)+0.05, ut, pch=2)

        for (i in 1:length(wi)) {
            lines(c(i,i)+0.05, c(ut[i]-eb*seu[i], ut[i]+eb*seu[i]))
            lines(c(i-0.02,i+0.02)+0.05, c(ut[i]-eb*seu[i], ut[i]-eb*seu[i]))
            lines(c(i-0.02,i+0.02)+0.05, c(ut[i]+eb*seu[i], ut[i]+eb*seu[i]))
        }

        ## If the availability is not known
        if (!x$avknown) {
            for (i in 1:length(wi)) {
                lines(c(i,i)-0.05, c(av[i]-eb*sea[i], av[i]+eb*sea[i]))
                lines(c(i-0.02,i+0.02)-0.05,
                      c(av[i]-eb*sea[i], av[i]-eb*sea[i]))
                lines(c(i-0.02,i+0.02)-0.05,
                      c(av[i]+eb*sea[i], av[i]+eb*sea[i]))
            }
        }
        axis(side=1, at=c(1:length(wi)), labels=names(wi),
             cex.axis=caxis, las=2)
        axis(side=2, cex.axis=caxis)
        box()
        legend(1,rgy[2], c("Available", "Used"), pch=c(16,2), cex=clab)
    } else {
        ## For the designs II and III
        if (noorder)
            wij<-x$wij[,order(x$wi, decreasing=TRUE)]
        else
            wij<-x$wij
        iii<-as.vector(wij)
        rgy<-range(iii[!is.na(iii)])
        plot(1, ty="n", ylim=rgy, xlim=c(1,ncol(wij)), xlab="",
             ylab=paste("Selection ratios"),
             cex.lab=clab, log=logy, axes=FALSE,
             main=main, ...)
        axis(side=1, at=c(1:length(wi)), labels=names(wi),
             cex.axis=caxis, las=2)
        axis(side=2, cex.axis=caxis)
        box()
        pt<-seq(-0.1, 0.1, by=0.2/nrow(wij))

        for (j in 1:nrow(wij)) {
            points(c(1:length(wi)), wij[j,], pch=16, col=j)
            lines(1:length(wi), wij[j,], col=j)
            abline(h=1, lwd=2)
        }
        rgx<-ncol(wij)/5
        legend(ncol(wij)-rgx, rgy[1]+19*(rgy[2]-rgy[1])/20,
               legend=row.names(wij), pch=16, col=1:nrow(wij),
               lwd=1, cex=clab)
    }
  }

