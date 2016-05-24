"kplot.kselect" <- function (object, xax = 1, yax = 2, csub = 2,
                             possub = c("topleft", "bottomleft",
                             "bottomright", "topright"),
                             addval=TRUE, cpoint=1, csize=1, clegend=2, ...)
{
    ## Verifications
    possub<-match.arg(possub)
    x<-object
    if (!inherits(x, "kselect"))
        stop("x should be a 'kselect' object")
    if (x$nf == 1) {
        hist.kselect(x)
        return(invisible())
    }

    ## Coordinates of the available points on the axes of the K-select
    ## (scores of the rows of initab)
    Xi<-x$initab
    Xrecalc<-t(as.matrix(apply(Xi, 1,
                               function(y) y*x$lw/sum(x$lw))))%*%as.matrix(x$l1)
    ## Preparation of the data needed for this graphs (spliting the "global"
    ## table into a list of multiple tables, idem for the utilization weights)
    rx<-range(Xrecalc[,xax])
    ry<-range(Xrecalc[,yax])
    li.Xi<-split(as.data.frame(Xrecalc), x$initfac)
    li.wei<-split(x$initwei, x$initfac)
    li.wei<-lapply(li.wei, function(x) x/sum(x) )
    maxsqrtw<-max(sqrt(unlist(li.wei)))

    ## The graphs
    csi<-0
    for (i in 1:length(li.wei))
        csi[i]<-csize*max(sqrt(li.wei[[i]]))/maxsqrtw
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    ngraph<-length(li.Xi)
    if (addval) {
        par(mfrow = n2mfrow(ngraph+1)) ## +1 for the legend
    } else {
        par(mfrow = n2mfrow(ngraph))
    }

    ## One graph per animal
    for (i in 1:ngraph) {
        Xtmp<-li.Xi[[i]]
        wgtmp<-li.wei[[i]]
        if (addval) {
            s.value(Xtmp, wgtmp, xax, yax,
                    sub=names(li.Xi)[i], cpoint=cpoint, xlim=rx,
                    ylim=ry, clegend=0,
                    csub=1.5, cgrid=1.5, csize=csi[i])
        }
        s.distri(Xtmp, wgtmp, xax, yax,
                 sub=names(li.Xi)[i], add.plot=addval,
                 cpoint=cpoint, xlim=rx, ylim=ry,
                 ...)
    }

    ## adds a legend if addval is TRUE
    if (addval) {
        coo <- scatterutil.base(dfxy = Xtmp, xax = xax, yax = yax,
                                xlim = rx, ylim = ry, grid = FALSE,
                                addaxes = FALSE,
                                cgrid = 0, include.origin = FALSE,
                                origin = c(0,0),
                                sub = "", csub = 0,
                                possub = "bottomleft", pixmap = NULL,
                                contour = NULL, area = NULL, add.plot = FALSE)

        coeff <- diff(range(coo$x))/15
        br0<-pretty(unlist(li.wei), 4)
        l0 <- length(br0)
        br0 <- (br0[1:(l0 - 1)] + br0[2:l0])/2
        sq0 <- sqrt(abs(br0))
        sq0 <- csize * coeff * sq0/max(sqrt(abs(wgtmp)))
        sig0 <- sign(br0)
        scatterutil.legend.bw.square(pretty(unlist(li.wei), 4),
                                     sq0, sig0, clegend=clegend)
    }
}
