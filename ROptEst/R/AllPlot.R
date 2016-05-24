setMethod("plot", signature(x = "IC", y = "missing"),
    function(x, ...,withSweave = getdistrOption("withSweave"),
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3],
             with.legend = FALSE, legend = NULL, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             withMBR = FALSE, MBRB = NA, MBR.fac = 2, col.MBR = par("col"),
             lty.MBR = "dashed", lwd.MBR = 0.8,  n.MBR = 10000,
             scaleX = FALSE, scaleX.fct, scaleX.inv,
             scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
             scaleN = 9, x.ticks = NULL, y.ticks = NULL,
             mfColRow = TRUE, to.draw.arg = NULL){

        mcl <- match.call(call = sys.call(sys.parent(1)), expand.dots = TRUE)

        L2Fam <- eval(x@CallL2Fam); trafO <- trafo(L2Fam@param)
        dims  <- nrow(trafO); to.draw <- 1:dims
        dimnms  <- c(rownames(trafO))
        if(is.null(dimnms))
           dimnms <- paste("dim",1:dims,sep="")
        if(! is.null(to.draw.arg)){
            if(is.character(to.draw.arg))
                 to.draw <- pmatch(to.draw.arg, dimnms)
            else if(is.numeric(to.draw.arg))
                 to.draw <- to.draw.arg
        }
        dims0 <- length(to.draw)

        MBRB <- matrix(rep(t(MBRB), length.out=dims0*2),ncol=2, byrow=T)
        if(withMBR && all(is.na(MBRB))){
           robModel <- InfRobModel(center = L2Fam, neighbor =
                             ContNeighborhood(radius = 0.5))
           ICmbr <- try(optIC(model = robModel, risk = asBias()), silent=TRUE)
           if(!is(ICmbr,"try-error"))
              MBRB <- .getExtremeCoordIC(ICmbr, distribution(L2Fam), to.draw,
                              n = n.MBR)
           else withMBR <- FALSE
        }
        mcl$MBRB <- MBRB
        mcl$withMBR <- withMBR
        plm <- getMethod("plot", signature(x = "IC", y = "missing"),
                           where="RobAStBase")
        do.call(plm, as.list(mcl[-1]), envir=parent.frame(2))
       return(invisible())
      })

.getExtremeCoordIC <- function(IC, D, indi, n = 10000){
    x <- q(D)(seq(1/2/n,1-1/2/n, length=n))
    y <- (matrix(evalIC(IC,matrix(x,ncol=1)),ncol=n))[indi,]
    return(cbind(min=apply(y,1,min),max=apply(y,1,max)))
}