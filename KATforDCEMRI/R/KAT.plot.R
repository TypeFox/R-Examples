##    KATforDCEMRI: a Kinetic Analysis Tool for DCE-MRI
##    Copyright 2014 Genentech, Inc.
##
##    For questions or comments, please contact
##    Gregory Z. Ferl, Ph.D.
##    Genentech, Inc.
##    Development Sciences
##    1 DNA Way, Mail stop 463A
##    South San Francisco, CA, United States of America
##    E-mail: ferl.gregory@gene.com

KAT.plot <-
function(F1=0, F2=0, F3=0, F4=0, F5=0, F6=0, F7=0, F8=0, plot.param="Ktrans", range.map = 1.5, cutoff.map=0.85, ...){


    dcemri.data <- NULL
    modeltype1 <- "xTofts"
    modeltype2 <- "Tofts"

    files <- c(F1, F2, F3, F4, F5, F6, F7, F8)

    MAP.list <- vector("list", length(files))
    AIC.list <- vector("list", length(files))
    maptimes.list <- vector("list", length(files))
    timestatus.list <- vector("list", length(files))
    aif.list <- vector("list", length(files))
    cvaif.list <- vector("list", length(files))
    aifparamsfitted.list <- vector("list", length(files))
    aiffitted.list <- vector("list", length(files))
    paramsxTofts.list <- vector("list", length(files))
    cvroi.list <- vector("list", length(files))
    ccmedian.list <- vector("list", length(files))
    roimedianfittedTofts.list <- vector("list", length(files))
    roimedianfittedxTofts.list <- vector("list", length(files))
    aifparamsfitted.list <- vector("list", length(files))
    nx.vec <- vector("numeric", length(files))
    ny.vec <- vector("numeric", length(files))
    nt.vec <- vector("numeric", length(files))
    slice.vec <- vector("numeric", length(files))
    modelstructure.vec <- vector("numeric", length(files))
    IDvisit.vec <- vector("character", length(files))
    aifsource.vec <- vector("character", length(files))


    xmin <- 10000
    xmax <- 0
    ymin <- 10000
    ymax <- 0
    MAPul <- 0

    number.of.visits <- 0

    cat("--------", "\n")

    for(i in 1:length(files)){
        if(files[i]!=0){
            cat("loading file for panel ", i, ": ", files[i], sep="", "\n")
            cat("--------", "\n")

            load(files[i])

            number.of.visits <- number.of.visits + 1

            if(plot.param=="Ktrans")
                MAP.list[[i]] <- dcemri.data$mapKtransxT

            if(plot.param=="kep")
                MAP.list[[i]] <- dcemri.data$mapkepxT

            if(plot.param=="vb")
                MAP.list[[i]] <- dcemri.data$mapvbxT

            if(plot.param=="ve")
                MAP.list[[i]] <- dcemri.data$mapvexT

            AIC.list[[i]] <- dcemri.data$mapAICcompare

            modelstructure.vec[i] <- 2
            nx.vec[i] <- dcemri.data$nx
            ny.vec[i] <- dcemri.data$ny
            nt.vec[i] <- dcemri.data$nt

            slice.vec[i] <- dcemri.data$args$slice

            IDvisit.vec[i] <- strsplit(dcemri.data$args$IDvisit, split=".RData")[[1]][1]

            aifsource.vec[i] <- "vec"

            aifparamsfitted.list[[i]] <- dcemri.data$aifparamsfitted
            maptimes.list[[i]] <- dcemri.data$maptimes
            timestatus.list[[i]] <- dcemri.data$timeStatus
            aif.list[[i]] <- dcemri.data$aif
            cvaif.list[[i]] <- dcemri.data$cvaif
            aifparamsfitted.list[[i]] <- dcemri.data$aifparamsfitted
            aiffitted.list[[i]] <- dcemri.data$aiffitted
            paramsxTofts.list[[i]] <- dcemri.data$paramestwholeroixTofts
            cvroi.list[[i]] <- dcemri.data$cvwholeroixTofts
            ccmedian.list[[i]] <- dcemri.data$ccmedian
            roimedianfittedxTofts.list[[i]] <- dcemri.data$roimedianfittedxTofts
            if(modelstructure.vec[i]==2)
                roimedianfittedTofts.list[[i]] <- dcemri.data$roimedianfittedTofts

            xmin <- min(xmin, dcemri.data$roiplotparams$xmin, na.rm=TRUE)
            xmax <- max(xmax, dcemri.data$roiplotparams$xmax, na.rm=TRUE)
            ymin <- min(ymin, dcemri.data$roiplotparams$ymin, na.rm=TRUE)
            ymax <- max(ymax, dcemri.data$roiplotparams$ymax, na.rm=TRUE)
            MAPul <- max(MAPul, dcemri.data$roiplotparams$MAPul, na.rm=TRUE)

            rm(dcemri.data)
        }
    }


    for(i in 1:length(files)){
        if(files[i]!=0)
            MAP.list[[i]][MAP.list[[i]]>=MAPul*0.99] <- MAPul*0.99
    }

    legend <- seq(0,MAPul,by=0.001)
    dim(legend) <- c(1,length(legend))

    if(number.of.visits <= 8){

        cat("generating PDF...", sep="")

        width.coeff <- number.of.visits+1

        pdf(file=paste(IDvisit.vec[1], "-multislice.pdf", sep=""), height=14, width=4*width.coeff)
        layout(matrix(c(1:(3*width.coeff)), 3, width.coeff, byrow=TRUE), widths=c(1,rep(4, number.of.visits),1,rep(4, number.of.visits),1,rep(4, number.of.visits)))
        par(omi=c(0.25, 0.25, 0.25, 0.25))

        image(legend, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), zlim=c(0,MAPul), xaxt="n", yaxt="n", xlab="0", main=format(MAPul, digits=3), ylab=expression(min^-1), cex.lab=1.5, cex.main=1.5, cex.axis=1.5)
        image(legend, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), zlim=c(0,MAPul), xaxt="n", yaxt="n", xlab="0", main=format(MAPul, digits=3), ylab=expression(min^-1), cex.lab=1.5, cex.main=1.5, cex.axis=1.5, add=TRUE)

        for(i in 1:length(files)){
            if(files[i]!=0){
                image(MAP.list[[i]], col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), xlim=c(xmin/nx.vec[i],xmax/nx.vec[i]), ylim=c(ymin/ny.vec[i],ymax/ny.vec[i]), zlim=c(0,MAPul), main=paste("Median ", plot.param,"=", format(median(MAP.list[[i]][MAP.list[[i]]>0], na.rm=TRUE), digit=2), "     ", IDvisit.vec[i], "     ROI=", slice.vec[i], sep=""), cex.axis=1.5)
                image(MAP.list[[i]], col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), xlim=c(xmin/nx.vec[i],xmax/nx.vec[i]), ylim=c(ymin/ny.vec[i],ymax/ny.vec[i]), zlim=c(0,MAPul), main=paste("Median ", plot.param,"=", format(median(MAP.list[[i]][MAP.list[[i]]>0], na.rm=TRUE), digit=2), "     ", IDvisit.vec[i], "     ROI=", slice.vec[i], sep=""), cex.axis=1.5, add=TRUE)
            }
        }

        image(array(1:2, dim=c(1,2)), col=palette(colorRampPalette(c("red","blue"))(2)), zlim=c(0,2), xaxt="n", yaxt="n", xlab=modeltype1, main=modeltype2, cex.lab=1.5, cex.main=1.5, cex.axis=1.5)
        image(array(1:2, dim=c(1,2)), col=palette(colorRampPalette(c("red","blue"))(2)), zlim=c(0,2), xaxt="n", yaxt="n", xlab=modeltype1, main=modeltype2, cex.lab=1.5, cex.main=1.5, cex.axis=1.5, add=TRUE)

        for(i in 1:length(files)){
            if(files[i]!=0){
                image(AIC.list[[i]], xlim=c(xmin/nx.vec[i],xmax/nx.vec[i]), ylim=c(ymin/ny.vec[i],ymax/ny.vec[i]), col=palette(colorRampPalette(c("black", "red","blue"))(3)), main="AIC test", cex.axis=1.5, cex.lab=1.5)
                image(AIC.list[[i]], xlim=c(xmin/nx.vec[i],xmax/nx.vec[i]), ylim=c(ymin/ny.vec[i],ymax/ny.vec[i]), col=palette(colorRampPalette(c("black", "red","blue"))(3)), main="AIC test", cex.axis=1.5, cex.lab=1.5, add=TRUE)
            }
        }

        frame()

        for(i in 1:length(files)){
            if(files[i]!=0){
                if(modelstructure.vec[i]==1){
                    plot(maptimes.list[[i]], ccmedian.list[[i]], xlab="min", ylab="contrast agent", main="median contrast agent conc; xTofts=red", cex=3, cex.axis=1.5, cex.lab=1.5)
                    lines(maptimes.list[[i]], roimedianfittedxTofts.list[[i]], col="red", lwd=5)
                }
                if(modelstructure.vec[i]==2){
                    plot(maptimes.list[[i]], ccmedian.list[[i]], xlab="min", ylab="contrast agent", main="median contrast agent conc; Tofts=blue, xTofts=red", cex=3, cex.axis=1.5, cex.lab=1.5)
                    lines(maptimes.list[[i]], roimedianfittedxTofts.list[[i]], col="red", lwd=5)
                    lines(maptimes.list[[i]], roimedianfittedTofts.list[[i]], col="blue", lwd=2)
                }
                text(x=0.6*max(maptimes.list[[i]]), y=0.30*max(ccmedian.list[[i]], na.rm=TRUE), labels="EXTENDED TOFTS PARAMS", cex=1.5)
                text(x=0.6*max(maptimes.list[[i]]), y=0.23*max(ccmedian.list[[i]], na.rm=TRUE), labels=paste(expression(K^trans), " = ", format(paramsxTofts.list[[i]][1], digits=3), " 1/min (", cvroi.list[[i]][1], "%)", sep=""), cex=1.5)
                text(x=0.6*max(maptimes.list[[i]]), y=0.16*max(ccmedian.list[[i]], na.rm=TRUE), labels=paste(expression(k_ep), " = ", format(paramsxTofts.list[[i]][2], digits=3), " 1/min (", cvroi.list[[i]][2], "%)", sep=""), cex=1.5)
                text(x=0.6*max(maptimes.list[[i]]), y=0.09*max(ccmedian.list[[i]], na.rm=TRUE), labels=paste(expression(v_b), " = ", format(paramsxTofts.list[[i]][3], digits=3), " dimensionless (", cvroi.list[[i]][3], "%)", sep=""), cex=1.5)
            }
        }

        dev.off()
    }

    cat("...done.", sep="", "\n")
    cat("--------", "\n")
}
