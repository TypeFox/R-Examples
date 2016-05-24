#' A function for identifying spectral outliers based on Multidimensional Scaling (MDS) approach
#'
#' @author Andrew Sila \email{asila@cgiar.org} and Tomislav Hengl \email{tom.hengl@wur.nl}
spectra.outliers<-function (file.name, spectra = "non.standard", plots = TRUE, smethod = "euclidean", k = 5) 
{
    message("Reading spectra...")
    spec <- read.csv(file.name)
    std.slots<-c("metadata","data","sp")
    std.slots.a<-which(slotNames(spec)%in%std.slots)
    
    if(length(std.slots.a)<1){
    	names(spec)[1] = "SAMPLEID"
	samples <- cbind(spec["SAMPLEID"], MID="AfSIS-IR", DateTime=Sys.time())
	spec<- Spectra(samples, spec, replace.prefix="X")}

    ss <- spec@ab[, -1]
    wl <- as.numeric(substr(colnames(ss), 2, 19))
    colnames(ss) <- wl
    message("Calculating first derivative using Savitzky_Golay Algorithm...")
    deriv.spec <- as.data.frame(trans(ss)$trans)
    message("First derivative computation completed successifuly!")
    x <- data.frame(lapply(deriv.spec, function(x) {
        as.numeric(paste(x))
    }))
    colnames(x) <- paste0(substr(colnames(spec@ab), 1, 1)[200], colnames(x))
   	der.ssn <- cbind(as.data.frame(cbind(as.vector(spec@ab[, 1]), x)))
    colnames(der.ssn) <- c(colnames(spec@ab)[1], paste0(substr(colnames(spec@ab), 
    1, 1)[200], colnames(deriv.spec)))
    message("Computing spectral similarity distance ")
    dist.spec <- dist(der.ssn[, -1], smethod)
    dist <- cmdscale(dist.spec, k = 5)
    mn <- mean(dist)
    sd <- sd(dist)
    if (spectra == "standard") {
        cols <- ifelse(dist > mn + sd | dist < mn - sd, "red", 
            "blue")[, 1]
        UL <- mn + sd
        LL <- mn - sd
    }
    if (spectra == "non.standard") {
        cols <- ifelse(dist > mn + (3 * sd) | dist < mn - (3 * 
            sd), "red", "grey")[, 1]
        UL <- mn + (3 * sd)
        LL <- mn - (3 * sd)
    }
    distt <- as.data.frame(dist)
    x <- data.frame(lapply(distt, function(x) {
        sd(x)
    }))
    out.select <- menu(c("Yes", "No - I want to select my own"), 
        graphics = TRUE, title = "Use default outliers")
    close.screen(all.screens = TRUE)
    screen(1)
    split.screen(c(1, 2))
    require(lava)
    if (out.select == 1) {
        plot(dist[, 1], pch = 19, xlab = "Spectral points", ylab = "Similarity distance", 
            ylim = c(min(dist) - 3 * sd, max(dist) + 3 * sd), 
            col = Col(cols, 0.5))
        abline(h = UL, lty = 3, col = 3)
        abline(h = LL, lty = 2, col = 3)
        legend("topleft",c("Upper Control Limit", "Mean similarity distance", " Lower Control Limit"),lty=c(3,1,2),col=4,bty="n")
        mtext("MDS computed outliers", cex = 1, col = "blue4")
        all.outliers <- which(cols == "red")
        split.screen(c(2, 1), screen = 2)
        spec.outliers <- spec@ab[all.outliers, ]
        ## Create Spectra objects for plotting
		samples <- cbind(spec.outliers ["SAMPLEID"], MID="AfSIS-IR", DateTime=Sys.time())
		xo <- Spectra(samples, spec.outliers , replace.prefix="X")
		plot(xo)
        mtext("Spectral outliers", side = 1, line = 3, cex = 1.3, col = "grey")
        screen(4)
        spec.normal <- spec@ab[-all.outliers, ]
        sel <- sample(1:nrow(spec.normal), round(1 * nrow(spec.normal), 
            0), replace = FALSE, prob = NULL)
         #Create Spectra objects for plotting
		samples <- cbind(spec.normal ["SAMPLEID"], MID="AfSIS-IR", DateTime=Sys.time())
		xn <- Spectra(samples, spec.normal , replace.prefix="X")
		plot(xn)
        mtext("Normal spectra", side = 1, line = 3, cex = 1.3, col = "grey")
    }
    if (out.select == 2) {
        close.screen(all.screens = TRUE)
        split.screen(c(1, 2))
        screen(1)
        plot(dist[, 1], pch = 19, xlab = "Spectral points", ylab = "Similarity distance", ylim = c(min(dist) - 3 * sd, max(dist) + 3 * sd), col = Col(cols, 0.5))
        abline(h = UL, lty = 3, col = 3)
        abline(h = LL, lty = 2, col = 3)
        outliers <- identify(dist[, 1], labels = spec@ab[, 1], cex = 0.6)
        mtext("Click to pick spectral outliers", cex = 1, col = "blue")
        ssno <- as.vector(spec@ab[outliers, 1])
        split.screen(c(2, 1), screen = 2)
        spec.outliers <- spec@ab[outliers, ]
        #Create Spectra objects for plotting
		samples <- cbind(spec.outliers ["SAMPLEID"], MID="AfSIS-IR", DateTime=Sys.time())
		xo <- Spectra(samples, spec.outliers , replace.prefix="X")
		plot(xo)
        mtext("Spectral outliers", side = 1, line = 3, cex = 1.3, col = "red")
        screen(4)
        spec.normal <- spec@ab[-outliers, ]
        sel <- sample(1:nrow(spec.normal), round(0.1 * nrow(spec.normal), 0), replace = FALSE, prob = NULL)
		#Create Spectra objects for plotting
		samples <- cbind(spec.normal ["SAMPLEID"], MID="AfSIS-IR", DateTime=Sys.time())
		xn <- Spectra(samples, spec.normal , replace.prefix="X")
		plot(xn)
        mtext("Normal spectra", side = 1, line = 3, cex = 1.3, col = "blue")
        more <- menu(c("Yes", "No"), graphics = TRUE, title = "Pick more outliers?")
        if (more == 1) {
            close.screen(all.screens = TRUE)
            screen(1)
            split.screen(c(1, 2))
            plot(dist[, 1], pch = 19, xlab = "Spectral points", ylab = "Similarity distance", ylim = c(min(dist) - 
             3 * sd, max(dist) + 3 * sd), col = Col(cols, 
                  0.6))
            abline(h = UL, lty = 3, col = 3)
            abline(h = LL, lty = 2, col = 3)
            legend("topleft",c("Upper Control Limit", "Mean similarity distance", " Lower Control Limit"),lty=c(3,1,2),col=4,bty="n")

            mtext("Click to pick spectral outliers", cex = 0.8, col = "blue")
            outliers <- identify(dist[, 1], labels = spec@ab[, 1], 
                cex = 0.6)
            split.screen(c(2, 1), screen = 2)
             spec.outliers <- spec@ab[outliers, ]
        	#Create Spectra objects for plotting
			samples <- cbind(spec.outliers ["SAMPLEID"], MID="AfSIS-IR", DateTime=Sys.time())
			xo <- Spectra(samples, spec.outliers , replace.prefix="X")
			plot(xo)
            mtext("Spectral outliers", side = 1, line = 3, cex = 1.3, col = "red")
            screen(4)
            spec.normal <- spec@ab[-outliers, ]
        	sel <- sample(1:nrow(spec.normal), round(0.1 * nrow(spec.normal), 0), replace = FALSE, prob = NULL)
			#Create Spectra objects for plotting
			samples <- cbind(spec.normal ["SAMPLEID"], MID="AfSIS-IR", DateTime=Sys.time())
			xn <- Spectra(samples, spec.normal , replace.prefix="X")
			plot(xn)
      
            mtext("Normal spectra", side = 1, line = 3, cex = 1.3, col = "blue")
        }
        if (more == 2) {
            outliers = outliers
        }
        all.outliers <- outliers
    }
   #Check for trend	 along the other dimensions
   	require(date)
	if(spectra=="standard"){
	dat<-spec@ab[,1]
	date<-as.Date(dat) 
		par(mfrow=c(1,1))
	#Define point colors using the ranges
	plot(date,dist[,1],ylim=c(min(dist[,1]),max(dist[,1])),ylab="",xlab="",col=Col(cols[1],0.3),pch=16)
	abline(h=UL,lty=3,col=4,lwd=1.3)
	abline(h=mn,lty=1,col=4,lwd=0.8)
	abline(h=LL,lty=2,col=4,lwd=1.3)
	mtext("Similarity distance trend over time",cex=1)
	mtext("Measurement dates",side=1,cex=1,line=2)
	mtext("Similarity distance",side=2,cex=1,line=2)
	legend("topleft",c("Upper Control Limit", "Mean similarity distance", " Lower Control Limit"),lty=c(3,1,2),col=4,bty="n")
	}
	
    if (spectra == "non.standard") {
        par(mfrow = c(3, 2))
        for (j in 1:4) {
            plot(dist[, j], ylim = c(min(dist[, j]), max(dist[, 
                j])), ylab = "Similarity distance", xlab = "Samples", 
                col = Col(cols, 0.4), pch = 16)
            mtext(paste("Spectral data trend along", j, " dimension", 
                sep = ""), cex = 0.7, col = "blue")
            abline(h = UL, lty = 3, col = 2)
            abline(h = LL, lty = 2, col = 2)
             legend("topleft",c("Upper Control Limit", "Mean similarity distance", " Lower Control Limit"),lty=c(3,1,2),col=4,bty="n")

        }
    }
    if (out.select == 1) {
        out.ids <- spec@ab[all.outliers, 1]
    }
    if (out.select == 2) {
        out.ids <- spec@ab[all.outliers, 1]
    }
    out <- as.data.frame(cbind(c(1:length(out.ids)), as.vector(out.ids)))
    colnames(out) <- c("No", "SSN")
    spec.out <- as.data.frame(spec@ab[all.outliers, ])
    spec.normal <- as.data.frame(spec@ab[-all.outliers, ])
    output <- list(dist = dist, all.outliers = all.outliers, 
        spec.out = spec.out, spec.normal = spec.normal)
    class(output) <- "spectra.outliers"
    return(output)
}
setMethod("spectra.outliers",signature("Spectra"),spectra.outliers)#end script
