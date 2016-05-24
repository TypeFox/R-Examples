# For making binary templates
# Modified: 2015 Sept 17

makeBinTemplate <-
function(
   clip,               # File path to wav or mp3 file, or a list of exactly two file paths
   t.lim=NA,           # Length two vector of time limits of spectrogram plot or template, or a list of exactly two such vectors
   frq.lim=c(0, 12),    # Frequency limits of spectrogram plot or template 
   select="auto",      # How should points be selected? Options are "cell" or "click" (the same), "rectangle", "auto"
   binary=TRUE,        # True for a binary plot with high and low amplitude cells
   buffer=0,           # Buffer around "on" points for "rectangle" and "auto"
   dens=1,             # Density of points included with "rectangle" and "auto" (fraction of 1)
   score.cutoff=12,    # Score cutoff for template
   name="A",           # Name of template
   comment="", 
   amp.cutoff="i",     # "i" for interactive, otherwise, a length-one numeric vector
   shift="i",          # Forward shift for the second clip, to align with the first, in time bins, or "i" for interactive
   high.pass=-Inf,     # Sets all amplitudes below this frequency to minimum
   spec.col=gray.3(),  # Color palette for spectrogram when binary=FALSE
   bin.col=c("white", "black"), # Colors for binary plot
   quat.col=c("white", "gray40", "gray75", "black"), # Colors for quaternary plot (two clips)
   sel.col=c("orange", "blue"), # Colors for selected points (on, off)
   legend.bg.col="#2E2E2E94", # Legend background color
   legend.text.col="black", # Legend text color
   wl=512,             # Window length for spectro
   ovlp=0,             # % overlap between windows for spectro
   wn="hanning",       # Window type for spectro
   write.wav=FALSE,    # Set to TRUE to allow writing clip wave objects to file
   ...                 # Additional arguments to spectro
){

   # Check some arguments
   if(!binary & sum(grepl(select, c("rectangle", "automatic"))) != 0) { 
      warning("binary adjusted to TRUE for select=\"rectangle\" or select=\"auto\"", immediate.=TRUE)
      binary <- TRUE
   }
   if(!binary & buffer>0) 
      warning("buffer argument ignored unless binary=TRUE", immediate.=TRUE)
   if(select%in%c("cell", "click") & dens<1) 
      warning("dens argument ignored for select=\"click\"", immediate.=TRUE)
   if(dens<0.0001 | dens>1) {
      warning("dens adjusted to 1.0", immediate.=TRUE)
      dens <- 1
   }
   if(length(clip) == 2 & !binary)
      warning("binary adjusted to TRUE for two clips", immediate.=TRUE)
   if(class(amp.cutoff) != "numeric" & amp.cutoff != "i") {
      warning("amp.cutoff value not recognized, so set to \"i\" for interactive", immediate.=TRUE)
      amp.cutoff="i"
   }
   if(class(shift) != "numeric" & shift != "i") {
      warning("shift value not recognized, so set to \"i\" for interactive", immediate.=TRUE)
      shift="i"
   }
   if(!select%in%c("cell", "click", "auto", "rectangle", "rect")) stop("select argument, ", select, ", not recognized")

   # Creates a wav file for any clip elements that are not already files
   clip <- getClip(clip, name=deparse(substitute(clip)), write.wav=write.wav)

##### Single clip ##### 
   if(length(clip) == 1) { 
      clip.path <- clip
      clip <- readClip(clip)

      # Trim clip
      if(is.na(t.lim[1])) {
	t.lim <- c(0, Inf) 
      } else {
	  clip <- cutWave(clip, from=t.lim[1], to=t.lim[2])
      }
      samp.rate <- clip@samp.rate
      first.t.bin <- t.lim[1]
      # Fourier transform
      fspec <- spectro(wave=clip, wl=wl, ovlp=ovlp, wn=wn, ...)
      
      # Sort out time and frequency bins
      t.bins <- fspec$time
      n.t.bins <- length(t.bins)
      which.t.bins <- 1:n.t.bins
      which.frq.bins <- which(fspec$freq >= frq.lim[1] & fspec$freq <= frq.lim[2])
      frq.bins <- fspec$freq[which.frq.bins]
      n.frq.bins <- length(frq.bins)

      # Filter amplitudes
      amp <- fspec$amp[which.frq.bins, ]
      amp[frq.bins<high.pass, ] <- min(c(amp))

      # Create empty amplitude matrices for plotting
      on.mat <- off.mat <- matrix(0, nrow=n.frq.bins, ncol=n.t.bins)

      # Bin steps
      t.step <- t.bins[2]-t.bins[1]
      frq.step <- frq.bins[2]-frq.bins[1]

      # Determine amp cutoff for binary plot
      if(binary && amp.cutoff == "i") {
         select.cutoff <- TRUE
         amp.cutoff <- round(quantile(amp, 0.7))
      } else 
         select.cutoff <- FALSE

      # Make plot
      oldpar <- par(mar=c(5, 4,4, 4))

      # Color ramp plot
      if(!binary) {
        image(x=which.t.bins, y=which.frq.bins, t(amp), col=spec.col, xlab="Time (s)", ylab="Frequency (kHz)", las=1, useRaster=TRUE, axes=FALSE, las=1)
        t.bin.ticks <- pretty(t.bins, n=6)
        axis(1, at=t.bin.ticks/t.step, labels=t.bin.ticks+t.lim[1])
        frq.bin.ticks <- pretty(frq.bins, n=6)
        axis(2, at=frq.bin.ticks/frq.step, labels=frq.bin.ticks, las=1)
        axis(3)
        axis(4, las=1)
        box()
      } else {
      # Binary plot, redrawn if amplitude is interactively selected
        y <- 0
        if(select.cutoff) {
           cat("\nInteractive amplitude cutoff selection.") 
           cat("\nEnter l, ll, ll, etc. for lower cutoff, \nh, hh, hhh, etc. for higher cutoff, \nor hit Enter to continue\n") 
        }
        while(y != "") {
          bin.amp <- matrix(0, nrow=n.frq.bins, ncol=n.t.bins)
          bin.amp[amp>amp.cutoff] <- 1
          image(x=which.t.bins, y=which.frq.bins, t(bin.amp), col=bin.col, xlab="Time (s)", ylab="Frequency (kHz)", las=1, useRaster=TRUE, axes=FALSE, las=1)
          legend("topleft", c(paste("Amplitude cutoff: ", amp.cutoff)), bg=legend.bg.col)
          t.bin.ticks <- pretty(t.bins, n=6)
          axis(1, at=t.bin.ticks/t.step, labels=t.bin.ticks+t.lim[1])
          frq.bin.ticks <- pretty(frq.bins, n=6)
          axis(2, at=frq.bin.ticks/frq.step, labels=frq.bin.ticks, las=1)
          axis(3)
          axis(4, las=1)
          box()
          # Amplitude cutoff selection
          if(select.cutoff) {
             cat("\nCurrent cutoff: ", amp.cutoff, "\n", sep="")
             #cat("\nCurrent cutoff: ", amp.cutoff, ":", sep="")
             #y <- tolower(scan(n=1, what="character", quiet=TRUE))
             y <- tolower(readLines(n=1))
             #if(length(y) == 1) 
             if(y != "") 
                amp.cutoff <- switch(y, 0,"l"=-1, "ll"=-3, "lll"=-6, "llll"=-10, "lllll"=-20, "llllll"=-30, "h"=1, "hh"=3, "hhh"=6, "hhhh"=10, "hhhhh"=20, "hhhhhh"=30) + amp.cutoff
          } else 
             y <- ""
        }
      } 

      # Point-by-point selection
      if(select%in%c("cell", "click")) {
	if(grepl('[Xx]11', .Device)) {
          cat("\nSelect \"on\" points with left mouse click. To continue, right click.\n")
	} else {
          cat("\nSelect \"on\" points with left mouse click. To continue, press \'ESC\'.\n")
	}
      }
      if(select%in%c("cell", "click")) {
         # Plot over legend
         if(!binary)
            image(x=which.t.bins, y=which.frq.bins, t(amp), col=spec.col, add=TRUE)
         else
            image(x=which.t.bins, y=which.frq.bins, t(bin.amp), col=bin.col, zlim=c(0, 1), add=TRUE)
         # Select "on" points
         i <- 0
         pts <- NULL
         pos <- NULL
         while(!is.null(pos)|i == 0) {
            i <- i + 1
            pos <- locator(n=1)
            if(!is.null(pos)) pos <- lapply(pos, round)
            if(!is.null(pos)) pos$y <- pos$y - min(which.frq.bins) + 1
            pts <- rbind(pts, as.numeric(pos))
            on.mat[pts[, 2:1, drop=FALSE]] <- 1
            # Add points to plot
            image(x=which.t.bins, y=which.frq.bins, t(on.mat), col=c("transparent", sel.col[1]), zlim=c(0, 1), add=TRUE)
            if(!is.null(pos)) 
               cat("\n", nrow(pts), " selected")
         } 
         pt.on <- pts
         pt.on <- unique(pt.on)
         colnames(pt.on) <- c("t", "frq")
         pt.on[, "frq"] <- pt.on[, "frq"] + min(which.frq.bins) - 1
 
         # Select "off" points
         i <- 0
	if(grepl('[Xx]11', .Device)) {
          cat("\nSelect \"off\" points with left mouse click. When done, right click.\n")
	} else {
          cat("\nSelect \"off\" points with left mouse click. When done, press \'ESC\'.\n")
	}
         pts <- NULL
         pos <- NULL
         while(!is.null(pos)|i == 0) {
            i <- i + 1
            pos <- locator(n=1)
            if(!is.null(pos)) pos <- lapply(pos, round)
            if(!is.null(pos)) pos$y <- pos$y - min(which.frq.bins) + 1
            pts <- rbind(pts, as.numeric(pos))
            off.mat[pts[, 2:1, drop=FALSE]] <- 1
            # Add points to plot
            image(x=which.t.bins, y=which.frq.bins, t(off.mat), col=c("transparent", sel.col[2]), zlim=c(0, 1), add=TRUE)
            if(!is.null(pos)) 
               cat("\n", nrow(pts), " selected")
         } 
         pt.off <- pts
         pt.off <- unique(pt.off)
         colnames(pt.off) <- c("t", "frq")
         pt.off[, "frq"] <- pt.off[, "frq"] + min(which.frq.bins) - 1

      } else if(select%in%c("rect", "rectangle")) {
         # Rectangular selection
         cat("\nRectangular selection\n")
         # Plot over legend
         image(x=which.t.bins, y=which.frq.bins, t(bin.amp), col=bin.col, zlim=c(0, 1), add=TRUE)
         i <- 0
         pos1 <- pt.on <- NULL
         while(!is.null(pos1)|i == 0) {
            # On cells first
	    if(grepl('[Xx]11', .Device)) {
              cat("\nSelect upper left corner of \"on\" rectangle with a left click.\nRight click to continue.\n")
	    } else {
              cat("\nSelect upper left corner of \"on\" rectangle with a left click.\nPress \'ESC\' to continue.\n")
	    }
            i <- i + 1
            pos1 <- locator(n=1)
            points(pos1$x, pos1$y, pch=22, cex=0.5, col="red", bg="red")
            if(!is.null(pos1)) {
               cat("\nSelect lower right corner of \"on\" rectangle with a left click.\n")
               pos2 <- locator(n=1)
               points(pos2$x, pos2$y, pch=22, cex=0.5, col="red", bg="red")

               # First find positions within the matrix that are within the rectangle
               frq.in.rect <- which.frq.bins<pos1$y & which.frq.bins>pos2$y
               x.in.rect <- which.t.bins>pos1$x & which.t.bins<pos2$x

               # Set cells that meet criteria to 1 in on.mat
               temp.mat <- on.mat
               temp.mat[frq.in.rect, x.in.rect] <- bin.amp[frq.in.rect, x.in.rect]
               temp.mat[temp.mat == 1] <- sample(c(1, 0), sum(temp.mat == 1), TRUE, c(dens, 1-dens))
               on.mat[temp.mat == 1] <- 1

               # Then find locations of all high amplitude cells within rectangle
               pts <- which(on.mat == 1, arr.ind=TRUE)
               pts <- pts[, 2:1]
               colnames(pts) <- c("t", "frq")
               pts[, "frq"] <- pts[, "frq"] + min(which.frq.bins) - 1
               pt.on <- rbind(pt.on, pts)

               # Plot over rectangle corners
               image(x=which.t.bins, y=which.frq.bins, t(bin.amp), col=bin.col, zlim=c(0, 1), add=TRUE)
               # Add points to plot
               image(x=which.t.bins, y=which.frq.bins, t(on.mat), col=c("transparent", sel.col[1]), zlim=c(0, 1), add=TRUE)
            }
            pt.on <- unique(pt.on)
         }
         # Off cells
         i <- 0
         pos1 <- pt.off <- NULL
         while(!is.null(pos1)|i == 0) {
            #message()
	    if(grepl('[Xx]11', .Device)) {
              cat("\nSelect upper left corner of \"off\" rectangle with a left click.\nRight click to continue.\n")
	    } else {
              cat("\nSelect upper left corner of \"off\" rectangle with a left click.\nPress \'ESC\' to continue.\n")
	    }
            i <- i + 1
            pos1 <- locator(n=1)
            points(pos1$x, pos1$y, pch=22, cex=0.5, col="red", bg="red")
            if(!is.null(pos1)) {
               #message()
               cat("\nSelect lower right corner of \"off\" rectangle with a left click.\n")
               pos2 <- locator(n=1)
               points(pos2$x, pos2$y, pch=22, cex=0.5, col="red", bg="red")

               # First find positions within the matrix that are within the rectangle
               frq.in.rect <- which.frq.bins<pos1$y & which.frq.bins>pos2$y
               x.in.rect <- which.t.bins>pos1$x & which.t.bins<pos2$x
               
               # Find cells with high binary cells within buffer distance--these will have 0 in buff.amp
               # This loop is a bit slow
               buff.amp <- matrix(TRUE, nrow=n.frq.bins, ncol=n.t.bins)
               if(buffer>0) {
                  for(i in 1:n.frq.bins) {
                     for(j in 1:n.t.bins) {
                        # Next line looks for cells with 1 (high amplitude) anywhere within a matrix 1 + 2*buffer square, including the center cell [i, j]
                        # The min expression prevents subscripts < 0
                        buff.amp[i, j] <- (sum(bin.amp[i+ -min(buffer, i-1):min(buffer, n.frq.bins-i), j+ -min(buffer, j-1):min(buffer, n.t.bins-j)] == 1) == 0)
                     }
                  }
               }

               # Set cells that meet all criteria to 1 in off.mat
               temp.mat <- off.mat + 1
               temp.mat[frq.in.rect, x.in.rect] <- bin.amp[frq.in.rect, x.in.rect]
               temp.mat[temp.mat == 0] <- sample(c(1, 0), sum(temp.mat == 0), TRUE, c(1-dens, dens)) # Note that probability order is reversed compared to "on" points
               off.mat[temp.mat == 0 & buff.amp] <- 1 

               # Then find locations of cells
               pts <- which(off.mat == 1, arr.ind=TRUE)
               pts <- pts[, 2:1]
               colnames(pts) <- c("t", "frq")
               pts[, "frq"] <- pts[, "frq"] + min(which.frq.bins) - 1
               pt.off <- rbind(pt.off, pts)

               # Plot over rectangle corners
               image(x=which.t.bins, y=which.frq.bins, t(bin.amp), col=bin.col, zlim=c(0, 1), add=TRUE)
               # Add points to plot
               image(x=which.t.bins, y=which.frq.bins, t(on.mat), col=c("transparent", sel.col[1]), zlim=c(0, 1), add=TRUE)
               image(x=which.t.bins, y=which.frq.bins, t(off.mat), col=c("transparent", sel.col[2]), zlim=c(0, 1), add=TRUE)
            }
            pt.off <- unique(pt.off)
         } 
      } else if(select %in% c("auto", "automatic")) {
         # Automatic point selection
         cat("\nAutomatic point selection.\n")

         # Plot over legend
         image(x=which.t.bins, y=which.frq.bins, t(bin.amp), col=bin.col, zlim=c(0, 1), add=TRUE)

         # On cells first
         # Set cells that meet criteria to 1 in on.mat, incorporating random cell selection from among these cells based on dens
         on.mat[bin.amp == 1] <- sample(c(1, 0), sum(bin.amp == 1), TRUE, c(dens, 1-dens))

         # Then find locations of all high binary cells 
         pts <- which(on.mat == 1, arr.ind=TRUE)
         pts <- pts[, 2:1]
         colnames(pts) <- c("t", "frq")
         pts[, "frq"] <- pts[, "frq"] + min(which.frq.bins) - 1
         pt.on <- pts

         # Off cells next
         # Find cells with high binary cells within buffer distance--these will have FALSE in buff.amp
         buff.amp <- matrix(TRUE, nrow=n.frq.bins, ncol=n.t.bins)
         if(buffer>0) {
            for(i in 1:n.frq.bins) {
               for(j in 1:n.t.bins) {
                  # Next line looks for cells with 1 anywhere within a matrix 1 + 2*buffer square, including the center cell [i, j]
                  # The min expression prevents subscripts < 0
                  buff.amp[i, j] <- (sum(bin.amp[i+ -min(buffer, i-1):min(buffer, n.frq.bins-i), j+ -min(buffer, j-1):min(buffer, n.t.bins-j)] == 1) == 0)
               }
            }
         }

         # Set cells that meet criteria to 1 in off.mat, incorporating random point selection based on dens
	 # Remember TRUE cells in buff.amp are outside the buffer (and so OK for "off" cells)
         off.mat[bin.amp == 0 & buff.amp] <- sample(c(1, 0), sum(bin.amp == 0 & buff.amp), TRUE, c(dens, 1-dens)) # Here 1 in off.mat means cell is "off" cell

         # Then find locations of all low binary cells 
         pts <- which(off.mat == 1, arr.ind=TRUE)
         pts <- pts[, 2:1]
         colnames(pts) <- c("t", "frq")
         pts[, "frq"] <- pts[, "frq"] + min(which.frq.bins) - 1
         pt.off <- pts

         # Add points to plot
         image(x=which.t.bins, y=which.frq.bins, t(on.mat), col=c("transparent", sel.col[1]), zlim=c(0, 1), add=TRUE)
         image(x=which.t.bins, y=which.frq.bins, t(off.mat), col=c("transparent", sel.col[2]), zlim=c(0, 1), add=TRUE)
      }
   } else if(length(clip) == 2) {
##### Two clips ##### 
      if(is.na(t.lim[1])) t.lim <- list(NA, NA)
      if(class(t.lim) != "list") stop("Supplied two clips, so t.lim must be a list of length 2.")
#      if(class(t.lim) != "list" || length(t.lim) != 2) stop("Supplied two clips, so t.lim must be a list of length 2.")
      clip.path <- clip[[2]]
      first.t.bin <- t.lim[[2]][2]
      if(is.na(first.t.bin)) first.t.bin <- 0
      samp.rate <- fspec <- t.bins <- n.t.bins <- t.step <- which.t.bins <- amp <- list()
      which.frq.bins <- frq.bins <- frq.step <- n.frq.bins <- list()
      for(i in 1:2) {
         clip[[i]] <- readClip(clip[[i]])

         # Trim clip[[i]]
         if(is.na(t.lim[[i]][1])) 
            t.lim[[i]] <- c(0, Inf) else  
            clip[[i]] <- cutWave(clip[[i]], from=t.lim[[i]][1], to=t.lim[[i]][2])
         samp.rate[[i]] <- clip[[i]]@samp.rate

         # Fourier transform
         fspec[[i]] <- spectro(wave=clip[[i]], wl=wl, ovlp=ovlp, wn=wn, ...)
         
         # Filter amplitudes 
         t.bins[[i]] <- fspec[[i]]$time
         t.step[[i]] <- t.bins[[i]][2] - t.bins[[i]][1]
         n.t.bins[[i]] <- length(t.bins[[i]])
         which.t.bins[[i]] <- 1:n.t.bins[[i]]
         which.frq.bins[[i]] <- which(fspec[[i]]$freq >= frq.lim[1] & fspec[[i]]$freq <= frq.lim[2])
         frq.bins[[i]] <- fspec[[i]]$freq[which.frq.bins[[i]]]
         frq.step[[i]] <- frq.bins[[i]][2] - frq.bins[[i]][1]
         amp[[i]] <- fspec[[i]]$amp[which.frq.bins[[i]], ]
         amp[[i]][frq.bins[[i]]<high.pass, ] <- min(c(amp[[i]]))
      }

      # Check sampling rate
      if(clip[[1]]@samp.rate != clip[[2]]@samp.rate) stop("Sampling rates do not match.")
      samp.rate <- clip[[1]]@samp.rate

      # Frequency
      frq.step <- frq.step[[1]]
      frq.bins <- frq.bins[[1]]
      which.frq.bins <- which.frq.bins[[1]]
      n.frq.bins <- length(frq.bins)

      # Adjust size of smaller amplitude matrix to match the sizes
      if(n.t.bins[[1]] > n.t.bins[[2]]) {
         t.bins <- t.bins[[1]]
         which.t.bins <- which.t.bins[[1]]
         zero.mat <- matrix(min(c(amp[[2]])), nrow=n.frq.bins, ncol=n.t.bins[[1]] - n.t.bins[[2]])
         amp[[2]] <- cbind(amp[[2]], zero.mat)
      } else if(n.t.bins[[2]] > n.t.bins[[1]]) {
         t.bins <- t.bins[[2]]
         which.t.bins <- which.t.bins[[2]]
         zero.mat <- matrix(min(c(amp[[1]])), nrow=n.frq.bins, ncol=n.t.bins[[2]] - n.t.bins[[1]])
         amp[[1]] <- cbind(amp[[1]], zero.mat)
      } else {
         t.bins <- t.bins[[1]]
         which.t.bins <- which.t.bins[[1]]
      }

      n.t.bins <- length(t.bins)
      t.lim <- c(0, max(diff(t.lim[[1]]), diff(t.lim[[1]])))
      t.step <- t.step[[1]]

      # Amp cutoff for binary plot
      if(amp.cutoff == "i") {
         select.cutoff <- TRUE
         amp.cutoff <- round(quantile(c(amp[[1]], amp[[2]]), 0.7))
      } else select.cutoff <- FALSE

      # Make plot
      oldpar <- par(mar=c(5, 4,4, 4))

      # First plot for selecting cutoff
      x <- 0
      if(select.cutoff) {
         cat("\nInteractive amplitude cutoff selection.") 
         cat("\nEnter l, ll, ll, etc. for lower cutoff, \nh, hh, hhh, etc. for higher cutoff, \nor hit Enter to continue\n") 
      }
      while(x != "") {
        # Create matrices of on/off data from amplitude data
        bin.amp <- list()
        for(i in 1:2) {
           bin.amp[[i]] <- matrix(0, nrow=nrow(amp[[i]]), ncol=ncol(amp[[i]]))
           bin.amp[[i]][amp[[i]]>amp.cutoff] <- i
        }
        # mat3 is main plotted matrix. Key to mat3: 1 = clip 1 above cutoff, 2 = clip 2 above cutoff, 3 = both clips above cutoff, 0 = no clips above cutoff
        mat3 <- bin.amp[[1]] + bin.amp[[2]]
        n.on <- sum(mat3 == 3)
        n.off <- sum(mat3 == 0)
        # First plot
        image(x=which.t.bins, y=which.frq.bins, t(mat3), col=quat.col, zlim=c(0, 3), xlab="Time (s)", ylab="Frequency (kHz)", useRaster=TRUE, axes=FALSE)
        legend("topleft", c(paste("Shift", ifelse(is.numeric(shift), shift, 0)), paste("No. overlapped cells", n.on), paste("No. empty cells", n.off), paste("Amplitude cutoff", amp.cutoff)), text.col=legend.text.col, bg=legend.bg.col)
        t.bin.ticks <- pretty(t.bins, n=6)
        axis(1, at=t.bin.ticks/t.step, labels=t.bin.ticks+t.lim[1])
        frq.bin.ticks <- pretty(frq.bins, n=6)
        axis(2, at=frq.bin.ticks/frq.step, labels=frq.bin.ticks, las=1)
        axis(3)
        axis(4, las=1)
        box()
        # Amplitude cutoff selection
        if(select.cutoff) {
          cat("\nCurrent cutoff: ", amp.cutoff, "\n", sep="")
          #cat("\nCurrent cutoff: ", amp.cutoff, "\n")
          #x <- tolower(scan(n=1, what="character", quiet=TRUE))
          x <- tolower(readLines(n=1))
          #if(length(x) == 1) amp.cutoff <- switch(x, 0,"l"=-1, "ll"=-3, "lll"=-6, "llll"=-10, "lllll"=-20, "llllll"=-30, "h"=1, "hh"=3, "hhh"=6, "hhhh"=10, "hhhhh"=20, "hhhhhh"=30) + amp.cutoff
          if(x != "") amp.cutoff <- switch(x, 0,"l"=-1, "ll"=-3, "lll"=-6, "llll"=-10, "lllll"=-20, "llllll"=-30, "h"=1, "hh"=3, "hhh"=6, "hhhh"=10, "hhhhh"=20, "hhhhhh"=30) + amp.cutoff
        } else x <- ""
      }

      # Second plot for aligning clips
      select.shift <- FALSE
      if(shift == "i") {
         select.shift <- TRUE
         shift <- 0
         cat("\nInteractive clip alignment.") 
         cat("\nEnter l, ll, ll, etc. for left shift, \nr, rr, rrr, etc. for right shift, \nor Enter to continue.\n")
      }
      x <- 0
      while(length(x)>0) {
         # Shift clips
         if(select.shift) {
            #x <- tolower(scan(n=1, what="character", quiet=TRUE))
            x <- tolower(readLines(n=1))
            #if(length(x) == 1) shift <- switch(x, 0,"l"=-1, "ll"=-3, "lll"=-6, "llll"=-10, "lllll"=-20, "r"=1, "rr"=3, "rrr"=6, "rrrr"=10, "rrrrr"=20) + shift else x <- NULL
            if(x != "") shift <- switch(x, 0,"l"=-1, "ll"=-3, "lll"=-6, "llll"=-10, "lllll"=-20, "r"=1, "rr"=3, "rrr"=6, "rrrr"=10, "rrrrr"=20) + shift else x <- NULL
         } else x <- NULL
         if(shift<0) {
            zero.mat <- matrix(0, nrow=n.frq.bins, ncol=-shift)
            mat1 <- cbind(zero.mat, bin.amp[[1]])
            mat2 <- cbind(bin.amp[[2]], zero.mat)
         } else if(shift>0) {
            zero.mat <- matrix(0, nrow=n.frq.bins, ncol=shift)
            mat1 <- cbind(bin.amp[[1]], zero.mat)
            mat2 <- cbind(zero.mat, bin.amp[[2]])
         } else if(shift == 0) {
            mat1 <- bin.amp[[1]]
            mat2 <- bin.amp[[2]]
         } else stop("Time shift is ", shift)
         
         # Adjust number of time bins for the shift
         n.t.bins <- length(t.bins)+abs(shift)
         which.t.bins <- 1:n.t.bins
         mat3 <- mat1 + mat2
         n.on <- sum(mat3 == 3)
         n.off <- sum(mat3 == 0)
         # Re-create plot
         image(x=which.t.bins, y=which.frq.bins, t(mat3), col=quat.col, zlim=c(0, 3), xlab="Time (s)", ylab="Frequency (kHz)", useRaster=TRUE, axes=FALSE)
         legend("topleft", c(paste("Shift", ifelse(is.numeric(shift), shift, 0)), paste("No. overlapped cells", n.on), paste("No. empty cells", n.off), paste("Amplitude cutoff", amp.cutoff)), text.col=legend.text.col, bg=legend.bg.col)
         t.bin.ticks <- pretty(t.bins, n=6)
         axis(1, at=t.bin.ticks/t.step, labels=t.bin.ticks+t.lim[1])
         frq.bin.ticks <- pretty(frq.bins, n=6)
         axis(2, at=frq.bin.ticks/frq.step, labels=frq.bin.ticks, las=1)
         axis(3)
         axis(4, las=1)
         box()
      }

      # Create empty amplitude matrices for plotting
      on.mat <- off.mat <- matrix(0, nrow=n.frq.bins, ncol=n.t.bins)

      # Point-by-point selection
      if(select%in%c("cell", "click")) {
	 if(grepl('[Xx]11', .Device)) {
           cat("\nSelect \"on\" points with left mouse click. To continue, right click.\n")
	 } else {
           cat("\nSelect \"on\" points with left mouse click. To continue, press \'ESC\'.\n")
	 }
      }
      if(select%in%c("cell", "click")) {
         # Plot over legend
         image(x=which.t.bins, y=which.frq.bins, t(mat3), col=quat.col, zlim=c(0, 3), add=TRUE)
         # Select "on" points
         i <- 0
         pts <- NULL
         pos <- NULL
         while(!is.null(pos)|i == 0) {
            i <- i + 1
            pos <- locator(n=1)
            if(!is.null(pos)) pos <- lapply(pos, round)
            if(!is.null(pos)) pos$y <- pos$y - min(which.frq.bins) + 1
            pts <- rbind(pts, as.numeric(pos))
            on.mat[pts[, 2:1, drop=FALSE]] <- 1
            # Add points to plot
            image(x=which.t.bins, y=which.frq.bins, t(on.mat), col=c("transparent", sel.col[1]), zlim=c(0, 1), add=TRUE)
            if(!is.null(pos)) cat("\n", nrow(pts), " selected")
         } 
         pt.on <- pts
         pt.on <- unique(pt.on)
         colnames(pt.on) <- c("t", "frq")
         pt.on[, "frq"] <- pt.on[, "frq"] + min(which.frq.bins) - 1
 
         # Select "off" points
         i <- 0
	 if(grepl('[Xx]11', .Device)) {
           cat("\nSelect \"off\" points with left mouse click. When done, right click.\n")
	 } else {
           cat("\nSelect \"off\" points with left mouse click. When done, press \'ESC\'.\n")
	 }
         pts <- NULL
         pos <- NULL
         while(!is.null(pos)|i == 0) {
            i <- i + 1
            pos <- locator(n=1)
            if(!is.null(pos)) pos <- lapply(pos, round)
            if(!is.null(pos)) pos$y <- pos$y - min(which.frq.bins) + 1
            pts <- rbind(pts, as.numeric(pos))
            off.mat[pts[, 2:1, drop=FALSE]] <- 1
            # Add points to plot
            image(x=which.t.bins, y=which.frq.bins, t(off.mat), col=c("transparent", sel.col[2]), zlim=c(0, 1), add=TRUE)
            if(!is.null(pos)) cat("\n", nrow(pts), " selected")
         } 
         pt.off <- pts
         pt.off <- unique(pt.off)
         colnames(pt.off) <- c("t", "frq")
         pt.off[, "frq"] <- pt.off[, "frq"] + min(which.frq.bins) - 1

      } else if(select%in%c("rect", "rectangle")) {
         # Rectangular selection
         # Plot over legend
         image(x=which.t.bins, y=which.frq.bins, t(mat3), col=quat.col, zlim=c(0, 3), add=TRUE)
         cat("\nRectangular selection\n")
         i <- 0
         pos1 <- pt.on <- NULL
         while(!is.null(pos1)|i == 0) {
            # On cells first
	    if(grepl('[Xx]11', .Device)) {
              cat("\nSelect upper left corner of \"on\" rectangle with a left click.\nRight click to continue.\n")
	    } else {
              cat("\nSelect upper left corner of \"on\" rectangle with a left click.\nPress \'ESC\' to continue.\n")
	    }
            i <- i + 1
            pos1 <- locator(n=1)
            points(pos1$x, pos1$y, pch=22, cex=0.5, col="red", bg="red")
            if(!is.null(pos1)) {
               cat("\nSelect lower right corner of \"on\" rectangle with a left click.\n")
               pos2 <- locator(n=1)
               points(pos2$x, pos2$y, pch=22, cex=0.5, col="red", bg="red")

               # First find positions within the matrix that are within the rectangle
               frq.in.rect <- which.frq.bins<pos1$y & which.frq.bins>pos2$y
               x.in.rect <- which.t.bins>pos1$x & which.t.bins<pos2$x

               # Set cells that meet criteria to 1 in on.mat
               temp.mat <- on.mat
               temp.mat[frq.in.rect, x.in.rect] <- mat3[frq.in.rect, x.in.rect]
               if(dens<1) temp.mat[temp.mat == 3] <- sample(c(3, 0), sum(temp.mat == 3), TRUE, c(dens, 1-dens))
               on.mat[temp.mat == 3] <- 1
               
               # Then find locations of all overlapped cells within rectangle
               pts <- which(on.mat == 1, arr.ind=TRUE)
               pts <- pts[, 2:1]
               colnames(pts) <- c("t", "frq")
               # Next line should not be necessary here
               pts[, "t"] <- pts[, "t"] + min(which.t.bins) - 1
               pts[, "frq"] <- pts[, "frq"] + min(which.frq.bins) - 1
               pt.on <- rbind(pt.on, pts)

               # Plot over rectangle corners
               image(x=which.t.bins, y=which.frq.bins, t(mat3), col=quat.col, zlim=c(0, 3), add=TRUE)
               # Add points to plot
               image(x=which.t.bins, y=which.frq.bins, t(on.mat), col=c("transparent", sel.col[1]), zlim=c(0, 1), add=TRUE)
            }
            pt.on <- unique(pt.on)
         }
         # Off cells
         i <- 0
         pos1 <- pt.off <- NULL
         while(!is.null(pos1)|i == 0) {
	    if(grepl('[Xx]11', .Device)) {
              cat("\nSelect upper left corner of \"off\" rectangle with a left click.\nRight click to continue.\n")
	    } else {
              cat("\nSelect upper left corner of \"off\" rectangle with a left click.\nPress \'ESC\' to continue.\n")
	    }
            i <- i + 1
            pos1 <- locator(n=1)
            points(pos1$x, pos1$y, pch=22, cex=0.5, col="red", bg="red")
            if(!is.null(pos1)) {
               cat("\nSelect lower right corner of \"off\" rectangle with a left click.\n")
               pos2 <- locator(n=1)
               points(pos2$x, pos2$y, pch=22, cex=0.5, col="red", bg="red")

               # First find positions within the matrix that are within the rectangle
               frq.in.rect <- which.frq.bins<pos1$y & which.frq.bins>pos2$y
               x.in.rect <- which.t.bins>pos1$x & which.t.bins<pos2$x
               
               # Find cells with high binary cells within buffer distance--these will have FALSE in buff.amp
               # This loop is a bit slow
               buff.amp <- matrix(TRUE, nrow=n.frq.bins, ncol=n.t.bins)
               if(buffer>0) {
                  for(i in 1:n.frq.bins) {
                     for(j in 1:n.t.bins) {
                        # Next line looks for cells with 1, 2, or 3 anywhere within a matrix 1 + 2*buffer square, including the center cell [i, j]
                        # The min expression prevents subscripts < 0
                        buff.amp[i, j] <- (sum(mat3[i+ -min(buffer, i-1):min(buffer, n.frq.bins-i), j+ -min(buffer, j-1):min(buffer, n.t.bins-j)]%in%1:3) == 0)
                     }
                  }
               }

               # Set cells that meet all criteria to 1 in off.mat
               temp.mat <- off.mat + 1
               temp.mat[frq.in.rect, x.in.rect] <- mat3[frq.in.rect, x.in.rect]
               if(dens<1) temp.mat[temp.mat == 0] <- sample(c(1, 0), sum(temp.mat == 0), TRUE, c(1-dens, dens)) # Here a 1 in temp.mat mean cell is not "off" cell, so prob values reversed from above
               off.mat[temp.mat == 0 & buff.amp] <- 1 
               
               # Then find locations of cells
               pts <- which(off.mat == 1, arr.ind=TRUE)
               pts <- pts[, 2:1]
               colnames(pts) <- c("t", "frq")
               pts[, "frq"] <- pts[, "frq"] + min(which.frq.bins) - 1
               pt.off <- rbind(pt.off, pts)

               # Plot over rectangle corners
               image(x=which.t.bins, y=which.frq.bins, t(mat3), col=quat.col, zlim=c(0, 3), add=TRUE)
               # Add points to plot
               image(x=which.t.bins, y=which.frq.bins, t(on.mat), col=c("transparent", sel.col[1]), zlim=c(0, 1), add=TRUE)
               image(x=which.t.bins, y=which.frq.bins, t(off.mat), col=c("transparent", sel.col[2]), zlim=c(0, 1), add=TRUE)
            }
            pt.off <- unique(pt.off)
         } 
      } else if(select %in% c("auto", "automatic")) {
         # Automatic point selection
         cat("\nAutomatic point selection.\n")

         # On cells first
         # Set cells that meet criteria to 1 in on.mat, incorporating random selection based on dens
         on.mat[mat3 == 3] <- sample(c(1, 0), sum(mat3 == 3), TRUE, c(dens, 1-dens))

         # Then find locations of all high cells 
         pts <- which(on.mat == 1, arr.ind=TRUE)
         pts <- pts[, 2:1]
         colnames(pts) <- c("t", "frq")
         pts[, "t"] <- pts[, "t"] + min(which.t.bins) - 1
         pts[, "frq"] <- pts[, "frq"] + min(which.frq.bins) - 1
         pt.on <- pts

         # Off cells next
         # Find cells with high binary cells within buffer distance--these will have FALSE in buff.amp
         # This loop is a bit slow
         buff.amp <- matrix(TRUE, nrow=n.frq.bins, ncol=n.t.bins)
         if(buffer>0) {
            for(i in 1:n.frq.bins) {
               for(j in 1:n.t.bins) {
                  # Next line looks for cells with 1, 2, or 3 anywhere within a matrix 1 + 2*buffer square, including the center cell [i, j]
                  # The min expression prevents subscripts < 0
                  buff.amp[i, j] <- (sum(mat3[i+ -min(buffer, i-1):min(buffer, n.frq.bins-i), j+ -min(buffer, j-1):min(buffer, n.t.bins-j)]%in%1:3) == 0)
               }
            }
         }

         # Set cells that meet criteria to 1 in off.mat
         off.mat[mat3 == 0 & buff.amp] <- sample(c(1, 0), sum(mat3 == 0 & buff.amp), TRUE, c(dens, 1-dens))
               
         # Then find locations of cells
         pts <- which(off.mat == 1, arr.ind=TRUE)
         pts <- pts[, 2:1]
         colnames(pts) <- c("t", "frq")
         pts[, "t"] <- pts[, "t"] + min(which.t.bins) - 1
         pts[, "frq"] <- pts[, "frq"] + min(which.frq.bins) - 1
         pt.off <- pts

         # Add points to plot
         image(x=which.t.bins, y=which.frq.bins, t(on.mat), col=c("transparent", sel.col[1]), zlim=c(0, 1), add=TRUE)
         image(x=which.t.bins, y=which.frq.bins, t(off.mat), col=c("transparent", sel.col[2]), zlim=c(0, 1), add=TRUE)
      }
   } else stop("Expected clip length of 1 or 2, but got ", length(clip))

   # Shift time bins to a minimum of 1
   t.shift <- min(pt.on[, "t"], pt.off[, "t"])
   # Adjust first.t.bin for shift
   first.t.bin <- first.t.bin + (t.shift - 1)*t.step
   pt.on[, "t"] <- pt.on[, "t"] - t.shift + 1
   pt.off[, "t"] <- pt.off[, "t"] - t.shift + 1


   # Note that duration is for the time bins (right of last time bin minus left of first time bin)
   n.t.bins <- diff(range(pt.on[, "t"], pt.off[, "t"]))
   n.frq.bins <- diff(range(pt.on[, "frq"], pt.off[, "frq"]))
   duration <- n.t.bins*t.step
   frq.lim <- range(pt.on[, "frq"], pt.off[, "frq"])*frq.step

   template <- list(new("binTemplate", clip.path=clip.path, samp.rate=samp.rate, pt.on=pt.on, pt.off=pt.off, t.step=t.step, frq.step=frq.step, n.t.bins=as.integer(n.t.bins), first.t.bin=first.t.bin, n.frq.bins=as.integer(n.frq.bins), duration=duration, frq.lim=frq.lim, wl=as.integer(wl), ovlp=as.integer(ovlp), wn=wn, score.cutoff=score.cutoff, comment=comment))
   names(template) <- name
   template <- new("binTemplateList", templates=template)

   par(oldpar)
   cat("\nDone.\n")
   return(template)
}
