CityPlot <- function(controlfile="", pdfmode="n", shownames=T) {
   cat("CityPlot\n\n")
   # read controlfile
   if (controlfile == "") controlfile <- file.choose()
   currentdir <- getwd()
   setwd(dirname(controlfile))
   cf <- read.table(file=controlfile, header=F, sep="\t", as.is=T)
   print(cf)
   if (dim(cf)[2] == 1) { files <- unlist(cf) }
   if (dim(cf)[2] > 1)  { files <- unique(unlist(cf[,1:2])) }

   # output as pdf?
   if (pdfmode == "") { 
	   cat("Output as pdf (y/n)")
	   pdfmode <- readLines(n=1)
	}

	# determine size of plotting region
	wtotal <- 0
	htotal <- 0
	pos <- vector()
	for (filename in files) {
	   label <- tail(unlist(strsplit(filename, "\\", fixed=T)),n=1)
	   label <- sub(".csv","",label)
         # determine csv-separator
         testline <- readLines(filename,n=2)[2]
         s1 <- length(unlist(strsplit(testline,";")))
         s2 <- length(unlist(strsplit(testline,",")))
         s3 <- length(unlist(strsplit(testline,"\t")))
         separator <- ";"
         if (s3 > s1 && s3 > s2) { separator <- "\t" }
         if (s2 > s1 && s2 > s3) { separator <- ","  }
	   # width and height of each table
	   df <- read.table(file=filename, skip=1, fill=T, sep=separator, as.is=T, row.names = NULL)
	   h <- dim(df)[1]
	   if (h > htotal) htotal <- h
	   w <- dim(df)[2]
	   wtotal <- wtotal + 3 + w
	   cat(label," ncol=",w," nrow=",h,"\n")
	}

	# cat("htotal=",htotal," wtotal=",wtotal,"\n")
	if (pdfmode == "y" || pdfmode == "Y") {
	   cat("Creating pdf ... ")
	   pdf( paste(controlfile,".pdf",sep="") ) 
	}

	split.screen(matrix(c(0,1,0,1), ncol=4))
	split.screen(matrix(c(0,1,0.2,1), ncol=4))
	screen(2)
	par(mar=c(0,0,0,0))
	plot(c(0,wtotal), c(1,htotal*1.15), xlab="", ylab="", log="y", axes=F)

	xoffset <- 0
	for (filename in files) {
	   overall_type <- matrix()
	   color <- matrix()
	   label <- tail(unlist(strsplit(filename, "\\", fixed=T)),n=1)
	   label <- sub(".csv","",label)
         # determine csv-separator
         testline <- readLines(filename,n=2)[2]
         s1 <- length(unlist(strsplit(testline,";")))
         s2 <- length(unlist(strsplit(testline,",")))
         s3 <- length(unlist(strsplit(testline,"\t")))
         separator <- ";"
         if (s3 > s1 && s3 > s2) { separator <- "\t" }
         if (s2 > s1 && s2 > s3) { separator <- ","  }
	   # read table data
	   df <- read.table(file=filename, skip=1, fill=T, sep=separator, as.is=T, row.names = NULL)
	   h <- dim(df)[1]
	   w <- dim(df)[2]
         select <- rep(F,h)
         for (i in 1:min(1000,h)) { select[i] <- T }
         if (h>1000) {
             for (i in seq(1000,min(10000,h), by=10) ) {
                select[i] <- T
             }
         }
         if (h>10000) {
             for (i in seq(10000,min(100000,h), by=100) ) {
                select[i] <- T
             }
         }
         if (h>100000) {
             for (i in seq(100000,min(1000000,h), by=1000) ) {
                select[i] <- T
             }
         }
         index <- 1:h
         index <- index[select]
         dselect <- df[select,]
	   dselect[is.na(dselect)] <- ""
         dselect[dselect == "?"] <- ""
         dselect[dselect == "\\N"] <- ""
	   for (i in 1:w) {
		  dvec <- dselect[,i]
		  empty <- dvec == ""

		  # time items
		  timevalues  <- suppressWarnings(strptime(dvec, "%H:%M:%S"))
		  timevalues2 <- suppressWarnings(strptime(dvec, "%H:%M"))
		  itemtype_time <- sum(is.na(timevalues[!empty])) == 0 || sum(is.na(timevalues2[!empty])) == 0
		  if (itemtype_time) { overall_type[i] <- "time"; color[i] <- "yellow"; next }

		  # date items
		  datevalues <- suppressWarnings(strptime(dvec, "%d.%m.%y"))
		  itemtype_date <- sum(is.na(datevalues[!empty])) == 0
		  if (itemtype_date) { overall_type[i] <- "date"; color[i] <- "gold"; next }

		  # numeric / categoric items
		  tmp <- sub(",",".",dvec)
		  numbers <- suppressWarnings(as.numeric(tmp))
		  if (sum(is.na(numbers[!empty])) == 0) {
			 factors <- as.factor(dvec[!empty])
			 levelcount <- length(levels(factors))
			 if ( length(factors) < 4 * levelcount || levelcount > 10 ) 
				{ overall_type[i] <- "numeric"; color[i] <- "blue" }
			 else { overall_type[i] <- "categoric"; color[i] <- "green" }
			 next
		  }

		  # categoric items
		  factors <- as.factor(dvec[!empty])
		  levelcount <- length(levels(factors))
		  # cat("item",i,"levelcount=",levelcount,"\n")
		  itemtype_categoric <- length(factors) > (4 * levelcount) && levelcount <= 10
		  if (itemtype_categoric) { overall_type[i] <- "categoric"; color[i] <- "green" }
		  else { overall_type[i] <- "other"; color[i] <- "red" }

	   }

	   # for (i in 1:w) { cat(colnames(df)[i], "\t", overall_type[i], "\t", paste(df[1:5,i]), " ...\n") }

	   rect(xoffset,1, xoffset + w, h+1, lwd=5)
	   if (shownames) text(xoffset + w/2, h+1, label, pos=3)
	   pos <- append(pos, xoffset + w/2)

         size <- 1
         for (i in 1:length(index)) {
              if (index[i] > 1000)   size <- 10
              if (index[i] > 10000)  size <- 100
              if (index[i] > 100000)  size <- 1000

		  for (k in 1:w) {
                   if (dselect[i,k] != "")
			 rect((xoffset + k-0.9),(index[i] + 0.1 ),(xoffset + k-0.1),(index[i]+size - 0.1), col=color[k], border=NA )
		  }
	   }
	   xoffset <- xoffset + w + 3
	}

	# plot relations
   if (dim(cf)[2] > 2) { 
      screen(1)
      split.screen(matrix(c(0,1,0,0.23), ncol=4))
      par(mar=c(0,0,0,0))
      cf2 <- cf
      for (i in 1:length(files)) {
         cf2[cf == files[i]] <- pos[i]
      }
      nlev <- 1
      if (nrow(cf2) > 1) {
         for (i in 2:nrow(cf2)) {
	      if (as.numeric(cf2[i-1,2]) > as.numeric(cf2[i,1])) nlev <- nlev + 1 
	   }
      }
      plot(c(0,wtotal), c(0,nlev +1), xlab="", ylab="", axes=F)
      
      y <- nlev
      for (i in 1:nrow(cf2)) {
         label <- unlist(strsplit(cf2[i,3], ":", fixed=T))
         a  <- as.numeric(cf2[i,1])
         b  <- as.numeric(cf2[i,2])
         m  <- (b+a) / 2
         a2 <- m - wtotal/50
         b2 <- m + wtotal/50
         if (i>1) if (as.numeric(cf2[i-1,2]) > a) y <- y - 1
         lines(c(a,a2), c(y,y), lwd=3)
         lines(c(a2,m), c(y,y+0.3), lwd=3)
         lines(c(a2,m), c(y,y-0.3), lwd=3)
         text(a2,y+0.3,label[1],adj=1.1)
         text(b2,y+0.3,label[2],adj=-0.3)
         lines(c(b2,m), c(y,y+0.3), lwd=3)
         lines(c(b2,m), c(y,y-0.3), lwd=3)
         lines(c(b2,b), c(y,y), lwd=3) 
         lines(c(a,a), c(y,nlev+1), lwd=3) 
         lines(c(b,b), c(y,nlev+1), lwd=3)
      }
   }   
   if (pdfmode != "y" && pdfmode != "Y") readline(prompt="set focus to console and press enter")
   else cat("completed.\n")
   dev.off()
   setwd(currentdir)
}