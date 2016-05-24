  HeatMapPlot <- function(V, row.names = FALSE, col.labels, col.classes, col.names = FALSE, main = " ", xlab=" ", ylab=" ") {

       n.rows <- length(V[,1])
       n.cols <- length(V[1,])
       row.mean <- apply(V, MARGIN=1, FUN=mean)
       row.sd <- apply(V, MARGIN=1, FUN=sd)
       row.n <- length(V[,1])
       for (i in 1:n.rows) {
	   if (row.sd[i] == 0) {
    	       V[i,] <- 0
           } else {
	       V[i,] <- (V[i,] - row.mean[i])/(0.5 * row.sd[i])
           }
           V[i,] <- ifelse(V[i,] < -6, -6, V[i,])
           V[i,] <- ifelse(V[i,] > 6, 6, V[i,])
        }

        mycol <- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000") # blue-pinkogram colors. The first and last are the colors to indicate the class vector (phenotype). This is the 1998-vintage, pre-miR cluster, original pinkogram color map

        mid.range.V <- mean(range(V)) - 0.1
        heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
        heatm[1:n.rows,] <- V[seq(n.rows, 1, -1),]
        heatm[n.rows + 1,] <- ifelse(col.labels == 0, 7, -7)
        image(1:n.cols, 1:(n.rows + 1), t(heatm), col=mycol, axes=FALSE, main=main, xlab= xlab, ylab=ylab)

        if (length(row.names) > 1) {
            numC <- nchar(row.names)
            size.row.char <- 35/(n.rows + 5)
            size.col.char <- 25/(n.cols + 5)
            maxl <- floor(n.rows/1.6)
            for (i in 1:n.rows) {
               row.names[i] <- substr(row.names[i], 1, maxl)
            }
            row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
            axis(2, at=1:(n.rows + 1), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=2, line=-1)
        }

        if (length(col.names) > 1) {
           axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
        }

        C <- split(col.labels, col.labels)
        class1.size <- length(C[[1]])
        class2.size <- length(C[[2]])
        axis(3, at=c(floor(class1.size/2),class1.size + floor(class2.size/2)), labels=col.classes, tick=FALSE, las = 1, cex.axis=1.25, font.axis=2, line=-1)

}

  
  #Results<-MsReport(MsNAME="KEGG_ERBB_SIGNALING_PATHWAY", input.ds, input.cls, weighted.score.type = 1) 
 #miRlist<-Results[[2]]
 PlotHeatMap <-function(miRlist,input.ds,input.cls){ 
      content <- input.ds
      content <- content[-(1:2)]
      col.names <- noquote(unlist(strsplit(content[1], "\t")))
      col.names <- col.names[c(-1, -2)]
      num.cols <- length(col.names)
      content <- content[-1]
      num.lines <- length(content)


      row.nam <- vector(length=num.lines, mode="character")
      row.des <- vector(length=num.lines, mode="character")
      m <- matrix(0, nrow=num.lines, ncol=num.cols)

      for (i in 1:num.lines) {
         line.list <- noquote(unlist(strsplit(content[i], "\t")))
         row.nam[i] <- noquote(line.list[1])
         row.des[i] <- noquote(line.list[2])
         line.list <- line.list[c(-1, -2)]
         for (j in 1:length(line.list)) {
            m[i, j] <- as.numeric(line.list[j])
         }
      }
      dataset <- data.frame(m)
      names(dataset) <- col.names
      row.names(dataset) <- row.nam
     
  
  miR.labels <- row.names(dataset)
  sample.names <- names(dataset)
  A <- data.matrix(dataset)
    cols <- length(A[1,])
  rows <- length(A[,1])
#############################################################
   cls.cont <- input.cls
  num.lines <- length(cls.cont)
      class.list <- unlist(strsplit(cls.cont[[3]], " "))
      s <- length(class.list)
      t <- table(class.list)
      l <- length(t)
      class.phen <- vector(length=l, mode="character")
      phen.label <- vector(length=l, mode="numeric")
      class.labels <- vector(length=s, mode="numeric")
      for (i in 1:l) {
         class.phen[i] <- noquote(names(t)[i])
         phen.label[i] <- i - 1
      }
      for (i in 1:s) {
         for (j in 1:l) {
             if (class.list[i] == class.phen[j]) {
                class.labels[i] <- phen.label[j]
             }
         }
      }
     
     phen1 <- class.phen[1]
     phen2 <- class.phen[2]

 col.index <- order(class.labels, decreasing=FALSE)
 class.labels <- class.labels[col.index]
 sample.names <- sample.names[col.index]
 for (j in 1:length(A[,1])) {
    A[j, ] <- A[j, col.index]
 }
 names(A) <- sample.names
 #    
 ###################################################
          miResult<-miRlist
         size.M<-as.numeric(miResult[5])
		 obs.index<-as.numeric(miResult[10][[1]])
		 obs.miR.labels<-as.character(miResult[11][[1]])
		  Obs.indicator<-as.numeric(miResult[7][[1]])
		
            kk <- 1
            pinko <- matrix(0, nrow = size.M, ncol = cols)
            pinko.miR.names <- vector(length = size.M, mode = "character")
            for (k in 1:rows) {
               if (Obs.indicator[k] == 1) {
                  pinko[kk,] <- A[obs.index[k],]
                  pinko.miR.names[kk] <-  obs.miR.labels[k]
                  kk <- kk + 1
               }
            }
			HeatMapPlot(V = pinko, row.names = pinko.miR.names, col.labels = class.labels, col.classes = class.phen, col.names = sample.names, main =" Heat Map for MiRs in MiR Set", xlab=" ", ylab=" ")
         
		}