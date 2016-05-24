plot.SPF.BinBin <- function(x, Type="All.Histograms", Specific.Pi="r_0_0", Col="grey", 
                            Box.Plot.Outliers=FALSE, Legend.Pos="topleft", Legend.Cex=1, ...){
                            
  Object <- x 

  if (missing(Col)) {Col = "grey"}
  
  if (Type=="All.Histograms"){
    
    #no mono
    if ((length(unique(Object$r_min1_min1))> 1) & (length(unique(Object$r_0_min1))> 1) &
      (length(unique(Object$r_min1_0))> 1)){
    
    plot(0:100, 0:100, axes=F, xlab="", ylab="", type="n", ..., xlim=c(0, 1))  
    par(mfrow=c(3, 3), mar = c(4.5, 7, 4, 1), oma=rep(0, times=4))  
    
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      hist(Object$r_min1_min1, main="", col=Col, xlim=c(0, 1),
           xlab=expression(r(-1,-1)), cex.lab=1.3)}
    if (length(unique(Object$r_min1_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    #mtext(side = 3, expression(paste(Delta, "T = -1")), cex=1.5, padj = -1.6)   
    #mtext(side = 2, expression(paste(Delta, "S = -1")), cex=1.5, padj = -3.6)   
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      hist(Object$r_0_min1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(0,-1)))}
    if (length(unique(Object$r_0_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    #mtext(side = 3, expression(paste(Delta, "T = 0")), cex=1.5, padj = -1.6)   
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      hist(Object$r_1_min1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(1,-1)))}
    if (length(unique(Object$r_1_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    #mtext(side = 3, expression(paste(Delta, "T = 1")), cex=1.5, padj = -1.6)   
    
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      hist(Object$r_min1_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(-1,0)))}
    if (length(unique(Object$r_min1_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    #mtext(side = 2, expression(paste(Delta, "S = 0")), cex=1.5, padj = -3.6)   
        
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      hist(Object$r_0_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(0,0)))}
    if (length(unique(Object$r_0_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      hist(Object$r_1_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(1,0)))}
    if (length(unique(Object$r_1_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      hist(Object$r_min1_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(-1,1)))}
    if (length(unique(Object$r_min1_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    #mtext(side = 2, expression(paste(Delta, "S = 1")), cex=1.5, padj = -3.6)   
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      hist(Object$r_0_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(0,1)))}
    if (length(unique(Object$r_0_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      hist(Object$r_1_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(1,1)))}
    if (length(unique(Object$r_1_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    }
    
    # mono S and T
    if ((length(unique(Object$r_min1_min1))==1) & (length(unique(Object$r_0_min1))==1) &
          (length(unique(Object$r_min1_0))==1)){
      
      plot(0:100, 0:100, axes=F, xlab="", ylab="", type="n", ..., xlim=c(0, 1))  
      par(mfrow=c(2, 2), mar = c(4.5, 7, 4, 1), oma=rep(0, times=4))  
      
      # T = 0, S = 0
      if (length(unique(Object$r_0_0)) > 1){
        hist(Object$r_0_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(0,0)))}
      if (length(unique(Object$r_0_0)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 1, S = 0
      if (length(unique(Object$r_1_0)) > 1){
        hist(Object$r_1_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(1,0)))}
      if (length(unique(Object$r_1_0)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 0, S = 1
      if (length(unique(Object$r_0_1)) > 1){
        hist(Object$r_0_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(0,1)))}
      if (length(unique(Object$r_0_1)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 1, S = 1
      if (length(unique(Object$r_1_1)) > 1){
        hist(Object$r_1_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(1,1)))}
      if (length(unique(Object$r_1_1)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
    }
    
    # mono S
    if ((length(unique(Object$r_min1_min1))==1) & (length(unique(Object$r_0_min1))== 1) &
          (length(unique(Object$r_min1_0))> 1)){
      
      plot(0:100, 0:100, axes=F, xlab="", ylab="", type="n", ..., xlim=c(0, 1))  
      par(mfrow=c(2, 3), mar = c(4.5, 7, 4, 1), oma=rep(0, times=4))  
      
      # T = -1, S = 0
      if (length(unique(Object$r_min1_0)) > 1){
        hist(Object$r_min1_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(-1,0)))}
      if (length(unique(Object$r_min1_0)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 0, S = 0
      if (length(unique(Object$r_0_0)) > 1){
        hist(Object$r_0_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(0,0)))}
      if (length(unique(Object$r_0_0)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 1, S = 0
      if (length(unique(Object$r_1_0)) > 1){
        hist(Object$r_1_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(1,0)))}
      if (length(unique(Object$r_1_0)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = -1, S = 1
      if (length(unique(Object$r_min1_1)) > 1){
        hist(Object$r_min1_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(-1,1)))}
      if (length(unique(Object$r_min1_1)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 0, S = 1
      if (length(unique(Object$r_0_1)) > 1){
        hist(Object$r_0_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(0,1)))}
      if (length(unique(Object$r_0_1)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 1, S = 1
      if (length(unique(Object$r_1_1)) > 1){
        hist(Object$r_1_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(1,1)))}
      if (length(unique(Object$r_1_1)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
    }
    
    
    # mono T
    if ((length(unique(Object$r_min1_min1))==1) & (length(unique(Object$r_0_min1))>1) &
          (length(unique(Object$r_min1_0))==1)){
      
      plot(0:100, 0:100, axes=F, xlab="", ylab="", type="n", ..., xlim=c(0, 1))  
      par(mfrow=c(3, 2), mar = c(4.5, 7, 4, 1), oma=rep(0, times=4))  
      
      # T = 0, S = -1
      if (length(unique(Object$r_0_min1)) > 1){
        hist(Object$r_0_min1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(0,-1)))}
      if (length(unique(Object$r_0_min1)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 1, S = -1
      if (length(unique(Object$r_1_min1)) > 1){
        hist(Object$r_1_min1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(1,-1)))}
      if (length(unique(Object$r_1_min1)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 0, S = 0
      if (length(unique(Object$r_0_0)) > 1){
        hist(Object$r_0_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(0,0)))}
      if (length(unique(Object$r_0_0)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 1, S = 0
      if (length(unique(Object$r_1_0)) > 1){
        hist(Object$r_1_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(1,0)))}
      if (length(unique(Object$r_1_0)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 0, S = 1
      if (length(unique(Object$r_0_1)) > 1){
        hist(Object$r_0_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(0,1)))}
      if (length(unique(Object$r_0_1)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
      # T = 1, S = 1
      if (length(unique(Object$r_1_1)) > 1){
        hist(Object$r_1_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
             xlab=expression(r(1,1)))}
      if (length(unique(Object$r_1_1)) <= 1){
        plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
      
    }   
    
    
    
    
    
    
    par(mfrow=c(1, 1), c(5, 4, 4, 2) + 0.1)
  }
  
  if (Type=="All.Densities"){
    
    plot(0:100, 0:100, axes=F, xlab="", ylab="", type="n", ..., xlim=c(0, 1))  #LS!
    par(mfrow=c(3, 3), mar = c(4.5, 7, 4, 1), oma=rep(0, times=4), xpd=FALSE)  #LS!
    
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      plot(density(Object$r_min1_min1, na.rm=T), main=" ", col=Col, cex.lab=1.3, ..., xlim=c(0, 1),
           xlab=expression(r(-1,-1)))
    }
    if (length(unique(Object$r_min1_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    mtext(side = 3, expression(paste(Delta, "T = -1")), cex=1.5, padj = -1.6)   
    mtext(side = 2, expression(paste(Delta, "S = -1")), cex=1.5, padj = -3.6)   
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      plot(density(Object$r_0_min1, na.rm=T), main=" ", col=Col, cex.lab=1.3, ...,xlim=c(0, 1),
           xlab=expression(r(0,-1)))}
    if (length(unique(Object$r_0_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    
    mtext(side = 3, expression(paste(Delta, "T = 0")), cex=1.5, padj = -1.6)   
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      plot(density(Object$r_1_min1, na.rm=T), main=" ", col=Col, cex.lab=1.3, ...,xlim=c(0, 1),
           xlab=expression(r(1,-1)))}
    if (length(unique(Object$r_1_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    
    mtext(side = 3, expression(paste(Delta, "T = 1")), cex=1.5, padj = -1.6)   
    
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      plot(density(Object$r_min1_0, na.rm=T), main=" ", col=Col, cex.lab=1.3, ...,xlim=c(0, 1),
           xlab=expression(r(-1,0)))}
    if (length(unique(Object$r_min1_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    
    mtext(side = 2, expression(paste(Delta, "S = 0")), cex=1.5, padj = -3.6)   
    
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      plot(density(Object$r_0_0, na.rm=T), main=" ", col=Col, cex.lab=1.3, ...,xlim=c(0, 1),
           xlab=expression(r(0,0)))}
    if (length(unique(Object$r_0_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      plot(density(Object$r_1_0, na.rm=T), main=" ", col=Col, cex.lab=1.3, ...,xlim=c(0, 1),
           xlab=expression(r(1,0)))}
    if (length(unique(Object$r_1_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      plot(density(Object$r_min1_1, na.rm=T), main=" ", col=Col, cex.lab=1.3, ...,xlim=c(0, 1),
           xlab=expression(r(-1,1)))}
    if (length(unique(Object$r_min1_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    mtext(side = 2, expression(paste(Delta, "S = 1")), cex=1.5, padj = -3.6)   
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      plot(density(Object$r_0_1, na.rm=T), main=" ", col=Col, cex.lab=1.3, ...,xlim=c(0, 1),
           xlab=expression(r(0,1)))}
    if (length(unique(Object$r_0_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      plot(density(Object$r_1_1, na.rm=T), main=" ", col=Col, cex.lab=1.3, ...,xlim=c(0, 1),
           xlab=expression(r(1,1)))}
    if (length(unique(Object$r_1_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ", ...)}
    
    par(mfrow=c(1, 1), c(5, 4, 4, 2) + 0.1)
    
  }
  

  if (Type=="Histogram"){
    
    par(mfrow=c(1, 1), c(5, 4, 4, 2) + 0.1)
    
    if (Specific.Pi == "r_min1_min1"){
      # T = -1, S = -1
      if (length(unique(Object$r_min1_min1)) > 1){
        hist(Object$r_min1_min1, main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(-1,-1)))}
      if (length(unique(Object$r_min1_min1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_0_min1"){
      # T = 0, S = -1
      if (length(unique(Object$r_0_min1)) > 1){
        hist(Object$r_0_min1, main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(0,-1)))}
      if (length(unique(Object$r_0_min1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_1_min1"){
      # T = 1, S = -1
      if (length(unique(Object$r_1_min1)) > 1){
        hist(Object$r_1_min1, main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(1,-1)))}
      if (length(unique(Object$r_1_min1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_min1_0"){
      # T = -1, S = 0
      if (length(unique(Object$r_min1_0)) > 1){
        hist(Object$r_min1_0, main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(-1,0)))}
      if (length(unique(Object$r_min1_0)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_0_0"){
      # T = 0, S = 0
      if (length(unique(Object$r_0_0)) > 1){
        hist(Object$r_0_0, main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(0,0)))}
      if (length(unique(Object$r_0_0)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_1_0"){
      # T = 1, S = 0
      if (length(unique(Object$r_1_0)) > 1){
        hist(Object$r_1_0, main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(1,0)))}
      if (length(unique(Object$r_1_0)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_min1_1"){
      # T = -1, S = 1
      if (length(unique(Object$r_min1_1)) > 1){
        hist(Object$r_min1_1, main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(-1,1)))}
      if (length(unique(Object$r_min1_1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_0_1"){
      # T = 0, S = 1
      if (length(unique(Object$r_0_1)) > 1){
        hist(Object$r_0_1, main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(0,1)))}
      if (length(unique(Object$r_0_1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_1_1"){
      # T = 1, S = 1
      if (length(unique(Object$r_1_1)) > 1){
        hist(Object$r_1_1, main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(1,1)))}
      if (length(unique(Object$r_1_1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
  }
  
  if (Type=="Density"){
    
    par(mfrow=c(1, 1), c(5, 4, 4, 2) + 0.1)
    
    if (Specific.Pi == "r_min1_min1"){
      # T = -1, S = -1
      if (length(unique(Object$r_min1_min1)) > 1){
        plot(density(Object$r_min1_min1, na.rm=T), main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(-1,-1)))
      }
      if (length(unique(Object$r_min1_min1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_0_min1"){
      # T = 0, S = -1
      if (length(unique(Object$r_0_min1)) > 1){
        plot(density(Object$r_0_min1, na.rm=T), main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(0,-1)))}
      if (length(unique(Object$r_0_min1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_1_min1"){
      # T = 1, S = -1
      if (length(unique(Object$r_1_min1)) > 1){
        plot(density(Object$r_1_min1, na.rm=T), main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(1,-1)))}
      if (length(unique(Object$r_1_min1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_min1_0"){
      # T = -1, S = 0
      if (length(unique(Object$r_min1_0)) > 1){
        plot(density(Object$r_min1_0, na.rm=T), main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(-1,0)))}
      if (length(unique(Object$r_min1_0)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_0_0"){
      # T = 0, S = 0
      if (length(unique(Object$r_0_0)) > 1){
        plot(density(Object$r_0_0, na.rm=T), main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(0,0)))}
      if (length(unique(Object$r_0_0)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_1_0"){
      # T = 1, S = 0
      if (length(unique(Object$r_1_0)) > 1){
        plot(density(Object$r_1_0, na.rm=T), main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(1,0)))}
      if (length(unique(Object$r_1_0)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_min1_1"){
      # T = -1, S = 1
      if (length(unique(Object$r_min1_1)) > 1){
        plot(density(Object$r_min1_1, na.rm=T), main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(-1,1)))}
      if (length(unique(Object$r_min1_1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_0_1"){
      # T = 0, S = 1
      if (length(unique(Object$r_0_1)) > 1){
        plot(density(Object$r_0_1, na.rm=T), main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(0,1)))}
      if (length(unique(Object$r_0_1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
    
    if (Specific.Pi == "r_1_1"){
      # T = 1, S = 1
      if (length(unique(Object$r_1_1)) > 1){
        plot(density(Object$r_1_1, na.rm=T), main=" ", col=Col, ...,xlim=c(0, 1),
             xlab=expression(r(1,1)))}
      if (length(unique(Object$r_1_1)) <= 1){
        cat("\nNo valid pi values were found. \n")}
    }
  }
  
  if (Type=="Box.Plot"){
    
    par(mfrow=c(1, 1), c(5, 4, 4, 2) + 0.1)
        
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      a <- cbind(Object$r_min1_min1, "a")}
    if (length(unique(Object$r_min1_min1)) <= 1){
      a <- cbind(NA, "a")}
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      b <- cbind(Object$r_0_min1, "b")}
    if (length(unique(Object$r_0_min1)) <= 1){
      b <- cbind(NA, "b")}
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      c <- cbind(Object$r_1_min1, "c")}
    if (length(unique(Object$r_1_min1)) <= 1){
      c <- cbind(NA, "c")}
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      d <- cbind(Object$r_min1_0, "d")}
    if (length(unique(Object$r_min1_0)) <= 1){
      d <- cbind(NA, "d")}
    
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      e <- cbind(Object$r_0_0, "e")}
    if (length(unique(Object$r_0_0)) <= 1){
      e <- cbind(NA, "e")}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      f <- cbind(Object$r_1_0, "f")}
    if (length(unique(Object$r_1_0)) <= 1){
      f <- cbind(NA, "f")}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      g <- cbind(Object$r_min1_1, "g")}
    if (length(unique(Object$r_min1_1)) <= 1){
      g <- cbind(NA, "g")}
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      h <- cbind(Object$r_0_1, "h")}
    if (length(unique(Object$r_0_1)) <= 1){
      h <- cbind(NA, "h")}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      i <- cbind(Object$r_1_1, "i")}
    if (length(unique(Object$r_1_1)) <= 1){
      i <- cbind(NA, "i")}
    
    data <- data.frame(rbind(a, b, c, d, e, f, g, h, i))
    
    as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
    
    boxplot(as.numeric.factor(data$X1) ~ data$X2, col=rep(c(2, 3, 4), times=3), names=rep(c(-1, 0, 1), each=3), 
            outline = Box.Plot.Outliers, xlab=expression(paste(Delta, S)), ...)
    
    abline(v = c(3.5, 6.5), col="blue", lty=3)
    
    legend(Legend.Pos, cex = Legend.Cex, c(expression(paste(Delta, "T=-1")), expression(paste(Delta, "T=0")), expression(paste(Delta, "T=1"))),
           fill = c(2, 3, 4))

  }  

  if (Type=="Lines.Mean"){
    
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      a <- cbind(mean(Object$r_min1_min1), "a")}
    if (length(unique(Object$r_min1_min1)) <= 1){
      a <- cbind(NA, "a")}
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      b <- cbind(mean(Object$r_0_min1), "b")}
    if (length(unique(Object$r_0_min1)) <= 1){
      b <- cbind(NA, "b")}
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      c <- cbind(mean(Object$r_1_min1), "c")}
    if (length(unique(Object$r_1_min1)) <= 1){
      c <- cbind(NA, "c")}
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      d <- cbind(mean(Object$r_min1_0), "d")}
    if (length(unique(Object$r_min1_0)) <= 1){
      d <- cbind(NA, "d")}
    
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      e <- cbind(mean(Object$r_0_0), "e")}
    if (length(unique(Object$r_0_0)) <= 1){
      e <- cbind(NA, "e")}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      f <- cbind(mean(Object$r_1_0), "f")}
    if (length(unique(Object$r_1_0)) <= 1){
      f <- cbind(NA, "f")}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      g <- cbind(mean(Object$r_min1_1), "g")}
    if (length(unique(Object$r_min1_1)) <= 1){
      g <- cbind(NA, "g")}
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      h <- cbind(mean(Object$r_0_1), "h")}
    if (length(unique(Object$r_0_1)) <= 1){
      h <- cbind(NA, "h")}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      i <- cbind(mean(Object$r_1_1), "i")}
    if (length(unique(Object$r_1_1)) <= 1){
      i <- cbind(NA, "i")}
    
    data <- data.frame(rbind(a, b, c, d, e, f, g, h, i))
    
    as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
    
    
    plot(as.numeric.factor(data$X1), col=rep(c(2, 3, 4), times=3), axes=FALSE, type="h", lwd=5,
         xlab=expression(paste(Delta, S)), ylab="Mean", ...)
    axis(1, at = c(1:9), labels = rep(c(-1, 0, 1), each=3))     
    axis(2)
    box()        
    
    legend(Legend.Pos, cex = Legend.Cex, c(expression(paste(Delta, "T=-1")), expression(paste(Delta, "T=0")), 
     expression(paste(Delta, "T=1"))), lty=rep(1, times=3), col=c(2, 3, 4), lwd=c(5, 5, 5))
    
    abline(v = c(3.5, 6.5), col="blue", lty=3)
    
  }

  if (Type=="Lines.Median"){
    
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      a <- cbind(median(Object$r_min1_min1), "a")}
    if (length(unique(Object$r_min1_min1)) <= 1){
      a <- cbind(NA, "a")}
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      b <- cbind(median(Object$r_0_min1), "b")}
    if (length(unique(Object$r_0_min1)) <= 1){
      b <- cbind(NA, "b")}
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      c <- cbind(median(Object$r_1_min1), "c")}
    if (length(unique(Object$r_1_min1)) <= 1){
      c <- cbind(NA, "c")}
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      d <- cbind(median(Object$r_min1_0), "d")}
    if (length(unique(Object$r_min1_0)) <= 1){
      d <- cbind(NA, "d")}
    
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      e <- cbind(median(Object$r_0_0), "e")}
    if (length(unique(Object$r_0_0)) <= 1){
      e <- cbind(NA, "e")}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      f <- cbind(median(Object$r_1_0), "f")}
    if (length(unique(Object$r_1_0)) <= 1){
      f <- cbind(NA, "f")}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      g <- cbind(median(Object$r_min1_1), "g")}
    if (length(unique(Object$r_min1_1)) <= 1){
      g <- cbind(NA, "g")}
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      h <- cbind(median(Object$r_0_1), "h")}
    if (length(unique(Object$r_0_1)) <= 1){
      h <- cbind(NA, "h")}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      i <- cbind(median(Object$r_1_1), "i")}
    if (length(unique(Object$r_1_1)) <= 1){
      i <- cbind(NA, "i")}
    
    data <- data.frame(rbind(a, b, c, d, e, f, g, h, i))
    
    as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
    
    
    plot(as.numeric.factor(data$X1), col=rep(c(2, 3, 4), times=3), axes=FALSE, type="h", lwd=5,
         xlab=expression(paste(Delta, S)), ylab="Median", ...)
    axis(1, at = c(1:9), labels = rep(c(-1, 0, 1), each=3))     
    axis(2)
    box()        
    
    legend(Legend.Pos, cex = Legend.Cex, c(expression(paste(Delta, "T=-1")), expression(paste(Delta, "T=0")), 
    expression(paste(Delta, "T=1"))), lty=rep(1, times=3), col=c(2, 3, 4), lwd=c(5, 5, 5))
    
    abline(v = c(3.5, 6.5), col="blue", lty=3)
  }
  
  if (Type=="Lines.Mode"){
    
    mode <- function(data) {
      x <- data
      if (unique(x[1])!=0){
        z <- density(x)
        mode_val <- z$x[which.max(z$y)]
        if (mode_val < 0){mode_val <- c(0)}
      }
      if (unique(x[1])==0){
        model_val <- c(0)
      }  
      fit <- list(mode_val= mode_val)  
    }
    
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      a <- cbind(mode(Object$r_min1_min1)$mode_val, "a")}
    if (length(unique(Object$r_min1_min1)) <= 1){
      a <- cbind(NA, "a")}
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      b <- cbind(mode(Object$r_0_min1)$mode_val, "b")}
    if (length(unique(Object$r_0_min1)) <= 1){
      b <- cbind(NA, "b")}
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      c <- cbind(mode(Object$r_1_min1)$mode_val, "c")}
    if (length(unique(Object$r_1_min1)) <= 1){
      c <- cbind(NA, "c")}
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      d <- cbind(mode(Object$r_min1_0)$mode_val, "d")}
    if (length(unique(Object$r_min1_0)) <= 1){
      d <- cbind(NA, "d")}
    
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      e <- cbind(mode(Object$r_0_0)$mode_val, "e")}
    if (length(unique(Object$r_0_0)) <= 1){
      e <- cbind(NA, "e")}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      f <- cbind(mode(Object$r_1_0)$mode_val, "f")}
    if (length(unique(Object$r_1_0)) <= 1){
      f <- cbind(NA, "f")}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      g <- cbind(mode(Object$r_min1_1)$mode_val, "g")}
    if (length(unique(Object$r_min1_1)) <= 1){
      g <- cbind(NA, "g")}
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      h <- cbind(mode(Object$r_0_1)$mode_val, "h")}
    if (length(unique(Object$r_0_1)) <= 1){
      h <- cbind(NA, "h")}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      i <- cbind(mode(Object$r_1_1)$mode_val, "i")}
    if (length(unique(Object$r_1_1)) <= 1){
      i <- cbind(NA, "i")}
    
    data <- 
      data.frame(rbind(a, b, c, d, e, f, g, h, i))
    
    
    as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
    
    plot(as.numeric.factor(data$X1), col=rep(c(2, 3, 4), times=3), axes=FALSE, type="h", ...,
         xlab=expression(paste(Delta, S)), ylab="Mode", lwd=5)
    axis(1, at = c(1:9), labels = rep(c(-1, 0, 1), each=3))     
    axis(2)
    box()        
    
    legend(Legend.Pos, cex = Legend.Cex, c(expression(paste(Delta, "T=-1")), expression(paste(Delta, "T=0")), 
                                           expression(paste(Delta, "T=1"))), lty=rep(1, times=3), col=c(2, 3, 4), lwd=c(5, 5, 5))
    
    abline(v = c(3.5, 6.5), col="blue", lty=3)
  }
  
  if (Type=="3D.Mean"){
    
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      a <- cbind(mean(Object$r_min1_min1), -1, -1)}
    if (length(unique(Object$r_min1_min1)) <= 1){
      a <- cbind(NA, -1, -1)}
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      b <- cbind(mean(Object$r_0_min1), 0, -1)}
    if (length(unique(Object$r_0_min1)) <= 1){
      b <- cbind(NA, 0, -1)}
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      c <- cbind(mean(Object$r_1_min1), 1, -1)}
    if (length(unique(Object$r_1_min1)) <= 1){
      c <- cbind(NA, 1, -1)}
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      d <- cbind(mean(Object$r_min1_0), -1, 0)}
    if (length(unique(Object$r_min1_0)) <= 1){
      d <- cbind(NA, -1, 0)}
    
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      e <- cbind(mean(Object$r_0_0), 0, 0)}
    if (length(unique(Object$r_0_0)) <= 1){
      e <- cbind(NA, 0, 0)}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      f <- cbind(mean(Object$r_1_0), 1, 0)}
    if (length(unique(Object$r_1_0)) <= 1){
      f <- cbind(NA, 1, 0)}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      g <- cbind(mean(Object$r_min1_1), -1, 1)}
    if (length(unique(Object$r_min1_1)) <= 1){
      g <- cbind(NA, -1, 1)}
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      h <- cbind(mean(Object$r_0_1), 0, 1)}
    if (length(unique(Object$r_0_1)) <= 1){
      h <- cbind(NA, 0, 1)}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      i <- cbind(mean(Object$r_1_1), 1, 1)}
    if (length(unique(Object$r_1_1)) <= 1){
      i <- cbind(NA, 1, 1)}
    
    data <- data.frame(rbind(a, b, c, d, e, f, g, h, i))
    names(data) <- c("Y", "Delta_T", "Delta_S")
    
    temp <- lattice::cloud(Y ~ as.factor(Delta_S) + as.factor(Delta_T), data, panel.3d.cloud=latticeExtra::panel.3dbars, 
                                xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1),  
                                par.settings = list(axis.line = list(col = "transparent")), 
                                xlab=expression(paste(Delta, "S")), ylab=expression(paste(Delta, T)), 
                                zlab="Mean", col.facet=rep(c(2:4), each=3))
    plot(temp) 
  }  

  if (Type=="3D.Median"){
    
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      a <- cbind(median(Object$r_min1_min1), -1, -1)}
    if (length(unique(Object$r_min1_min1)) <= 1){
      a <- cbind(NA, -1, -1)}
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      b <- cbind(median(Object$r_0_min1), 0, -1)}
    if (length(unique(Object$r_0_min1)) <= 1){
      b <- cbind(NA, 0, -1)}
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      c <- cbind(median(Object$r_1_min1), 1, -1)}
    if (length(unique(Object$r_1_min1)) <= 1){
      c <- cbind(NA, 1, -1)}
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      d <- cbind(median(Object$r_min1_0), -1, 0)}
    if (length(unique(Object$r_min1_0)) <= 1){
      d <- cbind(NA, -1, 0)}
    
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      e <- cbind(median(Object$r_0_0), 0, 0)}
    if (length(unique(Object$r_0_0)) <= 1){
      e <- cbind(NA, 0, 0)}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      f <- cbind(median(Object$r_1_0), 1, 0)}
    if (length(unique(Object$r_1_0)) <= 1){
      f <- cbind(NA, 1, 0)}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      g <- cbind(median(Object$r_min1_1), -1, 1)}
    if (length(unique(Object$r_min1_1)) <= 1){
      g <- cbind(NA, -1, 1)}
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      h <- cbind(median(Object$r_0_1), 0, 1)}
    if (length(unique(Object$r_0_1)) <= 1){
      h <- cbind(NA, 0, 1)}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      i <- cbind(median(Object$r_1_1), 1, 1)}
    if (length(unique(Object$r_1_1)) <= 1){
      i <- cbind(NA, 1, 1)}
    
    data <- data.frame(rbind(a, b, c, d, e, f, g, h, i))
    names(data) <- c("Y", "Delta_T", "Delta_S")
    
    temp <- lattice::cloud(Y ~ as.factor(Delta_S) + as.factor(Delta_T), data, panel.3d.cloud=latticeExtra::panel.3dbars, 
                                xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1),  
                                par.settings = list(axis.line = list(col = "transparent")), 
                                xlab=expression(paste(Delta, "S")), ylab=expression(paste(Delta, T)), 
                                zlab="Median", col.facet=rep(c(2:4), each=3))
    
    plot(temp)
  }

  if (Type=="3D.Mode"){
    
    mode <- function(data) {
      x <- data
      if (unique(x[1])!=0){
        z <- density(x)
        mode_val <- z$x[which.max(z$y)]
        if (mode_val < 0){mode_val <- c(0)}
      }
      if (unique(x[1])==0){
        model_val <- c(0)
      }  
      fit <- list(mode_val= mode_val)  
    }
    
    
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      a <- cbind(mode(Object$r_min1_min1)$mode_val, -1, -1)}
    if (length(unique(Object$r_min1_min1)) <= 1){
      a <- cbind(NA, -1, -1)}
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      b <- cbind(mode(Object$r_0_min1)$mode_val, 0, -1)}
    if (length(unique(Object$r_0_min1)) <= 1){
      b <- cbind(NA, 0, -1)}
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      c <- cbind(mode(Object$r_1_min1)$mode_val, 1, -1)}
    if (length(unique(Object$r_1_min1)) <= 1){
      c <- cbind(NA, 1, -1)}
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      d <- cbind(mode(Object$r_min1_0)$mode_val, -1, 0)}
    if (length(unique(Object$r_min1_0)) <= 1){
      d <- cbind(NA, -1, 0)}
    
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      e <- cbind(mode(Object$r_0_0)$mode_val, 0, 0)}
    if (length(unique(Object$r_0_0)) <= 1){
      e <- cbind(NA, 0, 0)}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      f <- cbind(mode(Object$r_1_0)$mode_val, 1, 0)}
    if (length(unique(Object$r_1_0)) <= 1){
      f <- cbind(NA, 1, 0)}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      g <- cbind(mode(Object$r_min1_1)$mode_val, -1, 1)}
    if (length(unique(Object$r_min1_1)) <= 1){
      g <- cbind(NA, -1, 1)}
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      h <- cbind(mode(Object$r_0_1)$mode_val, 0, 1)}
    if (length(unique(Object$r_0_1)) <= 1){
      h <- cbind(NA, 0, 1)}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      i <- cbind(mode(Object$r_1_1)$mode_val, 1, 1)}
    if (length(unique(Object$r_1_1)) <= 1){
      i <- cbind(NA, 1, 1)}
    
    data <- data.frame(rbind(a, b, c, d, e, f, g, h, i))
    names(data) <- c("Y", "Delta_T", "Delta_S")
    
    temp <- lattice::cloud(Y ~ as.factor(Delta_S) + as.factor(Delta_T), data, panel.3d.cloud=latticeExtra::panel.3dbars, 
                                xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1),  
                                par.settings = list(axis.line = list(col = "transparent")), 
                                xlab=expression(paste(Delta, "S")), ylab=expression(paste(Delta, T)), 
                                zlab="Mode", col.facet=rep(c(2:4), each=3))
    plot(temp)  
  }
  
  if (Type=="3D.Spinning.Mean"){
    
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      a <- cbind(mean(Object$r_min1_min1), -1, -1)}
    if (length(unique(Object$r_min1_min1)) <= 1){
      a <- cbind(NA, -1, -1)}
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      b <- cbind(mean(Object$r_0_min1), 0, -1)}
    if (length(unique(Object$r_0_min1)) <= 1){
      b <- cbind(NA, 0, -1)}
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      c <- cbind(mean(Object$r_1_min1), 1, -1)}
    if (length(unique(Object$r_1_min1)) <= 1){
      c <- cbind(NA, 1, -1)}
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      d <- cbind(mean(Object$r_min1_0), -1, 0)}
    if (length(unique(Object$r_min1_0)) <= 1){
      d <- cbind(NA, -1, 0)}
    
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      e <- cbind(mean(Object$r_0_0), 0, 0)}
    if (length(unique(Object$r_0_0)) <= 1){
      e <- cbind(NA, 0, 0)}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      f <- cbind(mean(Object$r_1_0), 1, 0)}
    if (length(unique(Object$r_1_0)) <= 1){
      f <- cbind(NA, 1, 0)}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      g <- cbind(mean(Object$r_min1_1), -1, 1)}
    if (length(unique(Object$r_min1_1)) <= 1){
      g <- cbind(NA, -1, 1)}
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      h <- cbind(mean(Object$r_0_1), 0, 1)}
    if (length(unique(Object$r_0_1)) <= 1){
      h <- cbind(NA, 0, 1)}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      i <- cbind(mean(Object$r_1_1), 1, 1)}
    if (length(unique(Object$r_1_1)) <= 1){
      i <- cbind(NA, 1, 1)}
    
    data <- data.frame(rbind(a, b, c, d, e, f, g, h, i))
    names(data) <- c("Y", "Delta_T", "Delta_S")
    
    x <- as.factor(data$Delta_S)
    y <- as.factor(data$Delta_T)
    z <- data$Y
    
    rgl::plot3d(x=x, y=y, z=z, col="red", size=3, type="s",
           xlab="Delta_S", ylab="Delta_T", 
           zlab="Mean")  
    
  }
  
  
  if (Type=="3D.Spinning.Median"){
    
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      a <- cbind(median(Object$r_min1_min1), -1, -1)}
    if (length(unique(Object$r_min1_min1)) <= 1){
      a <- cbind(NA, -1, -1)}
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      b <- cbind(median(Object$r_0_min1), 0, -1)}
    if (length(unique(Object$r_0_min1)) <= 1){
      b <- cbind(NA, 0, -1)}
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      c <- cbind(median(Object$r_1_min1), 1, -1)}
    if (length(unique(Object$r_1_min1)) <= 1){
      c <- cbind(NA, 1, -1)}
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      d <- cbind(median(Object$r_min1_0), -1, 0)}
    if (length(unique(Object$r_min1_0)) <= 1){
      d <- cbind(NA, -1, 0)}
    
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      e <- cbind(median(Object$r_0_0), 0, 0)}
    if (length(unique(Object$r_0_0)) <= 1){
      e <- cbind(NA, 0, 0)}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      f <- cbind(median(Object$r_1_0), 1, 0)}
    if (length(unique(Object$r_1_0)) <= 1){
      f <- cbind(NA, 1, 0)}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      g <- cbind(median(Object$r_min1_1), -1, 1)}
    if (length(unique(Object$r_min1_1)) <= 1){
      g <- cbind(NA, -1, 1)}
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      h <- cbind(median(Object$r_0_1), 0, 1)}
    if (length(unique(Object$r_0_1)) <= 1){
      h <- cbind(NA, 0, 1)}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      i <- cbind(median(Object$r_1_1), 1, 1)}
    if (length(unique(Object$r_1_1)) <= 1){
      i <- cbind(NA, 1, 1)}
    
    data <- data.frame(rbind(a, b, c, d, e, f, g, h, i))
    names(data) <- c("Y", "Delta_T", "Delta_S")
    
    x <- as.factor(data$Delta_S)
    y <- as.factor(data$Delta_T)
    z <- data$Y
    
    plot3d(x=x, y=y, z=z, col="red", size=3, type="s",
           xlab="Delta_S", ylab="Delta_T", 
           zlab="Median")  
    
  }
  
  if (Type=="3D.Spinning.Mode"){
    
    mode <- function(data) {
      x <- data
      if (unique(x[1])!=0){
        z <- density(x)
        mode_val <- z$x[which.max(z$y)]
        if (mode_val < 0){mode_val <- c(0)}
      }
      if (unique(x[1])==0){
        model_val <- c(0)
      }  
      fit <- list(mode_val= mode_val)  
    }
    
    
    # T = -1, S = -1
    if (length(unique(Object$r_min1_min1)) > 1){
      a <- cbind(mode(Object$r_min1_min1)$mode_val, -1, -1)}
    if (length(unique(Object$r_min1_min1)) <= 1){
      a <- cbind(NA, -1, -1)}
    
    # T = 0, S = -1
    if (length(unique(Object$r_0_min1)) > 1){
      b <- cbind(mode(Object$r_0_min1)$mode_val, 0, -1)}
    if (length(unique(Object$r_0_min1)) <= 1){
      b <- cbind(NA, 0, -1)}
    
    # T = 1, S = -1
    if (length(unique(Object$r_1_min1)) > 1){
      c <- cbind(mode(Object$r_1_min1)$mode_val, 1, -1)}
    if (length(unique(Object$r_1_min1)) <= 1){
      c <- cbind(NA, 1, -1)}
    
    # T = -1, S = 0
    if (length(unique(Object$r_min1_0)) > 1){
      d <- cbind(mode(Object$r_min1_0)$mode_val, -1, 0)}
    if (length(unique(Object$r_min1_0)) <= 1){
      d <- cbind(NA, -1, 0)}
    
    # T = 0, S = 0
    if (length(unique(Object$r_0_0)) > 1){
      e <- cbind(mode(Object$r_0_0)$mode_val, 0, 0)}
    if (length(unique(Object$r_0_0)) <= 1){
      e <- cbind(NA, 0, 0)}
    
    # T = 1, S = 0
    if (length(unique(Object$r_1_0)) > 1){
      f <- cbind(mode(Object$r_1_0)$mode_val, 1, 0)}
    if (length(unique(Object$r_1_0)) <= 1){
      f <- cbind(NA, 1, 0)}
    
    # T = -1, S = 1
    if (length(unique(Object$r_min1_1)) > 1){
      g <- cbind(mode(Object$r_min1_1)$mode_val, -1, 1)}
    if (length(unique(Object$r_min1_1)) <= 1){
      g <- cbind(NA, -1, 1)}
    
    # T = 0, S = 1
    if (length(unique(Object$r_0_1)) > 1){
      h <- cbind(mode(Object$r_0_1)$mode_val, 0, 1)}
    if (length(unique(Object$r_0_1)) <= 1){
      h <- cbind(NA, 0, 1)}
    
    # T = 1, S = 1
    if (length(unique(Object$r_1_1)) > 1){
      i <- cbind(mode(Object$r_1_1)$mode_val, 1, 1)}
    if (length(unique(Object$r_1_1)) <= 1){
      i <- cbind(NA, 1, 1)}
    
    data <- data.frame(rbind(a, b, c, d, e, f, g, h, i))
    names(data) <- c("Y", "Delta_T", "Delta_S")
    
    x <- as.factor(data$Delta_S)
    y <- as.factor(data$Delta_T)
    z <- data$Y
    
    rgl::plot3d(x=x, y=y, z=z, col="red", size=3, type="s",
           xlab="Delta_S", ylab="Delta_T", 
           zlab="Mode")  
  }
  
  
  
}


    