
#Function to read the data------
ReadData<-function(regions_file, Sample_list,copynumber450k=FALSE){

  #Functions for checking the headers:::
  #this function to check the header of the input file: ("Sample" "Chromosome" "bp.Start"  "bp.End" "Num.of.Markers" "Mean")
  checkHeaderRegions <- function(file){  
    colnamesFile<-c("Sample","Chromosome","bp.Start","bp.End","Num.of.Markers","Mean")
    if(sum(colnames(file)==colnamesFile)==6){
      print("The header of the regions is ... OK")
    }
    else{
      print("The header of the file is not OK...should be: Sample/Chromosome/bp.Start/bp.End/Num.of.Markers/Mean")
    }
    invisible(sum(colnames(file)==colnamesFile))  #if 6 then it is OK ekse there is problem in the header
  }
  
  #this function to check the Sample info
  #The file of the samples info shoud have the header ("Sample","Comment")
  checkHeaderSL <- function(file){  
    colnamesSL<-c("Number","Sample","Comment")
    if(sum(colnames(file)==colnamesSL)==3){
      print("The header of the samples list is ... OK")
    }
    else{
      print("The header of the sample file is not OK...should be: Number / Sample / Comment")
    }
    invisible(sum(colnames(file)==colnamesSL))  #if 3 then it is OK if else there is problem in the header
  }

  if(is.character(regions_file)){
    regions<-read.csv(regions_file,stringsAsFactors =FALSE)   #load the file that contains regions and means
  } else {
    regions <- regions_file
  }
  if (copynumber450k==TRUE){
    regions<-regions[,c("Sample","chrom","loc.start","loc.end","num.mark","seg.mean")]
    colnames(regions)<-c("Sample","Chromosome","bp.Start","bp.End","Num.of.Markers","Mean")
  }
  
  
  if(is.character(regions_file)){
    SL<-read.csv(Sample_list,stringsAsFactors =FALSE)         #load the file that contains the names of samples and the comments
  } else {
    SL <- Sample_list
  }
  
  checkHeaderRegions(regions)
  checkHeaderSL(SL)
  
  #regions[is.na(regions)] <- 0

  object <- list(
    mainDir = getwd(),
    regions = regions,
    regions_save = regions,
    regions_auto = regions,
    SL = SL)
  class(object) <- "CopyNumber450kCancer_data"

  mod<-as.data.frame(SL$Sample,stringsAsFactors =FALSE)  # copy to store the modification in it
  mod[,2:6]<-0
  mod[,6]<-"No"
  mod[is.na(mod)] <- 0
  colnames(mod)<-c("Sample","Lower_selected_level","Upper_selected_level","Mean_of_selected_regions","Shifting","Reviewed?")
  object$mod <- mod
  
  mod_auto<-mod[,c("Sample","Mean_of_selected_regions","Shifting","Reviewed?")]
  colnames(mod_auto)<-c("Sample","Auto_Maximum_Peak","Shifting","Auto_Corrected?")
  object$mod_auto <- mod_auto
  
  object
}

#---------------------------------------------------------------
print.CopyNumber450kCancer_data <- function(x, ...){
  cat(sprintf("CopyNumber450kCancer data with %i samples in regions and %i samples in sample list.\n",
              length(unique(x$regions$Sample)),
              length(unique(x$SL$Sample))
  ))
  cat("Contains:\n", paste("  ", names(x), collapse = "\n"), "\n", sep="")
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#functions for plotting:::
#---to plot-------------

#this uses the regions file: Chromosome should be in this format: "chr1"
#similar to the original function, this one use only the cutoff
plotRegions<-function(object, chr, start, end, cutoff=0.1,markers=20, ...) {   
  sample_segments <- object
  
  if(hasArg(markers)){
    sample_segments$Mean[which(sample_segments$Num.of.Markers<=markers)]<-0
  }
  
  segment_values <- as.numeric(sample_segments[,"Mean"])
  segment_colors <- rep("black", nrow(sample_segments))
  
  if (missing(cutoff)) {
    cutoff<-(0.1)
  }
  
  segment_colors[as.numeric(segment_values) >= cutoff] <- "green"
  segment_colors[as.numeric(segment_values) <= -cutoff] <- "red"
  
  if (missing(chr)) {
    # Plotting the whole genome
    chromosomes <- unique(sample_segments[, "Chromosome"])
    site_per_chr <- cumsum(c(0, sapply(chromosomes, function(chr) max(as.numeric(sample_segments[sample_segments[,"Chromosome"] == chr, "bp.End"])))))
    offset <- site_per_chr - min(as.numeric(sample_segments[sample_segments[, "Chromosome"] == "chr1", "bp.Start"])) # 1 instead of "chr1" #as.numeric(gsub("\\D", "", x))
    start <- 0
    end <- as.numeric(max(site_per_chr))
    x_axis_type <- "n"
  } else {
    # Plotting a region
    if (missing(start)) {
      start <- 0
    }
    
    if (missing(end)) {
      end <- as.numeric(max(sample_segments[sample_segments[, "Chromosome"] == chr, "bp.End"]))
    }
    
    chromosomes <- chr
    offset <- 0
    x_axis_type <- NULL
  }
  
  yMin <- (-1) #min(c(-1, as.numeric(sample_segments[significant_segments, "Mean"])))
  yMax <- 1 #max(c(1, as.numeric(sample_segments[significant_segments, "Mean"])))
  
  #if (missing(ylab)) {
  # ylab<-""
  #}
  
  myPlot <- plot(range(start, end), range(yMin, yMax), type = "n",axes=FALSE, xaxt = x_axis_type, xlab="",ylab="", ...) #ylab="Mean",
  
  #---this function to plot the tick on the right side
  tick.tick<-function (nx = 2, ny = 2, tick.ratio = 0.5) {
    ax <- function(w, n, tick.ratio) {
      range <- par("usr")[if (w == "x") 
        1:2
        else 3:4]
      tick.pos <- if (w == "x") 
        par("xaxp")
      else par("yaxp")
      distance.between.minor <- (tick.pos[2] - tick.pos[1])/tick.pos[3]/n
      possible.minors <- tick.pos[1] - (0:100) * distance.between.minor
      low.minor <- min(possible.minors[possible.minors >= range[1]])
      if (is.na(low.minor)) 
        low.minor <- tick.pos[1]
      possible.minors <- tick.pos[2] + (0:100) * distance.between.minor
      hi.minor <- max(possible.minors[possible.minors <= range[2]])
      if (is.na(hi.minor)) 
        hi.minor <- tick.pos[2]
      axis(if (w == "x") 
        1
        else 4, seq(low.minor, hi.minor, by = distance.between.minor), 
        labels = FALSE, tcl = par("tcl") * tick.ratio)
    }
    if (nx > 1) 
      ax("x", nx, tick.ratio = tick.ratio)
    if (ny > 1) 
      ax("y", ny, tick.ratio = tick.ratio)
    invisible()
  }
  
  if (missing(chr)) {
    xlabs <- sapply(2:length(site_per_chr), function(j) {
      ((site_per_chr[j] - site_per_chr[j - 1])/2) + site_per_chr[j - 1]
    })
    
    axis(1, at = xlabs, labels = chromosomes, lty = 0, las = 2, ...)
    axis(4)
    tick.tick(nx=0,ny=2, tick.ratio=1.6)
    tick.tick(nx=0,ny=10, tick.ratio=0.6)
    mtext("L-value", side = 4, line = 2, cex = par("cex.lab"))
    box()
    abline(v = site_per_chr, lty = 3)
    abline(h = c(0,-cutoff,cutoff), lty = 3)
  }
  
  lapply(1:length(chromosomes), function(i) {
    used_segments <- sample_segments[, "Chromosome"] == chromosomes[i]
    colors <- segment_colors[used_segments]
    starts <- as.numeric(sample_segments[used_segments, "bp.Start"]) + offset[i]
    ends <- as.numeric(sample_segments[used_segments, "bp.End"]) + offset[i]
    y <- as.numeric(sample_segments[used_segments, "Mean"])
    graphics::segments(starts, y, ends, y, col = colors, lwd = 2, lty = 1)
  })
  
  myPlot
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#The Auto correction based on the highest peak-------------------
#this function correct the mean of the regions based on the highest peak and plot them
AutoCorrectPeak<-function(object, cutoff=0.1,markers=20, ...){
  subDir<-"Auto_corrected_plots"
  dir.create(file.path(object$mainDir, subDir), showWarnings = FALSE)  #mainDir was saved in the readData function
  setwd(file.path(object$mainDir, subDir)) 
  
  if (missing(markers)) {markers<-(0)}
  
  #prepare the QC file
  QC<-object$SL[,1:2]
  QC[,3:7]<-0
  colnames(QC) <- c("Sample","Comment","peak.sharpness","number.of.regions","IQR","SD","MAPD")
  
  
  par(mfrow=c(2,2),mar= c(5.1,0,4.1,0),oma=c(2,0,0,4))
  layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), widths=c(3,21), heights=c(10,10), TRUE) 
  
  for (i in 1:length(object$SL[,"Sample"])){ #correction
    print(paste("Auto Correction....Sample number",i))
    
    sam <- object$regions_auto[which(object$regions_auto$Sample %in% as.character(object$SL[i,"Sample"])),]
    forDen<-sam[which(sam$Chromosome!="chrX" & sam$Chromosome!="chrY"),c("Num.of.Markers","Mean")]
    
    sam.original<-sam
    
    d<-density(forDen$Mean,weights=forDen$Num.of.Markers/sum(forDen$Num.of.Markers),na.rm=TRUE,kernel = c("gaussian"),adjust=0.15,n=512)
    max.peak.value<-d$x[which.max(d$y)]
    sam$Mean <- sam$Mean-max.peak.value
    
    #QC (peak sharpness) + (count the regions)
    point<-which(d$x==(max.peak.value))  #for peak sharpness
    QC.peak.sharpness <- ((d$y[point+20]+d$y[point-20])/2)/d$y[which(d$x==max.peak.value)]
    QC[i,"peak.sharpness"] <- as.numeric(QC.peak.sharpness)
    QC[i,"number.of.regions"]  <- length(sam[,1]) # for number of the regions
    
    QC[i,"IQR"]<-IQR(forDen[,"Mean"], na.rm = TRUE, type = 7) #for IQR
    QC[i,"SD"]<-sd(forDen[,"Mean"], na.rm = TRUE) #for SD
    QC[i,"MAPD"]<- median(abs(diff(forDen[,"Mean"],na.rm = TRUE)), na.rm = TRUE)# for Median Absolute Pairwise Difference (MAPD)
    
    object$regions_auto[which(object$regions_auto$Sample %in% as.character(object$SL[i,"Sample"])),"Mean"]<-sam$Mean #store the modifications
    object$mod_auto[which(object$mod_auto$Sample %in% as.character(object$SL[i,"Sample"])),2:4]<-c(round(max.peak.value, 3),round(-max.peak.value, 3),"Auto") #store the modifications
    
    print(paste("Plotting....Sample number",i)) #plotting
    
    png(filename = paste(i,"_",object$SL[i,"Sample"],"_auto_corrected_plot.png",sep=""), width = 1920, height = 1200)
    
    par(mfrow=c(2,2),mar= c(5.1,0,4.1,0),oma=c(2,0,0,4))
    layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), widths=c(3,21), heights=c(10,10), TRUE) 
    
    #original plots
    plot(d$y,d$x,ylim=c(-1,1),type='l',ylab="",xlab="",axes=FALSE,xlim=rev(range(d$y)))
    abline(h = c(0,-cutoff,cutoff), lty = 3)
    box()
    legend("bottomleft", legend="Density",cex=1)
    plotRegions(sam.original, cutoff=cutoff, markers=markers,
        main=paste("Sample::",object$SL[i,"Sample"],"       Info::",object$SL[i,"Comment"]), ...)
    
    
    #new plots
    plot(d$y,d$x-max.peak.value,ylim=c(-1,1),type='l',ylab="",xlab="",axes=FALSE,xlim=rev(range(d$y)))
    abline(h = c(0,-cutoff,cutoff), lty = 3)
    box()
    legend("bottomleft", legend="Density",cex=1)
    plotRegions(sam,cutoff=cutoff,markers=markers,main=paste("Autocorrected plot"))
    
    #mtext("L-value", side = 2, line = 2, cex = par("cex.lab"))
    dev.off()
  }
  
  setwd(file.path(object$mainDir))
  object$QC <- QC
  
  #added for the new version
  object$regions<-object$regions_auto
  
  print("Saving the file of the autocorrection files...")
  write.csv(object$regions_auto,file="autocorrected_regions.csv")
  write.csv(object$mod_auto,file="autocorrections.csv")
  write.csv(QC,file="QC.csv",row.names = FALSE)
  
  print("Auto Correction Done.")
  object
}


#AutoCorrectPeak()
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#functions for modifing and reviewing the plots:::

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#saving the reviewing results and the modifications and the plots #this only after the reviewing
# to run it saveOutput() or saveOutput(plots) to save with plots
# ReviewPlot()
ReviewPlot<-function(object, select,plots=TRUE,cutoff=0.1,markers=20,...){
 
  setwd(file.path(object$mainDir)) 
  
  if (missing(select)) {
    select<-c(1:length(object$SL[,"Sample"]))
  }
  
  to.stop <- 0
  
  # this function is to calculate the weighted median 
  weighted.median <- function (pp, w) {
    ox <- order(pp)
    pp <- pp[ox]
    w <- w[ox]
    k <- 1
    low <- cumsum(c(0, w))
    up <- sum(w) - low
    df <- low - up
    repeat {
      if (df[k] < 0) 
        k <- k + 1
      else if (df[k] == 0) 
        return((w[k] * pp[k] + w[k - 1] * pp[k - 1])/(w[k] + w[k - 1]))
      else return(pp[k - 1])
    }
  }
  
  # this function for shifting and review-------------
  Review<-function(name,cutoff=0.1,markers=20,...){ # sample name  
    to.stop <<- 0
    
    sam <- object$regions[which(object$regions$Sample %in% as.character(name)),]   #get the sample segments
    
    par(mfrow=c(2,2),mar= c(5.1,0,4.1,0),oma=c(2,0,0,4))
    layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), widths=c(3,21), heights=c(10,10), TRUE) 
    
    forDen<-sam[which(sam$Chromosome!="chrX" & sam$Chromosome!="chrY"),c("Num.of.Markers","Mean")]
    d<-density(forDen$Mean,weights=forDen$Num.of.Markers/sum(forDen$Num.of.Markers),na.rm=TRUE,kernel = c("gaussian"),adjust=0.15,n=512)
    plot(d$y,d$x,ylim=c(-1,1),type='l',ylab="",xlab="",axes=FALSE,xlim=rev(range(d$y)))
    abline(h = c(0,-cutoff,cutoff), lty = 3)
    box()
    legend("bottomleft", legend="Density",cex=1)
    
    plotRegions(sam,cutoff=cutoff,main=c(paste("Sample::",object$SL[which(object$SL$Sample %in% as.character(name)),"Sample"]),
                                         "             Info::",object$SL[which(object$SL$Sample %in% as.character(name)),"Comment"]), ...)
    legend("topleft", legend=c("2 Clicks on the plot to determine the range","Back/Next to skip this sample","Stop to quit"),cex=0.75)  
    
    
    
    
    A <- 0 #for save #for change between meadian to peak. 
    B <- 500000000  #for save #for change between meadian to peak.
    A2<- 600000000 #for modify
    B2<- 1100000000 #for modify
    
    A5<- 1200000000 #for back
    B5<- 1700000000 #for back
    A3<- 1800000000 #for skip
    B3<- 2300000000 #for skip
    A4<- 2400000000 #for stop
    B4<- 2900000000 #for stop

        
    Y1<- (-1)
    Y2<- (-0.85)
    
    rect(A3, Y1, B3, Y2, density = NULL, angle = 45, col = "green", border = NULL, lty = par("lty"), lwd = par("lwd"))
    rect(A4, Y1, B4, Y2, density = NULL, angle = 45, col = "red", border = NULL, lty = par("lty"), lwd = par("lwd"))
    rect(A5, Y1, B5, Y2, density = NULL, angle = 45, col = "green", border = NULL, lty = par("lty"), lwd = par("lwd"))
    
    text(((A5 + B5)/2), ((Y1 + Y2)/2), labels = paste("Back"), cex = 0.75)  
    text(((A3 + B3)/2), ((Y1 + Y2)/2), labels = paste("Next"), cex = 0.75)  
    text(((A4 + B4)/2), ((Y1 + Y2)/2), labels = paste("Stop"), cex = 0.75) 
    
    #the 1st click
    click1 <- locator(n = 1, type = "n")
        
    if ((click1$x > A5 & click1$x < B5 & click1$y > (Y1) & click1$y < (Y2))) { #back
      print("Back to the previous sample")

      to.stop <<- 3
      return()#
    }
    
    if ((click1$x > A3 & click1$x < B3 & click1$y > (Y1) & click1$y < (Y2))) { #skip sample
      print("Skipped")

      to.stop <<- 2
      return()#
    }
    
    if ((click1$x > A4 & click1$x < B4 & click1$y > (Y1) & click1$y < (Y2))) { #stop reviewing
      print("Stopped")
      to.stop <<- 1
      graphics.off()
      return()#
    }
    
    #the 2nd click
    click2 <- locator(n = 1, type = "n")
    
    if ((click2$x > A5 & click2$x < B5 & click2$y > (Y1) & click2$y < (Y2))) { #back
      print("Back to the previous sample")

      to.stop <<- 3   
      return()#
    }
    
    if ((click2$x > A3 & click2$x < B3 & click2$y > (Y1) & click2$y < (Y2))) { #skip sample
      print("Skipped")

      to.stop <<- 2   
      return()#
    }
    
    if ((click2$x > A4 & click2$x < B4 & click2$y > (Y1) & click2$y < (Y2))) { #stop reviewing
      print("Stopped")
      to.stop <<- 1
      graphics.off()
      return()#
    }
    
    #Draw the line for the range
    up <- c(click1$y, click2$y)
    up <- up[order(as.numeric(up[]))]
    abline(h = c(up[1], up[2]), lty = 4)
    
    
    rect(A, Y1, B, Y2, density = NULL, angle = 45, col = "green", 
         border = NULL, lty = par("lty"), lwd = par("lwd"))
    rect(A2, Y1, B2, Y2, density = NULL, angle = 45, col = "green", 
         border = NULL, lty = par("lty"), lwd = par("lwd"))
    text(((A + B)/2), ((Y1 + Y2)/2), labels = paste("Median"), 
         cex = 0.75)
    text(((A2 + B2)/2), ((Y1 + Y2)/2), labels = paste("Peak"), 
         cex = 0.75)
    
    #the 3rd click
    click3 <- locator(n = 1, type = "n")
    
    while ((click3$x > A & click3$x < B & click3$y > (Y1) & click3$y < (Y2)) == FALSE & 
             (click3$x > A2 & click3$x < B2 & click3$y > (Y1) & click3$y < (Y2)) == FALSE & 
             (click3$x > A3 & click3$x < B3 & click3$y > (Y1) & click3$y < (Y2)) == FALSE & 
             (click3$x > A4 & click3$x < B4 & click3$y > (Y1) & click3$y < (Y2)) == FALSE) {
      click3 <- locator(n = 1, type = "n")
    }
    
    if ((click3$x > A5 & click3$x < B5 & click3$y > (Y1) & click3$y < (Y2))) { #back
      print("Back to the previous sample")

      to.stop <<- 3   
      return()#
    }
    
    if ((click3$x > A3 & click3$x < B3 & click3$y > (Y1) & click3$y < (Y2))) { #skip sample
      print("Skipped")

      to.stop <<- 2 
      return()#
    }
        
    if ((click3$x > A4 & click3$x < B4 & click3$y > (Y1) & click3$y < (Y2))) { #stop reviewing
      print("Stopped")
      to.stop <<- 1 
      graphics.off()
      return()#
    }
    
    
    if((click3$x>A&click3$x<B&click3$y>(Y1)&click3$y<(Y2))){ #median
      print("median")
      
      forMedian    <- sam[which(sam$Chromosome!="chrX" & sam$Chromosome!="chrY"),c("Num.of.Markers","Mean")]  #removing the X and Y chrs from the mean
      subsetMeans  <- forMedian$Mean[(which(forMedian$Mean>up[1] & forMedian$Mean<up[2]))]
      subsetMarkers<- forMedian$Num.of.Markers[(which(forMedian$Mean>up[1] & forMedian$Mean<up[2]))]
      move <- weighted.median(subsetMeans,subsetMarkers)
      sam$Mean <- sam$Mean-move
      
      #plot the modified plot
      forDen<-sam[which(sam$Chromosome!="chrX" & sam$Chromosome!="chrY"),c("Num.of.Markers","Mean")]
      d<-density(forDen$Mean,weights=forDen$Num.of.Markers/sum(forDen$Num.of.Markers),na.rm=TRUE,kernel = c("gaussian"),adjust=0.15,n=512)
      plot(d$y,d$x,ylim=c(-1,1),type='l',ylab="",xlab="",axes=FALSE,xlim=rev(range(d$y)))
      abline(h = c(0,-cutoff,cutoff), lty = 3)
      box()
      legend("bottomleft", legend="Density",cex=1)
      
      plotRegions(sam,cutoff=cutoff,...)
      legend("topleft", legend=c(
        paste("Upper:", round(up[2], 3) ), 
        paste("Lower:", round(up[1], 3) ),
        paste("W.Median:",round(move, 4) )),      
        title="Selected",cex=0.75)
      
      rect(A, Y1, B, Y2,density = NULL , angle = 45, col = "green", border = NULL, lty = par("lty"), lwd = par("lwd"))
      rect(A2, Y1, B2, Y2,density = NULL, angle = 45, col = "red", border = NULL, lty = par("lty"), lwd = par("lwd"))  
      text(((A+B)/2),((Y1+Y2)/2),labels=paste("Save"),cex=0.75 )
      text(((A2+B2)/2),((Y1+Y2)/2),labels=paste("modify"),cex=0.75)
      
      click<-locator(n = 1,type = "l")
      
      while ((click$x>A&click$x<B&click$y>(Y1)&click$y<(Y2))==FALSE & (click$x>A2&click$x<B2&click$y>(Y1)&click$y<(Y2))==FALSE){
        click<-locator(n = 1,type = "l")
      }   
      if((click$x>A&click$x<B&click$y>(Y1)&click$y<(Y2))){ #save
        print(paste(name,"modification...saved"))
        
        object$regions[which(object$regions$Sample %in% as.character(name)),"Mean"]<<-sam$Mean
        object$mod[which(object$mod$Sample %in% as.character(name)),2:6]<<-c(round(up[1], 3),round(up[2], 3),round(move, 4),round(-move, 4),"Yes_using_median")
        
        if(plots==TRUE){
          
          print(paste("Saving the plot for ...",name))
          
          subDir<-"reviewed_plots"
          dir.create(file.path(object$mainDir, subDir), showWarnings = FALSE)
          setwd(file.path(object$mainDir, subDir)) 
          
          png(filename = paste(object$SL[which(object$SL$Sample %in% as.character(name)),"Number"],"_",name,"_reviewed_plot.png",sep=""), width = 1920, height = 1200)
          par(mfrow=c(2,2),mar= c(5.1,0,4.1,0),oma=c(2,0,0,4))
          layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), widths=c(3,21), heights=c(10,10), TRUE) 
          
          #original plots
          sam.original <- object$regions_save[which(object$regions_save$Sample %in% as.character(name)),]
          forDen<-sam.original[which(sam.original$Chromosome!="chrX" & sam.original$Chromosome!="chrY"),c("Num.of.Markers","Mean")]
          d<-density(forDen$Mean,weights=forDen$Num.of.Markers/sum(forDen$Num.of.Markers),na.rm=TRUE,kernel = c("gaussian"),adjust=0.20,n=1024)
          
          plot(d$y,d$x,ylim=c(-1,1),type='l',ylab="",xlab="",axes=FALSE,xlim=rev(range(d$y)))
          abline(h = c(0,-cutoff,cutoff), lty = 3)
          box()
          legend("bottomleft", legend="Density",cex=1)
          plotRegions(sam.original,cutoff=cutoff,markers=markers,main=paste("Sample:: ",name,"       Info:: ",object$SL[which(object$SL$Sample %in% as.character(name)),"Comment"]))
          
          #the reviewed plots
          draw <- object$regions[which(object$regions$Sample %in% as.character(name)),] 
          forDen<-draw[which(draw$Chromosome!="chrX" & draw$Chromosome!="chrY"),c("Num.of.Markers","Mean")]
          d<-density(forDen$Mean,weights=forDen$Num.of.Markers/sum(forDen$Num.of.Markers),na.rm=TRUE,kernel = c("gaussian"),adjust=0.20,n=1024)
          
          plot(d$y,d$x,ylim=c(-1,1),type='l',ylab="",xlab="",axes=FALSE,xlim=rev(range(d$y)))
          abline(h = c(0,-cutoff,cutoff), lty = 3)
          box()
          legend("bottomleft", legend="Density",cex=1)
          plotRegions(draw,cutoff=cutoff,markers=markers,main=paste("Reviewed plot"))
          
          dev.off()
          setwd(file.path(object$mainDir))
          
        }
      }
      if((click$x>A2&click$x<B2&click$y>(Y1)&click$y<(Y2))){ #modify
        print("modify...")
        Review(name)
        return()
      }
      
    }
    
    
    
    if((click3$x>A2&click3$x<B2&click3$y>(Y1)&click3$y<(Y2))){ #peak
      print("peak")
      d.mod.x<-d$x[which(d$x>up[1] & d$x<up[2])]
      d.mod.y<-d$y[which(d$x>up[1] & d$x<up[2])]
      i<-which.max(d.mod.y)
      max.peak.value<-d.mod.x[i]
      sam$Mean <- sam$Mean-max.peak.value
      
      
      #plot the modified plot
      forDen<-sam[which(sam$Chromosome!="chrX" & sam$Chromosome!="chrY"),c("Num.of.Markers","Mean")]
      d<-density(forDen$Mean,weights=forDen$Num.of.Markers/sum(forDen$Num.of.Markers),na.rm=TRUE,kernel = c("gaussian"),adjust=0.15,n=512)
      plot(d$y,d$x,type='l',ylim=c(-1,1),ylab="",xlab="",axes=FALSE,xlim=rev(range(d$y)))
      abline(h = c(0,-cutoff,cutoff), lty = 3)
      box()
      legend("bottomleft", legend="Density",cex=1)
      
      plotRegions(sam,cutoff=cutoff,...)
      legend("topleft", legend=c(
        paste("Upper:", round(up[2], 3) ), 
        paste("Lower:", round(up[1], 3) ),
        paste("Peak at:",round(max.peak.value, 4) )),      
        title="Selected",cex=0.75)
      
      rect(A, Y1, B, Y2,density = NULL , angle = 45, col = "green", border = NULL, lty = par("lty"), lwd = par("lwd"))
      rect(A2, Y1, B2, Y2,density = NULL, angle = 45, col = "red", border = NULL, lty = par("lty"), lwd = par("lwd"))  
      text(((A+B)/2),((Y1+Y2)/2),labels=paste("Save"),cex=0.75 )
      text(((A2+B2)/2),((Y1+Y2)/2),labels=paste("modify"),cex=0.75)
      
      click<-locator(n = 1,type = "l")
      while ((click$x>A&click$x<B&click$y>(Y1)&click$y<(Y2))==FALSE & (click$x>A2&click$x<B2&click$y>(Y1)&click$y<(Y2))==FALSE){
        click<-locator(n = 1,type = "l")
      }   
      if((click$x>A&click$x<B&click$y>(Y1)&click$y<(Y2))){ #save
        print(paste(name,"modification...saved"))
        
        object$regions[which(object$regions$Sample %in% as.character(name)),"Mean"]<<-sam$Mean
        object$mod[which(object$mod$Sample %in% as.character(name)),2:6]<<-c(round(up[1], 3),round(up[2], 3),round(max.peak.value, 4),round(-max.peak.value, 4),"Yes_using_peak")
        
        if(plots==TRUE){
          
          print(paste("Saving the plot for ...",name))
          
          subDir<-"reviewed_plots"
          dir.create(file.path(object$mainDir, subDir), showWarnings = FALSE)
          setwd(file.path(object$mainDir, subDir)) 
          
          png(filename = paste(object$SL[which(object$SL$Sample %in% as.character(name)),"Number"],"_",name,"_reviewed_plot.png",sep=""), width = 1920, height = 1200)
          par(mfrow=c(2,2),mar= c(5.1,0,4.1,0),oma=c(2,0,0,4))
          layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), widths=c(3,21), heights=c(10,10), TRUE) 
          
          #original plots
          sam.original <- object$regions_save[which(object$regions_save$Sample %in% as.character(name)),]
          forDen<-sam.original[which(sam.original$Chromosome!="chrX" & sam.original$Chromosome!="chrY"),c("Num.of.Markers","Mean")]
          d<-density(forDen$Mean,weights=forDen$Num.of.Markers/sum(forDen$Num.of.Markers),na.rm=TRUE,kernel = c("gaussian"),adjust=0.20,n=1024)
          
          plot(d$y,d$x,ylim=c(-1,1),type='l',ylab="",xlab="",axes=FALSE,xlim=rev(range(d$y)))
          abline(h = c(0,-cutoff,cutoff), lty = 3)
          box()
          legend("bottomleft", legend="Density",cex=1)
          plotRegions(sam.original,cutoff=cutoff,markers=markers,main=paste("Sample:: ",name,"       Info:: ",object$SL[which(object$SL$Sample %in% as.character(name)),"Comment"]))
          
          #the reviewed plots
          draw <- object$regions[which(object$regions$Sample %in% as.character(name)),] 
          forDen<-draw[which(draw$Chromosome!="chrX" & draw$Chromosome!="chrY"),c("Num.of.Markers","Mean")]
          d<-density(forDen$Mean,weights=forDen$Num.of.Markers/sum(forDen$Num.of.Markers),na.rm=TRUE,kernel = c("gaussian"),adjust=0.20,n=1024)
          
          plot(d$y,d$x,ylim=c(-1,1),type='l',ylab="",xlab="",axes=FALSE,xlim=rev(range(d$y)))
          abline(h = c(0,-cutoff,cutoff), lty = 3)
          box()
          legend("bottomleft", legend="Density",cex=1)
          plotRegions(draw,cutoff=cutoff,markers=markers,main=paste("Reviewed plot"))
          
          dev.off()
          setwd(file.path(object$mainDir))
          
        }
        
      }
      if((click$x>A2&click$x<B2&click$y>(Y1)&click$y<(Y2))){ #modify
        print("modify...")
        Review(name)  
        return()
        
      } 
    }
    
  }
  
  i<-1
  while (i<=length(select)){ 
    if(i<1) {i<-1}
    
    print(paste("Sample number:",i, "      name:",object$SL[select[i],"Sample"]))
    Review(object$SL[select[i],"Sample"])
    
    if(to.stop==3){i<-(i-2)}
    if(to.stop==1){
      rm(to.stop)
      break
    }
    
    i<-i+1
  }
  
  setwd(file.path(object$mainDir)) 
  print("saving the file of the reviewed regions...please wait...")
  write.csv(object$regions,file="reviewed_regions.csv")
  print("saving the file of the modification...please wait...")
  write.csv(object$mod,file="Manual_corrections.csv")
  
  print("Done.")
  object
  
}

#-  ------------------------------------------------------
#----------------------------------------------------------------
#This function to plot only oone plot per page
PlotCNV<- function(object, select,comment=FALSE,cutoff=0.1,markers=20){
  
  if (missing(select)) {
    select<-c(1:length(object$SL[,"Sample"]))
  }
  
  print("Plotting...")
  
  subDir<-"plots"
  dir.create(file.path(object$mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(object$mainDir, subDir)) 
  
  for(i in 1:length(select)){
    #the plots
    png(filename = paste(select[i],"_",object$SL[select[i],"Sample"],"_plot.png",sep=""), width = 1920, height = 1200)
    
    par(mfrow=c(1,2),mar= c(5.1,0,4.1,0),oma=c(2,0,0,4))
    layout(matrix(c(1,2),1,2,byrow=TRUE), widths=c(3,21), heights=c(10,10), TRUE) 
    
    #the density
    draw <- object$regions[which(object$regions$Sample %in% object$SL[select[i],"Sample"]),] 
    forDen<-draw[which(draw$Chromosome!="chrX" & draw$Chromosome!="chrY"),c("Num.of.Markers","Mean")]
    d<-density(forDen$Mean,weights=forDen$Num.of.Markers/sum(forDen$Num.of.Markers),na.rm=TRUE,kernel = c("gaussian"),adjust=0.20,n=1024)
    
    plot(d$y,d$x,ylim=c(-1,1),type='l',ylab="",xlab="",axes=FALSE,xlim=rev(range(d$y)))
    abline(h = c(0,-cutoff,cutoff), lty = 3)
    box()
    legend("bottomleft", legend="Density",cex=1)
    
    info<-""
    if(comment==TRUE){
      Comment <- object$SL[select[i],"Comment"]
      info<-paste("       Info:: ",Comment)
    }
    
    plotRegions(draw,cutoff=cutoff,markers=markers,main=paste("Sample::",object$SL[select[i],"Sample"],info))
    
    dev.off()
  }
  
  graphics.off()
  setwd(file.path(object$mainDir))
}

#----------------------------------------------------------------
#for merging based on the cutoff and number of markers in the region
PlotMerged<-function(object, cutoff=0.1, markers=20, ...){
  
  file <- object$regions
  
  file<-file[which(file$Chromosome!="chrX" & file$Chromosome!="chrY"),] #removing X and Y chr.s
  
  if(hasArg(markers)){
    file$Mean[which(file$Num.of.Markers<=markers)]<-0
  }
  
  forChange<-file$Mean
  
  forChange[forChange>cutoff]<-(cutoff)+0.01
  forChange[forChange<(-cutoff)]<-(-cutoff)-0.01
  forChange[forChange<cutoff&forChange>(-cutoff)]<-0
  
  file$Mean<-forChange
  
  subDir<-"Merged_plots"
  dir.create(file.path(object$mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(object$mainDir, subDir)) 
  
  for (i in 1:length(object$SL[,"Sample"])){  #plotting
    print(paste("Plotting....Sample number",i))
    sam<-file[which(file$Sample %in% as.character(object$SL[i,"Sample"])),]
    png(filename = paste(i,"_",object$SL[i,"Sample"],"_merged_plot_","cutoff_",cutoff,".png",sep=""), width = 1920, height = 1200)
    plotRegions(sam,main=paste("CNV_450K::",object$SL[i,"Sample"],object$SL[i,"Comment"]),cutoff=cutoff)
    dev.off()
  }
  
  setwd(file.path(object$mainDir))
  print("Plotting... Done.")
  
}

#----------------------------------------------------------------
#-----------------------------END--------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
