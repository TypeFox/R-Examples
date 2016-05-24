#   Mega2: Manipulation Environment for Genetic Analysis
#   Copyright (C) 1999-2013 Robert Baron, Charles P. Kollar,
#   Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
#   Daniel E. Weeks, and University of Pittsburgh
#  
#   This file is part of the Mega2 program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 3 of the License, or (at your option) any later
#   version.
#  
#   Mega2 is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#  
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
# 
# ===========================================================================
#NPL plot R package

## ltypes  = should be 1 to the number of statistics

## bw=TRUE

## bw=TRUE, plot.width=par()$pin[1]

## ymax=3
## ymin=0
## yfix=F

## draw.lgnd=T

nplplot<-function(plotdata=NULL, filename=NULL, yline=2.0, ymin=0, ymax=3.0,
                  header=TRUE, yfix=FALSE, title=NULL, draw.lgnd=TRUE,
                  xlabl="", ylabl="", lgndx=NULL, lgndy=NULL, lgndtxt=NULL,
                  cex.legend = 0.7, cex.axis=0.7, tcl=1,
                  bw=TRUE, my.colors=NULL, ltypes=NULL, ptypes=NULL,
                  na.rm=TRUE, plot.width=0, ...)
{

  leg = NULL; ltys = NULL; pchs = NULL
  
  # use a plot command instead of abline to set plot parameters
  # If the user has supplied titles, then use them.
  # There should be as many as there are files (number of plots).
  # If there are fewer titles, recycle the last one.
    
  if(!is.null(title)) {
    maint <- title
  }
  
  # Figure out how the data has been supplied

  if (is.null(plotdata) && is.null(filename)) {
    print("Both plotdata and file are NULL, no data to plot!")
    return(F);
  }

  if (!is.null(plotdata) && !is.null(filename)) {
    print("Both plotdata and file are specified, using plotdata.")
  }
  
  if (!is.null(plotdata)) {
    lods <- plotdata
    print("lods read in as table")
    
    if(ncol(lods) <= 2) {
      print(paste("plotdata does not have any data."));
      return(F)
    }
  }
  else {
    if (!is.null(filename)) {      
      lods <- read.table(filename, header=header, sep="", na.strings=c("NA", "."))
     
      if (ncol(lods) <= 2) {
        print(paste("File ", filename, "does not have any data."));
        return(F)
      }
    }
  }

# attach(lods)
  
  rows <- dim(lods)[1]
  cols <- dim(lods)[2]  
  dist <- lods[, 2]
  
  ydatamax <- max(lods[, 3:cols], na.rm=TRUE);
  ydatamin <- min(lods[, 3:cols], na.rm=TRUE);
  xdatamin <- min(dist, na.rm=TRUE);
  xdatamax <- max(dist, na.rm=TRUE);
  
  if(is.null(ymax)) {
    ymax1 <- ydatamax + 0.2;
  }
  else {
    if (yfix == TRUE) {
      ymax1 <- ymax
    }
    else {
      ymax1 <- max(ydatamax, yline, ymax);
    }
  }
  
  if(is.null(ymin)) {
    ymin1 <- ydatamin - 0.1;
  }
  else {
    if (yfix == TRUE) {
      ymin1 <- ymin;
    }
    else {
      ymin1 <- min(ymin, ydatamin);
    }
  }

  if ((yfix == TRUE) && ((yline > ymax1) || (yline < ymin1))) {
    if (yline > ymax1) {
      print(paste("Y-line exceeds Y-max fixed at ", as.character(ymax1)))
    }
    else if (yline < ymin1) {
      print(paste("Y-line falls below Y-min fixed at ", as.character(ymin1)))
    }
    print("Y-line will not be drawn.")
    yline <- NA
  }

  if (is.null(ltypes) && is.null(ptypes)) {
    ltypes = 1:(cols-2)
  }

  plot(c(xdatamin, xdatamax), c(ymin1, ymax1), type="n",
       lty=1,
       xlab=xlabl, ylab=ylabl, xlim=c(xdatamin, xdatamax), 
       ylim=c(ymin1, ymax1), tcl=tcl, ...)
  

  if(!is.na(yline)) {
    abline(h=yline, col="grey40", ...)
  }
  
  if(!is.null(title)) {
    title(sub=maint, line=2);
  }
  
  if (bw==FALSE) {
    if (is.null(my.colors)) {
      my.colors <- rainbow(cols-2)
    }
  }
  
  for (k in 3:cols) {
#   scores<-get(names(lods)[k])
    scores<-lods[,k]
    ifelse(!is.null(ltypes), lt<- ltypes[k-2], lt<-"none")
    ifelse(!is.null(ptypes), ch<- ptypes[k-2], ch<-"none")
    
    this.x<-dist
    this.y<-scores
    
    if (bw == FALSE) {
      color <- my.colors[k-2]
    }
    else {
      color <- "black"
    }
    
    
    if (na.rm) {
      this.x <- this.x[!is.na(this.y)];
      this.y <- this.y[!is.na(this.y)];
    }

    if (lt == 0) { lt = "none"}
    if (ch == 0) { ch = "none"}

    if(lt == "none" && ch != "none") {
      points(this.x, this.y, pch = ch, col=color, ...);
      pchs[k-2] <- ch;
      ltys[k-2] <- "blank";
    }
    if(lt != "none" && ch == "none") {
      lines(this.x, this.y, type = "l", lty = lt,
            xaxt="n", yaxt="n", col=color, ...);
      pchs[k-2] = " ";
      ltys[k-2] = lt;
    }
    if(lt != "none" && ch != "none") {
      lines(this.x, this.y, type="b", lty = lt, pch=ch,
            xaxt="n", yaxt="n", col=color, ...);
      pchs[k-2] <- ch;
      ltys[k-2] <- lt;
    }
    
    if(ch == "none" && lt == "none") {
      print("Skipping this line");
    }
  }
  
  markers<-as.character(lods[,1]);
  # check if too markers are very close together
  lbl<-markers[markers != "-"];
  lloc<-dist[markers != "-"];

  # If a plot width has been specified, check if labels are
  # too close

  if (plot.width == 0.0) {
    plot.width <- (par()$pin)[1]
  }
  last.pos<-lloc[1]/(lloc[length(lloc)]-lloc[1]) * plot.width;
  for (j in 2:length(lloc)) {
    diff<-lloc[j]/(lloc[length(lloc)]-lloc[1]) * plot.width - last.pos;
    if (diff < 0.15) {
      lbl[j]<- "    ";
    }
    else {
      last.pos<-last.pos+diff;
    }
  }
  
  axis(3, tck=0.05, at=lloc,labels=lbl,cex.axis=cex.axis,las=2);
  
  
  # If the user has supplied legends, then use these
  # There should be as many legends as the number of data columns.
  # If there are fewer, use the names of the table

  
  if (draw.lgnd) {
    if (is.null(lgndtxt) || length(lgndtxt) < (cols - 2)) {
      leg <- names(lods)[3:cols]
    }
    else {
      leg <- lgndtxt
    }
    
    if (! is.null(lgndx)) {
      lgx <- lgndx
    }
    else {
      lgx <- xdatamin + 0.05*(xdatamax-xdatamin)
    }
    if (! is.null(lgndy)) {
      lgy <- lgndy
    }
    else {
      lgy <- ymin1 + 0.9*(ymax1 - ymin1)
    }
    
    if (!is.null(my.colors)) {
      leg.color<- my.colors;
    }
    else {
      leg.color <- "black";
    }

    legend(lgx, lgy, leg, col=leg.color, lty=ltys, pch=pchs, cex=cex.legend);
  }
  
# detach(lods);
  return(TRUE)
  
}
