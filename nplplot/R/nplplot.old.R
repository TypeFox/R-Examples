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
nplplot.old<-function(files, col=2, row=2, mode="p", output="screen", 
                  yline=2.0, ymin=NULL, ymax=NULL, yfix=FALSE, batch=FALSE,
                  headerfiles=NULL,  titles=NULL,
                  xlabl="", ylabl="", lgnd="page", lgndx=NULL, lgndy=NULL,
                  bw=TRUE, na.rm=TRUE)
{

#   my.colors<-c("magenta", "navyblue", "green", "lightblue", "grey", "pink",
#                "black");
#   if (length(files) < 3) {
#     row = 1;
#   }
#   if (length(files) < 2) {
#     col = 1;
#   }
  my.colors<-c("red", "navyblue", "green", "lightblue", "grey", "pink",
               "black");
  
  if(batch == FALSE) {
    old.par <- par(no.readonly = TRUE);
    if(output == "screen") {
      on.exit(par(old.par));
    }
    else {
      on.exit(par(old.par), add=TRUE);
      on.exit(dev.off(dev.list()));
    }
  }

  numfiles<-length(files);
  
  # page size in inches
  
  if(mode=="p") {
    page.width<-7.0*2.54;
    page.height<-9.0*2.54;
    full.width<-7.5;
    full.height<-10.5;
    landscape<-FALSE;
  }
  else {
    page.width<-9.0*2.54;
    page.height<-7.0*2.54;
    full.width<-10.5;
    full.height<-7.5;
    landscape<-TRUE;
  }

  plot.width<-array(page.width/col, col);
  plot.height<-array(page.height/row, row);
  
  plots.inpage<-col*row;
  file.names<-as.character(files);
  
  pl.wd <- (plot.width[1]-2.0)/2.54;
  pl.ht <- (plot.height[1]-2.0)/2.54;

  if(batch == FALSE) {
    dev.off(dev.list());
  }

  if(output == "screen") {
#    x11(width=full.width, height=full.height);
#     par(omi=c(0.05, 0.05, 0.05, 0.05),pin=c(pl.wd, pl.ht), 
    par(mfrow = c(row, col), ask=FALSE);
  }
  else {
    output1 <- output;
    fileparts <- unlist(strsplit(output1, "\\."));
    if (fileparts[length(fileparts)] == "ps") {
      print("Creating postscript file");
      postscript(width=full.width, height=full.height,
                 file=output, paper="letter", horizontal=landscape);
    }
    else {
      print("Creating pdf file");
#      print(c(full.width, full.height));
      pdf(width=full.width, height=full.height,
          file=output, paper="special");
    }
#     par(omi=c(0.1, 0.1, 0.5, 0.1),
#         mfrow = c(row, col), pin=c(pl.wd, pl.ht), ask=FALSE);
    par(mfrow = c(row, col), ask=FALSE);
#    mai=c(0.25, 0.25, 0.7, 0.1),pin=c(pl.wd, pl.ht),
  }

  # for each file
   for (i in 1:numfiles) {
    # for each curve
    leg = NULL; ltys = NULL; pchs = NULL;
    # use a plot command instead of abline to set plot parameters

    # If the user has supplied titles, then use them.
    # There should be as many as there are files (number of plots).
    # If there are fewer titles, recycle the last one.
    
    if(!is.null(titles)) {
      if(i <= length(titles)) {
        ititl<- i;
      }
      else {
        ititl<-length(titles);
      }
      maint <- as.character(titles[ititl]);
    }

    # If the user has supplied headerfiles, then read in a
    # line from each file as header
    # There should be as many header files as there are data files
    # (number of plots). If there are fewer header files, recycle
    # the last one.

    if(!is.null(headerfiles)) {
      if(i <= length(headerfiles)) {
        leg1<-scan(headerfiles[i], what="string");
      }
      else {
        leg1<-scan(headerfiles[length(headerfiles)], what="string");
      }
      if(length(leg1) == 0) {
        print(paste("Empty header file ", headerfiles[i]));
        return(F);
      }
      skip.header<-1;
      read.header<-FALSE;
      leg<-leg1[3:length(leg1)];
    }
    else {
      skip.header<-0;
      read.header<-TRUE;
    }

    tmp<-scan(file.names[i], what="string");
    if(length(tmp)==0) {
      print(paste("File ", file.names[i], "is empty!"));
      return(F);
    }
    lods<-read.table(file.names[i], header=read.header, skip=skip.header, sep="",
                     na.strings=c("NA", "."));
    if(ncol(lods) <= 2) {
      print(paste("File ", file.names[i], "does not have any data."));
      return(F);
    }
#   attach(lods);

    rows <- dim(lods)[1];
    cols <- dim(lods)[2];
    if(read.header == TRUE) {
      leg <- names(lods)[3:cols];
      # unlist(strsplit(as.character(names(lods)[k]), "\\.")); 
      # leg[k-2] <- paste(leg1[1], "-", leg1[2]);
    }
    dist <- lods[1:(rows-2), 2];
    ydatamax <- max(lods[1:(rows-2), 3:cols], na.rm=TRUE);
    ydatamin <- min(lods[1:(rows-2), 3:cols], na.rm=TRUE);

    xdatamin <- min(dist, na.rm=TRUE);
    #xdatamin <- max(0, xdatamin - 5);

    xdatamax <- max(dist, na.rm=TRUE);
    
    if(is.null(ymax)) {
      ymax1 <- ydatamax + 0.2;
    }
    else {
      if (yfix == TRUE) {
        ymax1 <- ymax;
      }
      else {
        ymax1 <- max(ydatamax, ymax);
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
    if(yline > ymax1) {
      yline <- 0.0;
    }

    if(yline == 0.0) {
      plot(c(xdatamin, xdatamax), c(yline, yline), type="l", col="white",
           lty=1,
           xlab=xlabl, ylab=ylabl, xlim=c(xdatamin, xdatamax),xpd=FALSE,
           ylim=c(ymin1, ymax1),
           tcl=-0.3, cex.axis=0.9);
    }
    else {
      plot(c(xdatamin, xdatamax), c(yline, yline), type="l", col="grey40",
           lty=1,
           xlab=xlabl, ylab=ylabl, xlim=c(xdatamin, xdatamax),xpd=FALSE,
           ylim=c(ymin1, ymax1),
           tcl=-0.3, cex.axis=0.9);
    }
      
    if(ymax1 < 0.0) {
      abline(h=0.0, col = "grey72");
    }
    if(!is.null(titles)) {
      title(sub=maint, line=2);
    }
    for (k in 3:cols) {
#     scores<-get(names(lods)[k]);
      scores<-lods[,k];
      ch <- scores[rows];
      lt <- scores[rows-1];
      if(bw == FALSE) {
        # cycle through the first 6 colors
        color = (k-3 %% 5) + 1;
      }
      else {
        # 7th color is black
        color=7;
      }
      this.x<-dist[1:(rows-2)];
      this.y<-scores[1:(rows-2)];

      if (na.rm) {
        this.x <- this.x[!is.na(this.y)];
        this.y <- this.y[!is.na(this.y)];
      }
      if(lt == 0 && ch > 0) {
        points(this.x, this.y, pch = ch, col=my.colors[color]);
        pchs[k-2] <- ch;
        ltys[k-2] <- lt;
      }
      if(lt > 0 && ch == 0) {
        lines(this.x, this.y, type = "l", lty = lt,
              xaxt="n", yaxt="n", col=my.colors[color]);
        pchs[k-2] = " ";
        ltys[k-2] = lt;
      }
      if(lt > 0 && ch > 0) {
        lines(this.x, this.y, type="b", lty = lt, pch=ch,
              xaxt="n", yaxt="n", col=my.colors[color]);
        pchs[k-2] <- ch;
        ltys[k-2] <- lt;
      }
      
      if(ch == 0 && lt == 0) {
        print("Skipping this line");
      }
    }
    markers<-as.character(lods[1:(rows-2),1]);

    # check if too markers are very close together
    lbl<-markers[markers != "-"];
    lloc<-dist[markers != "-"];
    
    last.pos<-lloc[1]/(lloc[length(lloc)]-lloc[1]) * pl.wd;
    for (j in 2:length(lloc)) {
      diff<-lloc[j]/(lloc[length(lloc)]-lloc[1]) * pl.wd - last.pos;
      if(diff < 0.1) {
        lbl[j]<- "    ";
      }
      else {
        last.pos<-last.pos+diff;
      }
    }

    if (length(lgnd) > 1) {
      draw.lgnd = FALSE      
      for (lg in lgnd) {
        if (lg == "page") {
          draw.lgnd = ((i == 1)  || (i %% (col*row)) == 1)
        }
        else {
          if (i == lg) {
            draw.lgnd = TRUE
            break
          }
        }
      }
    }
    else {
      draw.lgnd = lgnd
    }

    if (draw.lgnd == "page") {
      draw.lgnd = ((i == 1)  || (i %% (col*row)) == 1)
    }
    
    if (draw.lgnd == TRUE) {
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

      if (bw == F) {
        leg.color<- my.colors[1:(ncol(lods)-2)];
      }
      else {
        leg.color <- "black";
      }
      legend(lgx, lgy, leg, col=leg.color, lty=ltys, pch=pchs, cex=0.7);
    }
    
    axis(3, tck=0.05, at=lloc,labels=lbl,cex.axis=0.7,las=2);

#   detach(lods);
      
    if(output == "screen") {
      if((i %% (col*row)) == 0) {
        par(ask = TRUE);
      }
    }
  }
  if(batch == TRUE && output != "screen") {
    dev.off(dev.list());
  }
}

