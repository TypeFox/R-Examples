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
# part of NPLPLOT library
# Calls nplplot()
# Allows the user to specify multiple plots within the same output file
# Usage:
#       nplplot.multi(files, col, row, mode, output, headerfiles, lgnd)
#
# Allows upto 6 colors for curves inside the plot.
#
# Taken out the plotting commands.
# This is a more hands-free, therefore restrictive function,
# The file format has been changed. No point or line-types, take these
# from the header file, or set default values.
#
# Allows user to select where to add legends, first plot on every page etc.
#
# headerfile contains :
# one value each for lgnd, title, ymin, ymax, xlabel, y label, 
# legends, linetypes, pointypes, colors (as many as plots).
# This function does not allow the ... argument ?
#

nplplot.multi<-function(filenames, plotdata=NULL, col=2, row=2, mode="l", output="screen",
                        headerfiles=NULL, lgnd="page",
                        customtracks=FALSE, mega2mapfile=NULL, 
                        pagewidth=NULL, pageheight=NULL, topmargin=0.25, ...)
{
  if (!is.null(plotdata))
    numfiles=length(plotdata)
  else
    numfiles<-length(filenames)
  
  # page size in inches
  if (is.null(pagewidth)) {
    pagewidth<- 8.0
  }
  
  if (is.null(pageheight)) {
    pageheight<- 10.5
  }
  # flip if necessary
  if(mode=="p") {
    full.width<-pagewidth
    full.height<-pageheight
    landscape<-FALSE
  }
  else {
    full.width<- pageheight
    full.height<- pagewidth
    landscape<-TRUE
  }

  pl.wd <- full.width/col
  
  if(output == "screen") {
    do.call(names(dev.cur()), list())
  }
  else {
    output1 <- output
    fileparts <- unlist(strsplit(output1, "\\."))
    if (fileparts[length(fileparts)] == "ps") {
      print("Creating postscript file")
      postscript(width=full.width, height=full.height,
                 file=output, paper="letter", horizontal=landscape)
    }
    else {
      print("Creating pdf file")
      pdf(width=full.width, height=full.height,
          file=output, paper="special")
    }
  }

  par(omi=c(0.05, 0.05, topmargin, 0.05));
  par(mfrow = c(row, col), ask=FALSE)

  retval <- TRUE;

  if (length(lgnd) > 1) {
    # test all elements
    for (lg in lgnd) {
      if (!is.numeric(lg)) {
        print("Invalid list value for lgnd, use plot numbers.")
        print("No legends will be drawn.")
        lgnd <- "page"
        break
      }
    }
  }
  else {
    if (!is.logical(lgnd) && (lgnd != "page") && !is.numeric(lgnd)) {
      print("Invalid value for lgnd, use PAGE, TRUE/FALSE or plot numbers.")
      print("Setting to PAGE")
      lgnd <- "page"
    }
  }
  
  if ((customtracks==TRUE) && !is.null(mega2mapfile)) {
    chrlist <- NULL;
  }
  
  for (i in 1:numfiles) {
     # read the header file to set the arguments to nplplot,
     # then call nplplot

    # If the user has supplied headerfiles, then read in a
    # line from each file as header
    # There should be as many header files as there are data files
    # (number of plots). If there are fewer header files, recycle
    # the last one.
    
    title<- ""
    yline<-2.0;     ymin<-0.0;     ymax<- 4.0;  yfix <- FALSE
    ltypes<-NULL;   ptypes<- NULL;     my.colors<-NULL
    bw <- TRUE 
    lgndx<-NULL;  lgndy <- NULL;    lgndtxt <- NULL;  cex.legend <- 0.7
    ylabel<- "Scores"
    tcl<- 0.3
    cex.axis <- 0.9
    draw.lgnd <- TRUE
    
    if(!is.null(headerfiles)) {
      if (i <= length(headerfiles)) {
        hdrfile <- headerfiles[i]
      }	     
      else {
        hdrfile <- headerfiles[length(headerfiles)]
      }
      
      # parse headerfile
      lines <- readLines(con=hdrfile, n=-1)
      
      for (l in 1: length(lines)) {
        eval(parse(text=lines[l]))
      }
    }

    if (length(lgnd) == 1) {
      if (i == lgnd) {
        draw.lgnd <- TRUE
      }
      else if (lgnd == "page") {
        draw.lgnd <- ((i == 1)  || (i %% (col*row)) == 1)
      }
      else if (is.logical(lgnd)) {
        draw.lgnd <- lgnd
      }
      else {
        draw.lgnd <- FALSE
      }
    }
    else {
      for (lg in lgnd) {
        if (is.numeric(lg) && i == lg) {
          draw.lgnd <- TRUE
          break
        }          
        else {
          draw.lgnd <- FALSE
        }
      }
    }

    # assume files with a header,

    if (!is.null(plotdata)) {
        title=paste("Chromosome", names(plotdata[i]))

        retval1 <- nplplot(filename=NULL, plotdata=plotdata[[i]], header=TRUE, yline=yline, ymin=ymin,
                           ymax=ymax, yfix=yfix, title=title, draw.lgnd=draw.lgnd,
                           xlabl="", ylabl=ylabel,
                           lgndx=lgndx, lgndy=lgndy, lgndtxt=lgndtxt, cex.legend=cex.legend,
                           bw=bw, my.colors=my.colors, ltypes=ltypes, ptypes=ptypes,
                           tcl=tcl, cex.axis=cex.axis,
                           plot.width=pl.wd, ...)
      } else {

        retval1 <- nplplot(filename=filenames[i], plotdata=NULL, header=TRUE, yline=yline, ymin=ymin,
                           ymax=ymax, yfix=yfix, title=title, draw.lgnd=draw.lgnd,
                           xlabl="", ylabl=ylabel,
                           lgndx=lgndx, lgndy=lgndy, lgndtxt=lgndtxt, cex.legend=cex.legend,
                           bw=bw, my.colors=my.colors, ltypes=ltypes, ptypes=ptypes,
                           tcl=tcl, cex.axis=cex.axis,
                           plot.width=pl.wd, ...)
      }
    
    retval <- retval && retval1
  
    if ((customtracks==TRUE) && !is.null(mega2mapfile)) {
      # first create the prefix by taking out the last extension
      # Check that each of these are numeric or 'X', 'XY', 'Y'
      if (is.null(plotdata)) {
          a<-unlist(strsplit(filenames[i], ".", fixed=TRUE))
          chrlist[i] <- switch(a[length(a)],
                               X = "X",
                               Y = "Y",
                               XY = "XY",
                               a[length(a)])
          if (i==1) {
              prefix <- paste(a[1:length(a)-1], sep=".")
          }
      } else {
          chrlist[i] = names(plotdata[i])
          if (i==1) {
              a<-unlist(strsplit(filenames[i], ".", fixed=TRUE))
              prefix <- paste(a[1:length(a)-1], sep=".")
          }
      }
    }
    
    if(output == "screen") {
      if((i %% (col*row)) == 0) {
        par(ask = TRUE)
      }
    }
      
  }
  
  
  if ((customtracks == TRUE) && !is.null(mega2mapfile)) {
#   prepareplot(prefix, chrlist,  mega2mapfile, output="both")
    prepareplot(plotdata, chrlist,  mega2mapfile, output="both")
    for (i in (1:numfiles)) {
      if (!is.null(plotdata))
        if (nchar(chrlist[i]) == 1) chrlist[i] = paste("0", chrlist[i], sep="")
      bedplot(paste("bed.data", chrlist[i], sep="."))
    }
    retval1 <- genomeplot("GG.data.all")
    if (retval1 == FALSE) {
      print("Unable to create BED plot files and Genome Graph plot files.")
    }
  }               


  if(names(dev.cur()) == "postscript" || names(dev.cur()) == "pdf") {
    dev.off()     
  }

  if (retval) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

  
