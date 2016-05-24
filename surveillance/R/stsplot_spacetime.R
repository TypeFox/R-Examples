################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Old implementation of (animated) maps of an sts-object
###
### Copyright (C) 2007-2013 Michael Hoehle, 2016 Sebastian Meyer
### $Revision: 1694 $
### $Date: 2016-04-02 22:37:52 +0200 (Sam, 02. Apr 2016) $
################################################################################


stsplot_spacetime <- function(
    x, type, legend=NULL, opts.col=NULL, labels=TRUE,
    wait.ms=250, cex.lab=0.7, verbose=FALSE, dev.printer=NULL, ...)
{
  #Extract the mappoly
  if (length(x@map) == 0)
    stop("The sts object has an empty map.")
  map <- x@map
  maplim <- list(x=bbox(map)[1,],y=bbox(map)[2,])

  #Check colnames, otherwise no need to continue
  if (is.null(colnames(x@observed)))
    stop("The sts observed slot does not have any colnames to match with the shapefile.")

  #Check for legend options
  if (is.null(legend)) {
    legend <- list(dx=0.4,dy=0.04,x=maplim$x[1],y=maplim$y[1],once=TRUE)
  }

  #Extract the data
  o <- x@observed
  alarm <- x@alarm
  
  #Formula is of type "observed ~ 1|unit" (i.e. no time)
  aggregate <- type[[3]][[3]] == "unit"
  if (aggregate) {
    o <- t(as.matrix(apply(o,MARGIN=2,sum)))
    alarm <- t(as.matrix(apply(alarm,MARGIN=2,sum)))>0
  }
  
  #Number of time points
  maxt <- dim(o)[1]

  #Process dev.printer options
  if (is.list(dev.printer)) {
    dev.printer <- modifyList(
      list(device = png, extension = ".png", width = 640, height = 480,
           name = "Rplot"),
      dev.printer)
    #filename format (padding with zeroes)
    fnfmt <- paste0("%s-%0", nchar(maxt), "i%s")
  }

  #Get color vector
  opts.col_default <- list(ncolors=length(o), use.color=TRUE)
  gyr <- do.call("hcl.colors", if (is.list(opts.col))
    modifyList(opts.col_default, opts.col) else opts.col_default)
  theCut <- cut(o, length(gyr))
  
  #Cut into specified number of colors
  o.cut <- matrix(as.numeric(theCut),nrow=nrow(o),ncol=ncol(o))
  o.col <- matrix(gyr[o.cut],ncol=ncol(o.cut))
  o.col[is.na(o.col)] <- gray(1)
  dimnames(o.col) <- dimnames(o)

  #Sort the o according to the names in the map
  region.id <- row.names(map)
  o.col.id <- dimnames(o.col)[[2]]

  #Make the columns of o as in the map object
  o.col <- o.col[,pmatch(region.id,o.col.id),drop=FALSE]
  alarm.col <- alarm[,pmatch(region.id,o.col.id),drop=FALSE]

  #Screen processing
  screen.matrix <- matrix(c(0,1,0,1,0,1,0.8,1),2,4,byrow=TRUE)
  split.screen(screen.matrix)

  #Loop over all time slices
  for (t in 1:maxt) {
    #Status information
    if (verbose) {
      cat(paste("Processing slice",t,"of",maxt,"\n"))
    }
    
    #Clean screen (title area)
    screen(n=2)
    par(bg=gray(1))
    erase.screen()
    par(bg="transparent")

    #Plot the map on screen 1
    screen(n=1)
    plot(map,col=o.col[t,],xlab="",ylab="",...)
    #Indicate alarms as shaded overlays
    if (!all(is.na(alarm.col))) {
      #Plotting using density "NA" does not appear to work
      #anymore in the new sp versions
      alarm.col[is.na(alarm.col)] <- 0
      plot(map,dens=alarm.col*15,add=TRUE)
    }
    

    if (labels)
      #getSpPPolygonsLabptSlots is deprecated. Use coordinates method insteas
      text(coordinates(map), labels=as.character(region.id), cex.lab=cex.lab)
  
    if (!aggregate) { title(paste(t,"/",maxt,sep="")) }

    #In case a legend is requested
    if (is.list(legend) && !(legend$once & t>1)  | (t==1)) {
      add.legend(legend, maplim,
                 list(col=gyr, min=min(o), max=max(o), trans=identity))
    }

    #Is writing to files requested?
    if (is.list(dev.printer)) {
      #Create filename
      fileName <- sprintf(fnfmt, dev.printer$name, t, dev.printer$extension)
      cat("Creating ",fileName,"\n")
      #Save the current device using dev.print
      if (inherits(try(
        dev.print(dev.printer$device, file=fileName,
                  width=dev.printer$width, height=dev.printer$height)
      ), "try-error")) {
        warning("disabling dev.print()", immediate. = TRUE)
        dev.printer <- NULL
      }
    }
    
    wait(wait.ms) 
  }
  close.screen(all.screens = TRUE)
}



#######################
### auxiliary functions
#######################


### wait a specific amount of milliseconds (via "while" and "proc.time")

wait <- function (wait.ms) # number of milliseconds to wait
{
  #Initialize
  start.time <- proc.time()[3]*1000
  ellapsed <- proc.time()[3]*1000 - start.time

  #Loop as long as required.
  while (ellapsed < wait.ms) {
    ellapsed <- proc.time()[3]*1000 - start.time
  }
}


### add the color key

add.legend <- function(legend, maplim, theColors)
{
  #Preproc
  dy <- diff(maplim$y) * legend$dy
  dx <- diff(maplim$x) * legend$dx
    
  #Add legend -- i.e. a slider
  xlu <- xlo <- legend$x
  xru <- xro <- xlu + dx 
  yru <- ylu <- legend$y
  yro <- ylo <- yru + dy 

  
  step <- (xru - xlu)/length(theColors$col)
  for (i in 0:(length(theColors$col) - 1)) {
    polygon(c(xlo + step * i, xlo + step * (i + 1), 
              xlu + step * (i + 1), xlu + step * i), c(ylo, 
                                                       yro, yru, ylu), col = theColors$col[i + 1], 
            border = theColors$col[i + 1])
  }
  
  
  #Write info about min and max on the slider.
  black <- grey(0)
  lines(c(xlo, xro, xru, xlu, xlo), c(ylo, yro, yru, ylu, ylo), col =   black)

  #Transformation function for data values, e.g., exp or identity
  trans <- theColors$trans

  text(xlu, ylu - 0.5*dy, formatC(trans(theColors$min)), cex = 1, col = black,adj=c(0,1))
  text(xru, yru - 0.5*dy, formatC(trans(theColors$max)), cex = 1, col = black,adj=c(1,1))
}
