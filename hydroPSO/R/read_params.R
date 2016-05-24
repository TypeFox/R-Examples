# File read_params.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2008-2011 Mauricio Zambrano-Bigiarini
# Distributed under GPL 2 or later

################################################################################
#                                'read_params'                                 #
################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #  
# Started: June 09th, 2010,                                                    #
# Updates: 03-Dec-2010 ; 20-Dec-2010;                                          #
#          13-Jan-2011 ; 01-Jul-2011 ; 04-Nov-2011                             # 
#          19-Jan-2012 ; 02-Feb-2012 ; 23-Feb-2012 ; 07-Mar-2012               #
################################################################################
# This function reads the files:
# -) 'modelpara.out'         created by the GLUE algorithm of SWAT-CUP', 
# -) 'modelpara.beh'         created by the GLUE algorithm of SWAT-CUP', 
# -) 'goal.sf2'              created by the SUFI-2 algorithm of SWAT-CUP (http://www.eawag.ch/organisation/abteilungen/siam/software/swat/index_EN)
# -) 'goal.pso'              created by the PSO algorithm of SWAT-CUP
# -) 'ParameterValues.log'   created by Nimbus (Lisflood model)

# file   : the name of the file which the data are to be read from. 
#          Each row of the table appears as one line of the file. 
#          If it does not contain an absolute path, the file name is relative to 
#          the current working directory, getwd(). Tilde-expansion is performed 
#          here supported. As from R  2.10.0 this can be a compressed file (see file). 
# skip   : integer: the number of lines of the data file to skip before beginning to read data.
# header : a logical value indicating whether the file contains the names of the 
#          variables as its first line. If missing, the value is determined from 
#          the file format: header is set to TRUE if and only if the first row 
#          contains one fewer field than the number of columns.
# param.cols: numerical vector, with the column position in 'param.fname' that 
#             represent the parameter values used during the calibrations
# param.names: character, with the names to be given to each one of the parameters
#              stored in the columns represented by 'param.cols' (usually, the 
#              most sensitive parameters in 'param.fname' )
# of.col : numeric, with the column position  in 'param.fname' that represents
#          the objective function values obtained during the calibration                   
# of.name: character, with the name that will be given to the column 'of.col'
# plot   : logical, indicating if dotty plots of the parameter sets and its
#          objective function have to be plotted or not
# beh.thr: OPTIONAL, only used when 'plot=TRUE'. \cr
#          numeric, with a behavioural threshold to be used for plotting 
#          a horizontal line 
# beh.col: OPTIONAL, only used when 'plot=TRUE'. \cr
#          color to be used for plotting the horizontal line on 'beh.thr'
# beh.lty: OPTIONAL, only used when 'plot=TRUE'. \cr
#          type of line to be used for plotting the horizontal line on 'beh.thr'
# beh.lwd: OPTIONAL, only used when 'plot=TRUE'. \cr
#          line width to be used for plotting the horizontal line on 'beh.thr'
# nrow   : OPTIONAL, only used when 'plot=TRUE'. \cr
#          Number of rows to be used in the plotting window. The number 
#          of columns is automatically computed depending on the number 
#          of columns of 'x'
# main   : OPTIONAL, only used when 'plot=TRUE'. \cr
#          an overall title for the plot: see 'title'. 
# ylab   : OPTIONAL, only used when 'plot=TRUE'. \cr
#          a title for the y axis: see 'title'.
# col    : OPTIONAL, only used when 'plot=TRUE'. \cr
#          color to be used for plotting the points corresponding to the 
#          value of each parameter 
# cex    : OPTIONAL, only used when 'plot=TRUE'. \cr
#          A numerical value giving the amount by which plotting text
#          and symbols should be magnified relative to the default. See '?par'
# pch    : OPTIONAL, only used when 'plot=TRUE'. \cr
#          Either an integer specifying a symbol or a single character
#          to be used as the default in plotting points. See 'points'
#          for possible values and their interpretation. 
# ...    : ONLY used when 'plot=TRUE'. Arguments to be passed to methods, such as graphical parameters (see '?par').     

read_params <- function(file, ...) UseMethod("read_params")
                      
read_params.default <- function(file,
                       header=TRUE,
                       skip=0,                       
                       param.cols=NULL,
                       param.names=NULL,
                       of.col=NULL, 
                       of.name="GoF",
                       na.strings="-9999",
                       ####################
                       plot=TRUE,
                       ptype=c("histogram", "dottyplot", "boxplot", "vioplot", "pairs"),
                       MinMax=NULL,
                       beh.thr=NA, 
                       beh.col="red", 
                       beh.lty=1, 
                       beh.lwd=2, 
                       nrows="auto",
                       #col="black", 
                       col="#00000030",
                       ylab=of.name, 
                       main=NULL,
                       pch=19, 
                       cex=0.5,  
                       cex.main=1.5, 
                       cex.axis=1.5, 
                       cex.lab=1.5,
                       breaks="Scott",
                       freq=TRUE,
                       verbose=TRUE,
                       ...,
                       do.png=FALSE, 
                       png.width=1500, 
                       png.height=900,
                       png.res=90, 
                       png.fname="Parameters.png" 
                       ) {
                         
                         
  ##############################################################################
  # 1)                            Checkings                                    #
  ##############################################################################

  # Setting 'ptype' 
  ptype <- match.arg(ptype) 

  # Checking that 'file' exists
  if ( missing(file) ) {
    stop( "Missing argument: 'file'" )
  } else if ( !file.exists(file) )
     stop( "Invalid argument value: The file '", basename(file), "' doesn't exist" )

  # Checking that 'param.cols' exists
  if ( is.null(param.cols) ) stop( "Missing argument: 'param.cols' mmust be provided" )
  
  if (header==FALSE) {
    # Checking 'param.names' 
    if ( is.null(param.names) ) 
      stop( "Missing argument: 'param.names' has to be provided when 'header=FALSE' !!" )

    # Checking 'of.name' 
    if ( missing(of.name) ) 
      stop( "Missing argument: 'of.name' has to be provided when 'header=FALSE' !!" )
  } # IF end

  # Checking that 'param.cols' and 'param.names' have the same lenght
  if ( ( !is.null(param.cols) ) &  !is.null(param.names) ) {
    if ( length(param.cols) != length(param.names) )
      stop( "Invalid argument: length(param.cols) != length(param.names)' (", 
            length(param.cols), "!=", length(param.names), ")"  )
  } # IF end

  # Checking 'beh.thr'
  if ( !is.na(beh.thr) ) {
     if ( is.null(MinMax) )
        stop("Missing argument: 'MinMax' has to be provided before using 'beh.thr' !!")        
      if ( is.null(of.col) )
        stop("Missing argument: 'of.col' has to be provided before using 'beh.thr' !!")
  } # IF end
  
  # Checking 'MinMax'
  if ( !is.null(MinMax) ) {
     if ( !(MinMax %in% c("min", "max")) )
        stop("Invalid argument: 'MinMax' must be in c('min', 'max')")
  } # IF end

  # If the file is one of c("goal.sf2", "modelpara.beh", "modelpara.out", "goal.pso", "ParameterValues.log")
  if ( !is.na(match(basename(file), c("goal.sf2", "modelpara.beh", "modelpara.out", "goal.pso", "ParameterValues.log") ) ) ) {     
   
     if ( basename(file) %in% c("goal.sf2", "goal.pso") ) { 
       lskip   <- 3
       lheader <- TRUE
       #of.name  = "goal_value"
       #algorithm= "SUFI2" OR "PSO"
     } else if ( basename(file) %in% c("modelpara.beh", "modelpara.out") ) { # For GLUE, when basename(file) in c("modelpara.beh", "modelpara.out")
           lskip   <- 0
           lheader <- TRUE
           #of.name  = "objfun"
           #algorithm= "GLUE"
       } else if ( basename(file)=="ParameterValues.log" ) { # For NIMBUS
           lskip   = 1
           lheader <- FALSE
           #algorithm= "NIMBUS"
       } #IF end
  } else { # If the file IS NOT in c("goal.sf2", "modelpara.beh", "modelpara.out", "goal.pso", "ParameterValues.log")

       lskip   <- skip
       lheader <- header

    } # ELSE end
  
  ##############################################################################
  # 2)                Reading ALL the PARAMETER SETS                           #
  ##############################################################################
  if (verbose) message( "                                                     ")  
  if (verbose) message( "[ Reading the file '", basename(file), "' ... ]" )  
  params <- read.table(file=file, header=lheader, skip=lskip) 

  # Amount of total parameter sets 
  nparamsets <- nrow(params)
  if (verbose) message( "[ Total number of parameter sets: ", nparamsets, " ]" )     
  
  # Giving the right name to the columns with the parameter values
  if (!is.null(param.names)) {
    colnames(params)[param.cols] <- param.names
  } else param.names <- colnames(params)[param.cols]

  # Keeping only the columns with parameter values
  params.only <- params[, param.cols] 

  # computing the number of parameters
  nparams <- ncol(params.only)

  # Keeping only the columns with the goodness-of-fit measure, if provided
  if ( !is.null(of.col) ) {
    gofs   <- params[, of.col] 
  } else gofs   <- NULL
  
  # Filtering out those parameter sets above/below a certain threshold
  if (!is.na(beh.thr)) {  
     # Checking 'beh.thr'
     mx <- max(gofs, na.rm=TRUE)
     if (beh.thr > mx)
       stop("Invalid argument: 'beh.thr' must be lower than ", mx ,"!!")
    
    # Computing the row index of the behavioural parameter sets
    ifelse(MinMax=="min", beh.row.index <- which(gofs < beh.thr), 
                          beh.row.index <- which(gofs > beh.thr) )
    
    # Removing non-behavioural parameter sets and gofs
    params.only <- params.only[beh.row.index, ]
    gofs        <- gofs[beh.row.index]
   
    # Amount of behavioural parameter sets 
    nbeh <- nrow(params.only)
    if (verbose) message( "[ Number of behavioural parameter sets: ", nbeh, " ]" )
    
    # To avoid problems with 'plot_params'
    if (plot) beh.thr <- NA
  } # IF end
  
  ##############################################################################
  # 3)                            Plotting                                     #
  ##############################################################################
  if (plot) {
  
    plot_params(params=params.only, 
                gofs=gofs,
                ptype=ptype,
                MinMax=MinMax,
                param.cols=1:nparams,
                param.names=param.names,
                ##of.col=of.col,
                of.name=of.name, 
                beh.thr=beh.thr, 
                beh.col=beh.col, 
                beh.lty=beh.lty, 
                beh.lwd=beh.lwd, 
                nrows=nrows,
                col=col, 
                ylab=ylab, 
                main=main, 
                pch=pch, 
                cex=cex, 
                cex.axis=cex.axis, 
                cex.main=cex.main, 
                cex.lab=cex.lab, 
                breaks=breaks, 
                freq=freq, 
                verbose=verbose,
                ..., 
                do.png=do.png, 
                png.width=png.width, 
                png.height=png.height,
                png.res=png.res, 
                png.fname=png.fname
                )  
  
  } # IF end
  
  # Creating the output
  out <- list(params=params.only, gofs=gofs)  
  return(out)
  
}  # 'read_param.default' END


################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #  
################################################################################
# Started: 23-Feb-2012                                                         #
# Updates:                                                                     #
################################################################################
#read_params.data.frame <- function(file,
#                       header=TRUE,
#                       skip=0,                       
#                       param.cols=NULL,
#                       param.names=NULL,
#                       of.col=NULL, 
#                       of.name="GoF",
#                       na.strings="-9999",
#                       ####################
#                       plot=TRUE,
#                       ptype=c("histogram", "dottyplot", "boxplot", "vioplot", "pairs"),
#                       MinMax=NULL,
#                       beh.thr=NA, 
#                       beh.col="red", 
#                       beh.lty=1, 
#                       beh.lwd=2, 
#                       nrows="auto",
#                       #col="black", 
#                       col="#00000030",
#                       ylab=of.name, 
#                       main=NULL,
#                       pch=19, 
#                       cex=0.5,  
#                       cex.main=1.5, 
#                       cex.axis=1.5, 
#                       cex.lab=1.5,
#                       breaks="FD",
#                       freq=FALSE,
#                       verbose=TRUE,
#                       ...,
#                       do.png=FALSE, 
#                       png.width=1500, 
#                       png.height=900,
#                       png.res=90, 
#                       png.fname="Parameters.png" 
#                       ) {
#                        
#    read_params.default(file=file,
#                       header=header,
#                       skip=skip,                       
#                       param.cols=param.cols,
#                       param.names=param.names,
#                       of.col=of.col, 
#                       of.name=of.name,
#                       na.strings=na.strings,
#                       ####################
#                       plot=plot,
#                       ptype=ptype,
#                       MinMax=MinMax,
#                       beh.thr=beh.thr, 
#                       beh.col=beh.col, 
#                       beh.lty=beh.lty, 
#                       beh.lwd=beh.lwd, 
#                       nrows=nrows,
#                       #col="black", 
#                       col=col,
#                       ylab=ylab, 
#                       main=main,
#                       pch=pch, 
#                       cex=cex,  
#                       cex.main=cex.main, 
#                       cex.axis=cex.axis, 
#                       cex.lab=cex.lab,
#                       breaks=breaks,
#                       freq=freq,
#                       verbose=verbose,
#                       ...,
#                       do.png=do.png, 
#                       png.width=png.width, 
#                       png.height=png.height,
#                       png.res=png.res, 
#                       png.fname=png.fname 
#                        )
#                        
#} # 'read_params.data.frame' END



