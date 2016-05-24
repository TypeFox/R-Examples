#
# OPI for Octopus 900 
# 
# This would all have been nicer in an OO style, with each implementation
# being a subclass of an opi class, but I don't think it can be in R.
# The OPI standard doesn't want users employing exactly the same function 
# no matter what the underlying implementation, and so there cannot be 
# extra parameters to create method signatures for different classes.
# Similarly, some implementations use exactly the same method signatures,
# again which will confuse R, I think. Anyway, if I am wrong, sorry about that.
# What I've done (use a list of implementations and then use a global
# integer to index them) works and makes sense to the non-OO person.
#
# This version creates a socket to O900Server.java and sends/receives 
# commands. It requires the server to be running, and the server 
# to know where the H-S jar files are (ie in the CLASSPATH). 
# 
# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
# Author: David Lawson    (XXX)
# Date: July 2014
#
# Copyright 2012 Andrew Turpin and David Lawson
# This program is part of the OPI (http://perimetry.org/OPI).
# OPI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


###################################################################
# .Octopus900Env$octopusObject is the java Opi object set in opiInitialize
# .Octopus900Env$... are lots of colors and Fixation constants set in setupBackgroundConstants()
###################################################################
if (!exists(".Octopus900Env"))
    .Octopus900Env <- new.env()

###########################################################################
# Get values for fixation, color and bg intensity constants
# from EyeSuite classes, and set globals
#       .Octopus900Env$* 
# to the values of those constants.
# INPUT: None.
# OUTPUT: None.
# SIDE EFFECTS: sets .Octopus900Env$* if possible.
#               store the names and values in .Octopus900Env$constList 
###########################################################################
setupBackgroundConstants <- function() {

    constList <- NULL

    getC <- function(cName) {
        writeLines(paste("OPI_GET_CONSTANT", cName), .Octopus900Env$socket)
        res <- readLines(.Octopus900Env$socket, n=1)
        if (res == "OZ900Fail") {
            warning(paste("Cannot set",cName,"constant for the O900."))
        } else {
            assign(cName, as.double(res), envir = .Octopus900Env) 
            assign("constList", c(.Octopus900Env$constList, list(list(cName, as.double(res)))), 
                    envir = .Octopus900Env)
        }
    }

    getC("FIX_CENTRE")
    getC("FIX_CROSS")
    getC("FIX_RING")
    getC("BG_OFF")
    getC("BG_1")        # 1.27 cd/m2 == 127 passed to MsgInitializePerimUnit
    getC("BG_10")       # 10 cd/m2 == 1000
    getC("BG_100")      # 100 cd/m2 == 10000

    assign("FIX_CENTER", .Octopus900Env$FIX_CENTRE, envir = .Octopus900Env) # help Americans

        # get the color fields from OCTO900
    getC("STIM_WHITE")
    getC("STIM_BLUE")
    getC("STIM_RED")
    getC("BG_WHITE")
    getC("BG_YELLOW")
    getC("MET_COL_WW")
    getC("MET_COL_BY")
    getC("MET_COL_RW")
    getC("MET_COL_BLUE_WHITE")
    getC("MET_COL_RED_YELLOW")
    getC("MET_COL_WHITE_YELLOW")
    getC("MET_COL_USER")

    assign("MET_COL_BW", .Octopus900Env$MET_COL_BLUE_WHITE,   envir = .Octopus900Env)
    assign("MET_COL_RY", .Octopus900Env$MET_COL_RED_YELLOW,   envir = .Octopus900Env)
    assign("MET_COL_WY", .Octopus900Env$MET_COL_WHITE_YELLOW, envir = .Octopus900Env)
}

#######################################################################
# INPUT: 
#   serverPort               = port number on which server is listening
#   eyeSuiteSettingsLocation = dir name containing EyeSuite settings
#   eye                      = "right" or "left"
#   gazeFeed                 = 0 (none), 1 (single frame), 2 (all frames with *)
#   bigWheel                 = FALSE (standard machine), TRUE for modified apeture wheel
#   buzzer                   = 0 (no buzzer),1,2, 3 (max volume)
#
#   Both input dirs should INCLUDE THE TRAILING SLASH.
#
# @return NULL if succeed
# @return 1 if already initialised
# @return 2 if failed to make ready
#
#######################################################################
octo900.opiInitialize <- function(serverPort=50001,eyeSuiteSettingsLocation=NA, 
                                  eye=NA, gazeFeed=0, bigWheel=FALSE, buzzer=0) {
    assign("gazeFeed", gazeFeed, envir=.Octopus900Env)

    if (!bigWheel) {
        assign("GOLDMANN", c(6.5, 13, 26, 52, 104) / 60, envir=.Octopus900Env)
    } else {
        mm <- c(0.125,0.25,0.5,1,1.41,2,2.83,4,5.66,8,11.3,16,22.6,32,64,128,256)
        ind <- c(32,28,31,26,30,29,27,24,25,23,21,22,39,38,20,37,36)
        GOLDMANN <- rep(NA,39)
        GOLDMANN[ind] <- (sqrt(mm/pi)*180/pi/149.1954)
        assign("GOLDMANN", GOLDMANN, envir=.Octopus900Env)
    }

    if (is.na(buzzer) || buzzer < 0) buzzer <- 0
    if (buzzer > 3) buzzer <- 3

    cat("Looking for server... ")
    suppressWarnings(tryCatch(    
        v <- socketConnection(host = "localhost", serverPort,
                      blocking = TRUE, open = "w+b",
                      timeout = 10)
        , error=function(e) { 
            stop(paste(" cannot find a server on port",serverPort))
        }
    ))
    close(v)
    
    print("found server :)")

    if (is.na(eyeSuiteSettingsLocation))
        stop("You must specify the EyeSuite settings folder in your call to opiInitialize")
    if (is.na(eye))
        stop("You must specify which eye ('left' or 'right') in your call to opiInitialize")
    if (eye != "left" && eye != "right")
        stop("The eye argument of opiInitialize must be 'left' or 'right'")

    socket <- tryCatch(
        socketConnection(host="localhost", serverPort, open = "w+b", blocking = TRUE, timeout = 1000), 
        error=function(e) stop(paste("Cannot connect to Octopus 900 on port", serverPort))
    )

    assign("socket", socket, envir = .Octopus900Env)
    msg <- paste0("OPI_INITIALIZE \"",eyeSuiteSettingsLocation,"\"\ ",eye, " ", gazeFeed, " ", buzzer)
    writeLines(msg, socket)
    res <- readLines(socket, n=1)
    
    setupBackgroundConstants()

    if (res == "0")
        return(NULL)
    else
        return(res)
}

###########################################################################
# INPUT: 
#   As per OPI spec
#   stim$color must be same as that initialised by opiSetBackground or opiInitialize
#   If F310 is FALSE, response is taken from internal button 
#        (sends OPI_PRESENT_STATIC to server)
#   If F310 is TRUE , response is taken from external controller
#        (sends OPI_PRESENT_STATIC_F310 to server)
#
# Return a list of 
#    err  = string message
#    seen = 1 if seen, 0 otherwise
#    time = reaction time
# 
# If stim is null, always return 0 status.
###########################################################################
octo900.presentStatic <- function(stim, nextStim, F310=FALSE) {
    if (is.null(stim)) 
        return(list(err=0))

    if(min(abs(.Octopus900Env$GOLDMANN - stim$size), na.rm=TRUE) != 0)
        warning("opiPresent: Rounding stimulus size to nearest Goldmann size")

    if (F310)
        msg <- "OPI_PRESENT_STATIC_F310 "
    else
        msg <- "OPI_PRESENT_STATIC "
    msg <- paste(msg, stim$x * 10.0, stim$y * 10.0, cdTodb(stim$level, 4000/pi) * 10.0)
    msg <- paste(msg, (which.min(abs(.Octopus900Env$GOLDMANN - stim$size))))
    msg <- paste(msg, stim$duration)
      msg <- paste(msg, stim$responseWindow)
    if (!is.null(nextStim)) {
        msg <- paste(msg, nextStim$x * 10.0, nextStim$y * 10.0)
    }

    writeLines(msg, .Octopus900Env$socket)
    #Sys.sleep(1)
    res <- readLines(.Octopus900Env$socket, n=1)
    s <- strsplit(res, "|||", fixed=TRUE)[[1]]
    if (s[1] == "null") {
      err <- NULL
    } else {
      err <- s[1]
    }



    if (.Octopus900Env$gazeFeed==0) {
      return(list(
        err=err,
        seen=strtoi(s[2]),
        time=strtoi(s[3]),
        frames=NA,
        numFrames=NA,
        width=NA,
        height=NA
      ))
    }#gazeFeed=0



    if (.Octopus900Env$gazeFeed==1) {
      return(list(
        err=err,
        seen=strtoi(s[2]),
        time=strtoi(s[3]),
        frames=NA,
        numFrames=NA,
        width=NA,
        height=NA
      ))
    }#gazeFeed=1



    if (.Octopus900Env$gazeFeed==2) {
      #frames <- strsplit(s[14], "###", fixed=TRUE)[[1]]
      #for (i in 1:length(frames)) {
      #  frames[i] <- strtoi(strsplit(frames[i],",",fixed=T)[[1]])
      #}
      return(list(
                   err=err, 
                  seen=strtoi(s[2]),
                  time=strtoi(s[3]),
             numFrames=strtoi(s[4]),
                 times=strtoi(strsplit(s[5],",",fixed=T)[[1]]),
                   ids=strtoi(strsplit(s[6],",",fixed=T)[[1]]),
                 stars=strtoi(strsplit(s[7],",",fixed=T)[[1]]),
                 rings=strtoi(strsplit(s[8],",",fixed=T)[[1]]),
              diamters=strtoi(strsplit(s[9],",",fixed=T)[[1]]),
                 pupilX=strtoi(strsplit(s[10],",",fixed=T)[[1]]),
                pupilY=strtoi(strsplit(s[11],",",fixed=T)[[1]]),
        pupilMajorAxis=strtoi(strsplit(s[12],",",fixed=T)[[1]]),
        pupilMinorAxis=strtoi(strsplit(s[13],",",fixed=T)[[1]]),
          firstFrame=strtoi(strsplit(s[14],",",fixed=T)[[1]])
                #frames=frames
      ))
    }#gazeFeed=2
}

###########################################################################
# Set up generic calls based on type of stim
###########################################################################
octo900.opiPresent <- function(stim, nextStim=NULL) { UseMethod("octo900.opiPresent") }
setGeneric("octo900.opiPresent")

octo900.opiPresentF310 <- function(stim, nextStim=NULL) { UseMethod("octo900.opiPresentF310") }
setGeneric("octo900.opiPresentF310")

octo900.opiPresent.opiStaticStimulus <- function(stim, nextStim) {
    return(octo900.presentStatic(stim, nextStim, FALSE))
}

octo900.opiPresentF310.opiStaticStimulus <- function(stim, nextStim) {
    return(octo900.presentStatic(stim, nextStim, TRUE))
}

 
###########################################################################
# INPUT: 
#   As per OPI spec
#   stim$color must be same as that initialised by opiSetBackground or opiInitialize
#
# Return a list of 
#    err  = string message
#    seen = 1 if seen, 0 otherwise
#    time = reaction time
#
# If stim is null, always return 0 status.
###########################################################################
octo900.opiPresent.opiTemporalStimulus <- function(stim, nextStim=NULL, ...) {
    if (is.null(stim)) 
        return(list(err=0))

    if(min(abs(.Octopus900Env$GOLDMANN - stim$size)) != 0)
        warning("opiPresent: Rounding stimulus size to nearest Goldmann size")

    msg <- "OPI_PRESENT_TEMPORAL "
    msg <- paste(c(msg, stim$x * 10.0, stim$y * 10.0, stim$rate), collapse=" ")
    msg <- paste(msg, (which.min(abs(.Octopus900Env$GOLDMANN - stim$size))))
    msg <- paste(msg, stim$duration)
    msg <- paste(msg, stim$responseWindow)
    if (!is.null(nextStim)) {
        msg <- paste(msg, nextStim$x * 10.0, nextStim$y * 10.0)
    }

    writeLines(msg, .Octopus900Env$socket)
    res <- readLines(.Octopus900Env$socket, n=1)
    s <- strsplit(res, "|||", fixed=TRUE)[[1]]

    if (s[1] == "null") {
      err <- NULL
    } else {
      err <- s[1]
    }

    return(list(
        err =err, 
        seen=strtoi(s[2]),
        time=strtoi(s[3]),
        frames=NA,
        numFrames=NA,
        width=NA,
        height=NA
    ))

}#opiPresent.opiTemporalStimulus()

########################################## 
# Present kinetic stim, return values 
#
# If stim is null, always return 0 status.
########################################## 
octo900.opiPresent.opiKineticStimulus <- function(stim, ...) {
    if (is.null(stim)) 
        return(list(err=0))

        # convert sizes to GOLDMANN
     stim$sizes <- sapply(stim$sizes, function(s) {
         i <- which.min(abs(.Octopus900Env$GOLDMANN - s))
         if(abs(.Octopus900Env$GOLDMANN[i] - s) > 0.000001) {
             warning(paste("opiPresent: Rounding stimulus size",s,"to nearest Goldmann size"))
         } 
         return(i)
     })

    msg <- "OPI_PRESENT_KINETIC "
    xs <- xy.coords(stim$path)$x
    ys <- xy.coords(stim$path)$y
    msg <- paste(c(msg, length(xs), xs, ys), collapse=" ")
    msg <- paste(c(msg, sapply(stim$levels, cdTodb, maxStim=4000/pi)), collapse=" ")
    msg <- paste(c(msg, stim$sizes), collapse=" ")
    
      # convert seconds/degree into total time for path segment in seconds
    pathLengths <- NULL
    for(i in 2:length(xs)) {
      d <- sqrt((xs[i]-xs[i-1])^2 + (ys[i]-ys[i-1]^2))
      stim$speeds[i-1] <- d/stim$speeds[i-1]
    }
    msg <- paste(c(msg, stim$speeds), collapse=" ")  
    
    writeLines(msg, .Octopus900Env$socket)
    res <- readLines(.Octopus900Env$socket, n=1)
    s <- strsplit(res, "|||", fixed=TRUE)[[1]]

    if (s[1] == "null") {
      err <- NULL
    } else {
      err <- s[1]
    }

    return(list(
        err =err, 
        seen=strtoi(s[2]),
        time=strtoi(s[3]),
        x=strtoi(s[4]),     
        y=strtoi(s[5])
    ))
}

###########################################################################
#
# Input paras are the Octopus900Env$* constants
# lum is in cd/m^2 (as per OPI spec) * 100 == .Octopus900Env$BG_{OFF | 1 | 10 | 100 }
# color is .Octopus900Env$MET_COL_{WW | BY | RW | BLUE_WHITE | RED_YELLOW | WHITE_YELLOW }
# fixation is .Octopus900Env$FIX_{RING | CROSS | CENTRE}
# fixIntensity is 0..100 %
#
# @return NULL is succeed.
# @return -1 if opiInitialize has not been successfully called
# @return -2 trouble setting backgound color
# @return -3 trouble setting fixation
###########################################################################
octo900.opiSetBackground <- function(lum="NA", color="NA", fixation="NA", fixIntensity=50) {

    msg <- paste("OPI_SET_BACKGROUND", color, lum, fixation, fixIntensity)
    writeLines(msg, .Octopus900Env$socket)
    ret <- strtoi(readLines(.Octopus900Env$socket, n=1))

    if (ret == 0) {
        return(NULL)
    } else {
        return(ret)
    }
}

###########################################################################
# return NULL on success (in fact, always!)
###########################################################################
octo900.opiClose <- function() {
    writeLines("OPI_CLOSE", .Octopus900Env$socket)
    close(.Octopus900Env$socket)
    return(NULL)
}

###########################################################################
# Lists defined constants
###########################################################################
octo900.opiQueryDevice <- function() {
    cat("Defined constants\n")
    cat("-----------------\n")
    lapply(.Octopus900Env$constList, function(x) {
      lapply(x, cat, " ")
      cat("\n")
    })

    return(NULL)
}
