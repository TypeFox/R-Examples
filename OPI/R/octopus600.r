#
# OPI for Octopus 600 
# 
# Authors: 
#   Andrew Turpin    (aturpin@unimelb.edu.au)
#   David Lawson     (david.lawson@unimelb.edu.au)
# Date: July 2014
#
# Copyright 2014 Andrew Turpin and David Lawson
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

if (!exists(".Octopus600Env"))
    .Octopus600Env <- new.env()

#######################################################################
# Networking helper functions
#######################################################################

send <- function(packet, socket) { 
  writeBin(packet, socket, size = 4, endian = "big")
}

pad <- function(data) { 
  padsize = 128 - length(data)
  data = append(data, rep(0, padsize)) 
  return(data)
}

sendCommand <- function(socket, id, ...) {
  packet = c(id, 0, 0, 10, 0, 0, 0)
  packet = append(packet, c(...))
  packet = pad(packet)
  packet = as.integer(packet)
  
  #print("Sending:")
  #print(packet)
  
  send(packet, socket)
  
  response = readBin(socket, integer(), n = 128, size = 4, signed = TRUE, endian = "big")
  
  #print("Received:")
  #print(response)
  
  if (length(response) == 0)
    stop("Octopus 600 not responding")

  if (response[3] != 0) {
    warning(paste("Octopus 600 returned error ID", response[3]))
  }
  
  return(list(response, response[3]))
}

#######################################################################
# INPUT: 
#   ipAddress                = IP address of Octopus 600 machine
#   eye                      = "right" or "left"
#   pupilTracking
#   pulsar
#   eyeControl
#
#
# @return NULL if succeed
# @return ...
#
#######################################################################
octo600.opiInitialize <- function(ipAddress, eye, pupilTracking=FALSE, pulsar=FALSE, eyeControl=0) {
  
  if (missing(ipAddress))
    stop("You must specify an IP address in opiInitialize()")
  
  if (missing(eye) || (eye != "left" && eye != "right"))
    stop("You must set eye=left or eye=right in opiInitialize()")
  
  if (!is.element(eyeControl, 0:3))
    stop("eyeControl must be equal to 0, 1, 2 or 3")
  
  socket = tryCatch(
    socketConnection(host = ipAddress, 50000, open = "w+b", blocking = TRUE, timeout = 5), 
    error = function(e) stop(paste("Cannot connect to Octopus 600 on", ipAddress))
  )
  
  assign("socket", socket, envir = .Octopus600Env)
  
  print("Connected to Octopus 600")
  
  # set_eyecontrol()
  res = sendCommand(.Octopus600Env$socket, 2005, eyeControl, 60, 47, 136, 75)
  if (res[[2]] != 0)
    return(res[[2]])
  
  # initialise_perimeter()
  res = sendCommand(.Octopus600Env$socket, 2001, ifelse(pulsar, 3183, 1000))
  if (res[[2]] != 0)
    return(res[[2]])
  
  print(paste("initialise_perimeter returned freqLeft =", res[[1]][9], "and freqRight =", res[[1]][10]))
  
  if (pupilTracking) {
    # set_ir_illumination()
    res = sendCommand(.Octopus600Env$socket, 2007, eye=="left", eye=="right")
    if (res[[2]] != 0)
      return(res[[2]])
  }
  
  # set_fixationmark()
  if (pulsar) {
    res = sendCommand(.Octopus600Env$socket, 2003, eye=="left", 1, 2, 255) # yellow dot
    if (res[[2]] != 0)
      return(res[[2]])
  } else {
    res = sendCommand(.Octopus600Env$socket, 2003, eye=="left", 2, 2, 255) # yellow cross
    if (res[[2]] != 0)
      return(res[[2]])
  }
  
  assign("pupilTrackingEnabled", pupilTracking, envir = .Octopus600Env)
  assign("pupilBlackLevelSet", !pupilTracking, envir = .Octopus600Env)
  assign("eye", eye, envir = .Octopus600Env)
  assign("pulsar", pulsar, envir = .Octopus600Env)
  
	return(NULL)
}

###########################################################################
# INPUT: 
#   As per OPI spec
#   stim$color must be same as that initialised by opiSetBackground or opiInitialize
#
# Return a list of 
#	err  = string message
#	seen = 1 if seen, 0 otherwise
#	time = reaction time
###########################################################################
octo600.opiPresent <- function(stim, nextStim=NULL) { UseMethod("octo600.opiPresent") }
setGeneric("octo600.opiPresent")

octo600.opiPresent.opiStaticStimulus <- function(stim, nextStim) {
  
  leftEye = .Octopus600Env$eye == "left"
  
  if (!.Octopus600Env$pupilBlackLevelSet) {
    # adjustPupilBlackLevel()
    res = sendCommand(.Octopus600Env$socket, 2029, leftEye)
    if (res[[2]] != 0)
      return(list(err = res[[2]], seen=NA, time=NA))
    else
      assign("pupilBlackLevelSet", TRUE, envir = .Octopus600Env)
  }
  
  # display_stimulus()
  res <- sendCommand(
    .Octopus600Env$socket, 2000,
    0, #checkBGIllumi [Do always set to 0]
    0, #BGIntensity [If checkBGIllumi is set to 0, don't care]
    stim$x*10, #positionX [in 1/10deg]
    stim$y*10, #positionY [in 1/10deg]
    .Octopus600Env$pulsar*5, #method [0 = White-On-White, 5 = pulsar]
    0, #color [don't care]
    3, #stimulusSize [don't care] (has to be 3 for W-on-W, don't care for pulsar)
    cdTodb(stim$level, 4000/pi)*10, #dLog (intensity) [in 1/10 dB]
    stim$duration, #duration [stimulus presentation duration in ms, for W/W 100ms, for pulsar 500ms]
    leftEye, #selectedEye [0 = OD, 1 = OS]
    stim$responseWindow, #maxAllowedReactionTime (maximal allowed reaction time in ms, >=500ms and <4s)
    0, #type [0 = present normal stimulus, 1 = present positive catch trial, 2 = present negative catch trial]
    #sound [0 = no sound; Bit0=1 sound for stimulus presentation ON;
    #Bit1=1 sound for patient response button ON; Bit2=1 sound for fixation lost ON]
    ifelse(is.element("sound", names(stim)), stim$sound, 0)
  )
  
  if (res[[2]] != 0)
    return(list(err = res[[2]], seen=NA, time=NA))
  
  pupilSize = res[[1]][21]
  reactionTimePAK = res[[1]][22]
  answer = res[[1]][23]

  return(list(
    err = 0, 
    seen = answer == 1,
    time = reactionTimePAK,
    answer = answer
  ))
  
}

octo600.opiPresent.opiTemporalStimulus <- function(stim, nextStim=NULL, ...) {
  stop("Temporal stimulus not supported by Octopus 600")
}

octo600.opiPresent.opiKineticStimulus <- function(stim, nextStim=NULL, ...) {
  stop("Kinetic stimulus not supported by Octopus 600")
}

###########################################################################
#
# Input paras are the Octopus600Env$* constants
# lum is in cd/m^2 (as per OPI spec) * 100 == .Octopus600Env$BG_{OFF | 1 | 10 | 100 }
# color is .Octopus600Env$MET_COL_{WW | BY | RW | BLUE_WHITE | RED_YELLOW | WHITE_YELLOW }
# fixation is .Octopus600Env$FIX_{RING | CROSS | CENTRE}
# fixIntensity is 0..100 %
#
# @return NULL is succeed.
# @return -1 if opiInitialize has not been successfully called
# @return -2 trouble setting backgound color
# @return -3 trouble setting fixation
###########################################################################
octo600.opiSetBackground <- function(bgColor=NA, fixType=NA, fixColor=NA, fixIntensity=255) {
  
  if (!is.element(fixType, 1:4))
    stop("Fixation type must be 1, 2, 3 or 4 in opiSetBackground()")
  
  if (!is.element(fixColor, c(0, 2, 4, 5)))
    stop("Fixation color must be 0, 2, 4 or 5 in opiSetBackground()")
  
  if (!is.na(bgColor) && !is.element(bgColor, 0:255))
    stop("Background color must be NA or an integer between 0 and 255")
  
  #todo return -1 if opiInitialize has not been successfully called
  
  # setBackground()
  res = sendCommand(.Octopus600Env$socket, 2021, ifelse(is.na(bgColor), 0, bgColor), ifelse(is.na(bgColor), .Octopus600Env$pulsar*5, -1))
  if (res[[2]] != 0)
    return(-2)
  
  # set_fixationmark()
  res = sendCommand(.Octopus600Env$socket, 2003, .Octopus600Env$eye == "left", fixType, fixColor, fixIntensity)
  if (res[[2]] != 0)
    return(-3)
  
  return(NULL)
}

###########################################################################
# return NULL on success (in fact, always!)
###########################################################################
octo600.opiClose <- function() {
    close(.Octopus600Env$socket)
    return(NULL)
}

###########################################################################
# Call opiPresent with a NULL stimulus
###########################################################################
octo600.opiQueryDevice <- function() {
  res <- sendCommand(.Octopus600Env$socket, 3004)
  
  ret <- list(
    answerButton        = res[[1]][8],
    headSensor          = res[[1]][9],
    eyeLidClosureLeft   = res[[1]][10],
    eyeLidClosureRight  = res[[1]][11],
    fixationLostLeft    = res[[1]][12],
    fixationLostRight   = res[[1]][13],
    pupilPositionXLeft  = res[[1]][14],
    pupilPositionYLeft  = res[[1]][15],
    pupilPositionXRight = res[[1]][16],
    pupilPositionYRight = res[[1]][17]
  )
  return(ret)
}
