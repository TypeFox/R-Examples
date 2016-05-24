#### #
#### # Glue to allow any algorithms to be  algorithm for a single location.
#### #
#### # Author: Andrew Turpin    (aturpin@unimelb.edu.au)
#### #         Tony Redmond
#### # Date: July 2012
#### #
#### # Copyright 2012 Andrew Turpin and Jonathan Denniss
#### # This program is part of the OPI (http://perimetry.org/OPI).
#### # OPI is free software: you can redistribute it and/or modify
#### # it under the terms of the GNU General Public License as published by
#### # the Free Software Foundation, either version 3 of the License, or
#### # any later version.
#### #
#### # This program is distributed in the hope that it will be useful,
#### # but WITHOUT ANY WARRANTY; without even the implied warranty of
#### # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#### # GNU General Public License for more details.
#### #
#### # You should have received a copy of the GNU General Public License
#### # along with this program.  If not, see <http://www.gnu.org/licenses/>.
#### #
#### 
#### ################################################################################
#### # Input parameters
#### #   params    A matrix where each row is 
#### #                   x y luminance-level number-of-presentations
#### #
#### #   interStimMethod Control inter-stimulus interval. Next stim begins
#### #               "none" - as soon as possible after 
#### #                        button press or responseWindow expires
#### #               "speed" - after an average of the last 'speedHistory' 
#### #                         response times, with a minimum of 'speedFloor'
#### #               "constant" - after 'constantWait' ms
#### #               "atLeast"  - after at least 'atLeastWait' ms including reaction time 
#### #
#### #   jitter    Range of random wait time to add to inter-stimulus interval.
#### #
#### #   makeStim  A helper function to create the required
#### #             OPI data type for passing to opiPresent
#### #   ...       Parameters for opiPresent
#### # Returns a list containing
#### #   npres    Total number of presentations
#### #   respSeq  Response sequence stored as a list of (seen,dB) pairs
#### #   first    First staircase estimate in dB
#### #   final    Final threshold estimate in dB
#### ################################################################################
#### MOCS <- function(params=NA, interStimMethod="none", 
####                  speedHistory=5, speedFloor=200, 
####                  constantWait=1000,
####                  atLeastWait=1000,
####                  jitter=150,
####                  makeStim, ...) {
#### 
####         ################################################
####         # expand the matrix to every presentation
####         # and randomise the order of rows in the matrix
####         ################################################
####     mocs <- NULL
####     for(i in 1:nrow(params))
####         mocs <- rbind(mocs, matrix(params[i,1:3], ncol=3, nrow=params[i,4], byrow=T))
####     mocs <- mocs[order(runif(nrow(mocs))), ]
####     mocs <- rbind(mocs, c(NA,NA,NA)) # dummy last presentation for nextStim
#### 
####         ####################################################
####         # loop through every presentation except last
####         ####################################################
####     results <- data.frame()
####     for(i in 1:(nrow(mocs)-1)) {
####         stim     <- makeStim(as.double(mocs[i,]))
####         nextStim <- makeStim(as.double(mocs[i+1,]))
####         class(stim)     <- "opiStaticStimulus"
####         class(nextStim) <- "opiStaticStimulus"
#### 
####         ret <- opiPresent(stim, nextStim, ...)
#### 
####         counter <- 1
####         while ((ret$err == "NotValid") && (counter < 5)) {
####             ret <- opiPresent(stim, nextStim)
####             counter <- counter + 1
####         }
#### 
####         if(counter == 5){
####             plot(0:10, 0:10, axes=FALSE, type="n", xlab="", ylab="")
####             text(5,5,"TEST PAUSED - CLICK MOUSE HERE TO RESUME")
####             locator(1)
####         }
#### 
####         print(paste(stim$x, stim$y, stim$level, ret$err, ret$seen, ret$time))
####         results[i,] <- c(stim$x, stim$y, stim$level, ret$seen, ret$time, ret$err)
####     }
#### }
#### 
#### #    return(list(
#### #        npres=length(fullResponseSeq),  # number of presentations
#### #        respSeq=fullResponseSeq,        # reposnse sequence (list of pairs)
#### #        first=first$final,              # estimate from first staircase
#### #        final=final                     # final threshold estimate
#### #    ))
#### }#MOCS()
#### 
#### t <- matrix(c(
####     9,9, 3145  ,  9,
####     6,6,  314  ,  6,
####     3,3,   31.4,  3
#### ), ncol=4, byrow=TRUE)
#### 
#### makeStim <- function(p) {
####     return(list(
####         x=p[1],
####         y=p[2],
####         level=p[3],
####         size=0.43,
####         duration=400, 
####         responseWindow=1500))
#### }
