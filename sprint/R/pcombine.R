##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008-2011 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.org/licenses/>.  #
#                                                                        #
##########################################################################

`pcombine` <- function(...) {
  rflist <- list(...)
  if ( length(rflist) == 0 ) {
    return(NULL)
  }
  rf <- randomForest::combine(...)
  classRF <- rf$type == "classification"
  haveTest <- !any(sapply(rflist, function(x) is.null(x$test)))
  nforest <- length(rflist)
  if(classRF) {
    haveConfusion <- !any(sapply(rflist, function(x) is.null(x$confusion)))
    if(haveConfusion) {
      rf$confusion <- 0
      for(i in 1:nforest) {
        rf$confusion <- rf$confusion + rflist[[i]]$confusion
      }
      rf$confusion <- rf$confusion / nforest
    }
    rf$err.rate <- NULL
    for(i in 1:nforest) {
      rf$err.rate <- rbind(rf$err.rate, rflist[[i]]$err.rate)
    }
    if(haveTest) {
      if(!any(sapply(rflist, function(x) is.null(x$test$confusion)))) {
        rf$test$confusion <- 0
        for(i in 1:nforest) {
          rf$test$confusion <- rf$test$confusion + rflist[[i]]$test$confusion
        }
        rf$test$confusion <- rf$test$confusion / nforest
      }
      rf$test$err.rate <- NULL
      for(i in 1:nforest) {
        rf$test$err.rate <- rbind(rf$test$err.rate, rflist[[i]]$test$err.rate)
      }
    }
  } else {
    rf$mse <- rf$rsq <- NULL
    for ( i in 1:nforest ) {
      rf$mse <- rbind(rf$mse, rflist[[i]]$mse)
      rf$rsq <- rbind(rf$rsq, rflist[[i]]$rsq)
    }
    if(haveTest) {
      rf$test$mse <- rf$test$rsq <- NULL
      for ( i in 1:nforest ) {
        rf$test$mse <- rbind(rf$test$mse, rflist[[i]]$test$mse)
        rf$test$rsq <- rbind(rf$test$rsq, rflist[[i]]$test$rsq)
      }
    }
  }
  rf
}
