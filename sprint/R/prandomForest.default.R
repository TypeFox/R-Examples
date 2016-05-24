##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
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
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################

## Because we end up calling back into R to do the work, using the
## serial randomForest library, we don't need to do any argument
## checking here.  Everything is done for us.
`prandomForest.default` <-
    function(x, y=NULL,  xtest=NULL, ytest=NULL, ntree=500,
             mtry=if (!is.null(y) && !is.factor(y))
             max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
             replace=TRUE, classwt=NULL, cutoff=NULL, strata=NULL,
             sampsize = if (replace) nrow(x) else ceiling(.632*nrow(x)),
             nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
             maxnodes=NULL,
             importance=FALSE, localImp=FALSE, nPerm=1,
             proximity=NULL, oob.prox=proximity,
             norm.votes=TRUE, do.trace=FALSE,
             keep.forest=!is.null(y) && is.null(xtest), corr.bias=FALSE,
             keep.inbag=FALSE, ...) {
    ret <- .External("prandomForest",
                     x = x,
                     y = y,
                     xtest = xtest,
                     ytest = ytest,
                     ntree = ntree,
                     mtry = mtry,
                     replace = replace,
                     classwt = classwt,
                     cutoff = cutoff,
                     strata = strata,
                     sampsize = sampsize,
                     nodesize = nodesize,
                     maxnodes = maxnodes,
                     importance = importance,
                     localImp = localImp,
                     nPerm = nPerm,
                     proximity = proximity,
                     oob.prox = oob.prox,
                     norm.votes = norm.votes,
                     do.trace = do.trace,
                     keep.forest = keep.forest,
                     cor.bias = corr.bias,
                     keep.inbag = keep.inbag)
    if ( is.numeric(ret) ) {
      if ( ret == -1 ) {
        warning(paste("MPI is not initialized, aborted prandomForest call.\n"))
        ret <- NULL
      }
      if ( ret == -2 ) {
        warning(paste("No worker processes exist, aborted prandomForest call.\n"))
        ret <- NULL
      }
    }
    ## Either fewer trees than processes, or MPI not working.
    if ( is.null(ret) ) {
      return(ret);
    }

    ret$call <- match.call()
    return(ret)
  }
