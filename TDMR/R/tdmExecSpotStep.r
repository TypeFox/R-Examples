require(SPOT);

######################################################################################
# tdmExecSpotStep:
#
#' Execute the tuning process with the given \code{spotStep}
#'
#' If \code{spotStep=="auto"}, construct \code{envT} from the general TDMR settings in \code{tdm}.\cr
#' If \code{spotStep=="rep"} or \code{"report"}, load \code{envT} from file \code{tdm$filenameEnvT} (from a prior tuning run).\cr
#' Then call \code{\link{tdmBigLoop}(envT,spotStep)}.
#'
#' If \code{tdm$filenameEnvT} is NULL, set it to  \code{sub(".conf",".RData",tdm$runList[1])}.
#'
#' @param tdm       a list with general settings for TDMR, see \code{\link{tdmDefaultsFill}}
#' @param spotStep  a string, either "auto", "rep" or "report"
#'
#' @return Environment \code{envT} containing (among others) the elements
#'      \item{\code{runList}}{ \code{=tdm$runList}  }
#'      \item{\code{spotList}}{ \code{=tdm$spotList}  }
#'      \item{\code{tdm}}{ \code{=tdm}  }
#'      \item{\code{getBst}}{ accessor function(confFile,nExp,theTuner) into \code{envT$bstGrid}   }
#'      \item{\code{getRes}}{ accessor function(confFile,nExp,theTuner) into \code{envT$resGrid}   }
#'      \item{\code{sCList}}{ list of spotConfig objects, as many as \code{envT$runList} has elements. Each spotConfig object
#'          \code{sCList[[k]]} contains a list \code{opts} as element, which is read from .apd file specified in \code{envT$runList[k]}.  }
#'
#' @seealso   \code{\link{tdmBigLoop}}, \code{\link{tdmEnvTMakeNew}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de})
#' @export
######################################################################################
tdmExecSpotStep <- function(tdm,spotStep) {
    if (spotStep == "auto") {
      #
      # perform a complete tuning + unbiased eval
      #
      
      ## construct an initial environment envT from the given TDMR settings in tdm
      ## (this contains also the fill-in of other defaults for tdm via
      ##      envT$tdm <- tdmDefaultsFill(tdm);
      ## )
      envT <- tdmEnvTMakeNew(tdm);

      ## the call to tdmBigLoop will start the whole TDMR process:
      ## - for each file in tdm$runList a complete DM tuning is started with each tuning
      ##   method tdm$tuneMethod  (if spotStep=="auto")
      ## - the best result from tuning is fed into an unbiased model build and evaluation run
      ## - results are printed and returned in envT$theFinals
      ## - more detailed results are in other elements of environment envT
      ## - two plots:
      ##      a) the progression of the response variable Y and the parameter variables during tuning
      ##      b) the sensitivity plot for each parameter in the vicinity of the best solution found
      envT <- tdmBigLoop(envT,spotStep);
    }
    else        # i.e. spotStep == "rep" or == "report"
    {
      #
      # re-use prior tuning result from tdm$filenameEnvT:
      # do only spot report and unbiased eval on best tuning solution
      #
      if (is.null(tdm$filenameEnvT)) tdm$filenameEnvT=sub(".conf",".RData",tdm$runList[1],fixed=TRUE);
      envT <- tdmEnvTLoad(tdm$filenameEnvT);
      envT <- tdmBigLoop(envT,"rep");
    }
    envT;
}