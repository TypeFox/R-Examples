##' Design Matrix
##' 
##' This function returns the design matrix for an ordinal
##' Bradley-Terry-Luce model.
##' 
##' @usage design(X, var1, var2, use.vars=NULL, reference=NULL, 
##'               prefix="GAMMA", prefix.home="ALPHA", 
##'               home.advantage=c("no","specific","yes"))
##' 
##' @param X a data frame in long format (see \code{\link{wide2long}}).
##' @param var1 a character of the column name from \code{X} 
##' specifying the first object to be compared (or, in sport context, the home team).
##' @param var2 a character of the column name from \code{X} 
##' specifying the second object to be compared (or, in sport context, the away team).
##' @param use.vars a character vector with the names of additional covariates 
##' of \code{X} that will additionally be included into the design matrix. 
##' (example: \code{use.vars = c("ENG", "SEX")}) if the covariates \code{ENG} and \code{SEX} should be included. 
##' If all covariates of \code{X} should be included, you can use \code{use.vars = "ALL"}. 
##' The default is \code{use.vars = NULL} for no additional covariates.
##' @param reference a character specifying the reference object. 
##' @param prefix (optional) a character added in the names of the estimated object parameters
##' @param prefix.home (optional) a character added in the names of the estimated home advantage parameters 
##' @param home.advantage Note that the home advantage is equivalent to an order effect
##'  \code{home.advantage="no"} uses no home advantage (order effect),
##'  \code{home.advantage="specific"} uses one home advantage (order effect) for each object and
##'  \code{home.advantage="yes"} uses one home advantage (order effect) for any object.
##' 
##' @return A data frame where each row refers to a pair comparison and each column corresponds to an object.
##' 
##' @author Giuseppe Casalicchio
##' 
##' @example inst/examples/design_ex.R
##' @export 
##' @importFrom caret dummyVars

design <-
  function(X, var1, var2, use.vars=NULL, reference=NULL,
           prefix="GAMMA", prefix.home="ALPHA",
           home.advantage=c("no","specific","yes")){
    dummy <-
      function(form, data, ...){
        dum <- dummyVars(form, data, ...)
        pred <- predict(dum, data)
        return(pred)
      }
    
    ### input 
    # X: dataset
    # var1: home Team
    # var2: guest Team
    # use.vars: character vector for additional (subject-specific) variables
    #           use "ALL" for all variables.
    # home.advantage: logical if home advantage should be considered
    ### output
    # Transformed dataset where 1: heimmannschaft, -1: gastmannschaft
    home.advantage <- match.arg(home.advantage)
    if(!is.null(use.vars)) {
      if(length(use.vars) == 1 && use.vars=="ALL") 
        append <- subset(X, select=-c(get(var1),get(var2))) else{
          append <- X[, use.vars, drop=FALSE]
        }
    } else append <- NULL
    if(!inherits(X[, var1], "factor") | 
         !inherits(X[, var2], "factor") | 
         length(unique(X[, var1]))!=length(levels(X[, var1])) | 
         length(unique(X[, var2]))!=length(levels(X[, var2]))){
      X[, var1] <- as.factor(as.character(X[, var1]))
      X[, var2] <- as.factor(as.character(X[, var2]))
    }  
    if(length(reference)>1) stop("you can set only one reference object")
    
    missA <- levels(X[,var2])[!levels(X[,var2])%in%levels(X[,var1])]
    missB <- levels(X[,var1])[!levels(X[,var1])%in%levels(X[,var2])]
    
    levels(X[,var1]) <- c(levels(X[,var1]),missA)
    levels(X[,var2]) <- c(levels(X[,var2]),missB)
    
    indA <- dummy(as.formula(paste("~",var1)), X, levelsOnly=TRUE)
    indB <- dummy(as.formula(paste("~",var2)), X, levelsOnly=TRUE)
    
    level <- levels(X[,var1]) #unique(c(as.character(X[,var1]), as.character(X[,var2])))
    ret <- indA[,level]-indB[,level]
    colnames(ret) <- #gsub(" ",".",gsub("|[[:punct:]]", "", colnames(ret)))
      paste(prefix,gsub(" ",".",gsub("|[[:punct:]]", "", colnames(ret))), sep=".")
    
    if(home.advantage=="specific"){
      colnames(indA) <- 
        paste(prefix.home,gsub(" ",".",gsub("|[[:punct:]]", "", colnames(indA))), sep=".")
      if(is.null(use.vars)) append <- indA else append <- cbind(indA,append)
    }
    if(home.advantage=="yes"){
      if(!is.null(use.vars)) append$ALPHA <- 1 else {
        append <- data.frame(rep(1, nrow(ret)))
        colnames(append) <- prefix.home
      }
    }
    if(!is.null(reference)){
      if(grepl(prefix, reference)){
        ret <- ret[,!grepl(reference,colnames(ret))]
      } else{
        ret <- ret[,!grepl(paste(prefix,reference,sep="."),colnames(ret))]
      }
      
    } 
    output <- cbind(ret,append)
    if(inherits(output, "data.frame")) return(output) else return(as.data.frame(output))
  }

##' Reshapes Paired Comparison Data
##' 
##' This function reshapes a data frame that contains pair comparison data in 
##' 'wide' format into a data frame in 'long' format (see 'Details').
##' 
##' @usage wide2long(data, paircomp, names=NULL, ...)
##' 
##' @param data a data frame in 'wide' format.
##' @param paircomp a character vector specifying the columns from \code{data} 
##' that corresponds to a pair comparison.
##' @param names a character vector of the same length as \code{paircomp} specifying the names (separated with a space) of the two objects 
##' that are compared in a pair comparison. 
##' @param ... arguments to be passed to \code{\link[stats]{reshape}}.
##' 
##' @details 
##' In the 'wide' format each row reflects a certain subject/judge and a column contains for example the results of a pair comparison or is a subject-specific covariate.
##' In the 'long' format each row represents a single pair comparison. 
##' 
##' @return The reshaped data frame.
##' 
##' @seealso
##' \code{\link[stats]{reshape}}, 
##' \code{\link{design}}
##' 
##' @author Giuseppe Casalicchio
##' 
##' @example inst/examples/wide2long_ex.R
##' @export 
##' @import wikibooks

wide2long <-
  function(data, paircomp, names=NULL, ...){
    # input:
    ## paircomp: character vector specifying the pair comparison columns from data
    ## names: character vector with object names (separated with space)
    ## ...
    # output:
    ## dataset long format
    if(is.null(names)) names <- paircomp
    if(length(grep(" ",names))!=length(names))
      stop("some 'names' do not contain space character")
    if(!any(paircomp%in%colnames(data)))
      stop("some 'paircomp' do not match with colnames from dataset")
    
    long <- reshape(data, varying=list(paircomp), direction="long", ...)
    playernames <- strsplit(unlist(lapply(names, rep, nrow(data))), " ")
    
    long$object1 <- as.factor(sapply(playernames, function(X){X[1]}))
    long$object2 <- as.factor(sapply(playernames, function(X){X[2]}))
    
    long$time <- NULL #long$id
    return(long)
  }


#' CEMS Data
#' 
#' The Community of European management schools (CEMS) data used in the paper by 
#' Dittrich et al. (1998, 2001), where the responses (the first 15 variables) 
#' are coded as 1 if the first university was preferred, 2 if no decision was 
#' made and 3 if the second university was preferred.
#' 
#' @usage data(CEMSwide)
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 303 observations on the following 23 variables:
#' \itemize{
#'   \item \code{V1} London vs. Paris
#'   \item \code{V2} London vs. Milano
#'   \item \code{V3} London vs. Paris
#'   \item \code{V3} Paris vs. Milano
#'   \item \code{V4} London vs. St.Gallen
#'   \item \code{V5} Paris vs. St.Gallen
#'   \item \code{V6} Milano vs. St.Gallen
#'   \item \code{V7} London vs. Barcelona
#'   \item \code{V8} Paris vs. Barcelona
#'   \item \code{V9} Milano vs. Barcelona
#'   \item \code{V10} St.Gallen vs. Barcelona
#'   \item \code{V11} London vs. Stockholm
#'   \item \code{V12} Paris vs. Stockholm
#'   \item \code{V13} Milano vs. Stockholm
#'   \item \code{V14} St.Gallen vs. Stockholm
#'   \item \code{V15} Barcelona vs. Stockholm
#'   \item \code{STUD} Main discipline of study : 1= commerce, 0= other (economics, business administration, business education)
#'   \item \code{ENG} Knowledge of English : 0= good, 1= poor
#'   \item \code{FRA} Knowledge of French : 0= good, 1= poor
#'   \item \code{SPA} Knowledge of Spanish : 0= good, 1= poor
#'   \item \code{ITA} Knowledge of Italian : 0= good, 1= poor
#'   \item \code{WOR} Full-time employment while studying: 0= no, 1= yes
#'   \item \code{DEG} Intention to take an international degree: 0= no, 1= yes
#'   \item \code{SEX} Gender : 0= female, 1= male
#' }
#' @name CEMSwide
#' @source The Royal Statistical Society Datasets Website (\url{http://www.blackwellpublishing.com/rss/Readmefiles/dittrich.htm})
#' 
#' @references Dittrich R, Hatzinger R and Katzenbeisser W (1998). 
#' "Modelling the effect of subject-specific covariates in paired 
#' comparison studies with an  application to university rankings." 
#' Applied Statistics, *47*(2), pp. 511-525.
#' 
#' @references Dittrich R, Hatzinger R and Katzenbeisser W (2001). 
#' "Corrigendum: Modelling the effect of subject-specific covariates 
#' in paired comparison studies with an application to university rankings." 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 
#' *50*(2), pp. 247-249.
NULL

#' Typewriter Ribbon Data
#' 
#' The typwriter ribbon data used in the paper by Agresti (1992), where 5 
#' different typewriter ribbon brands were compared pairwise by 30 secretaries.
#' 
#' @usage data(ribbon)
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with columns:
#' \itemize{
#'   \item \code{V1} absolute frequency reflecting strong preference for corresponding \code{obj1}
#'   \item \code{V2} absolute frequency reflecting moderate preference for corresponding \code{obj1}
#'   \item \code{V3} absolute frequency reflecting weak preference for corresponding \code{obj1}
#'   \item \code{V4} absolute frequency reflecting no preference
#'   \item \code{V5} absolute frequency reflecting weak preference for corresponding \code{obj2}
#'   \item \code{V6} absolute frequency reflecting moderate preference for corresponding \code{obj2}
#'   \item \code{V7} absolute frequency reflecting strong preference for corresponding \code{obj2}
#'   \item \code{obj1} integer 1, ..., 5 referring to the 5 ribbon brands
#'   \item \code{obj2} integer 1, ..., 5 referring to the 5 ribbon brands
#' }
#' @name ribbon
#' @source Agresti A (1992). "Analysis of ordinal paired comparison data." 
#' Applied Statistics, pp. 287-297. Table 1.
#' 
#' @references  Agresti A (1992). "Analysis of ordinal paired comparison data." 
#' Applied Statistics, pp. 287-297.
NULL
