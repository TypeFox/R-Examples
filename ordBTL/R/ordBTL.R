##' ordinal Bradley-Terry-Luce model (ordBTL)
##' 
##' Fits ordinal regression models to paired comparison data.
##' 
##' @usage ordBTL(formula, data, family=c("cumulative","acat"), 
##'               family.control = list(), restrict=NULL, ...)
##' 
##' @param formula a formula describing the model to be fitted.
##' @param data a data frame containing the design matrix for the model 
##' (See also \code{\link{design}} to generate such an design matrix).
##' @param family a character specifying which ordinal BTL model should be fitted. 
##' Can be either \code{"cumulative"} for the cumulative link model or \code{"acat"} for the adjacent categories model.
##' @param family.control a list with arguments passed to the corresponding \code{family},
##' either \code{\link[VGAM]{cumulative}} for the cumulative link model or \code{\link[VGAM]{acat}} for the adjacent categories model.
##' @param restrict (optional) a character vector specifying the covariates from \code{formula} that should be fitted with a symmetry constraint (can be used to fit threshold covariates).
##' @param ... further arguments for fitting function (see \code{\link[VGAM]{vglm}}).
##'
##' @author Giuseppe Casalicchio
##' 
##' @return An object of class \code{vglm}.
##' 
##' @references Dittrich R, Hatzinger R and Katzenbeisser W (2001). 
##' "Corrigendum: Modelling the effect of subject-specific covariates 
##' in paired comparison studies with an application to university rankings." 
##' Journal of the Royal Statistical Society: Series C (Applied Statistics), 
##' *50*(2), pp. 247-249.
##' 
##' @references Agresti A (1992). "Analysis of ordinal paired comparison data." 
##' Applied Statistics, pp. 287-297. Table 1.
##' 
##' @seealso 
##' \code{\link[VGAM]{vglm}}, 
##' \code{\link[ordBTL]{design}}, 
##' \code{\link[VGAM]{plotvgam}}
##' 
##' @example inst/examples/ordBTL_ex.R
##' @export 
##' @import VGAM
# passing arguments like 'etastart' does not work
# ordBTL <- function (formula, data, family = c("cumulative", "acat"), family.control = list(), 
#                     restrict = NULL, ...) 
# {
#   family <- match.arg(family)
#   mod <- vglm(formula, data, family = do.call(family, args = family.control), 
#               constraints = getConstr(formula = formula, data = data, 
#                                       restrict = restrict), qr.arg = FALSE, ...)
#   return(mod)
# }
#
# old ordBTL function from ordBTL v0.7
# ordBTL <- function(formula, data, family=c("cumulative","acat"), 
#                    family.control = list(), restrict=NULL, ...){
#   
#   family <- match.arg(family)
#   mf <- match.call()
#   mf$constraints <- substitute(getConstr(formula=formula, data=data, restrict=restrict))
#   mf$family <- substitute(do.call(family, args=family.control))
#   mf$family.control <- NULL
#   mf$restrict <- NULL
#   mf$qr.arg <- FALSE
#   
#   if(length(grep("s\\(",attr(terms(formula, data=data),"term.labels")))>0){
#     mf[[1]] <- as.name("vgam")
#     mod <- eval(mf)
#   } else{
#     mf[[1]] <- as.name("vglm")
#     mod <- eval(mf)
#   }
#   return(mod)
# }
#
# same as bevore but improved
# ordBTL <- function(formula, data, family=c("cumulative","acat"), 
#                    family.control = list(), restrict=NULL, ...){
#   
#   family <- match.arg(family)
#   mf <- match.call()
#   mf[[1]] <- as.name("vglm")
# 
#   mf$constraints <- substitute(getConstr(formula=formula, data=data, restrict=restrict))
#   mf$family <- substitute(do.call(family, args=family.control))
#   mf$qr.arg <- FALSE
#   mf$family.control <- mf$restrict <- NULL
#   
#   mod <- eval(mf)
# 
#   mod@call <- match.call()
#   return(mod)
# }

ordBTL <- function(formula, data, family=c("cumulative","acat"), 
                    family.control = list(), restrict=NULL, ...){
  
  getConstr <- function(formula, data, restrict=NULL){
    #   if(!is.null(reference)) {
    #     formula <- get.formula(formula, data, reference)
    #   }
    getVars <-
      function(formula, data){
        attr(terms(formula, data=data),"term.labels")
      }
    response <- as.character(formula)[2]#all.vars(formula)[1]
    if(length(grep("cbind",response))!=0){
      nthresholds <- length(gregexpr(",", response)[[1]])
    } else{
      if(!inherits(data[,response], "ordered")){
        #warning("response variable will be transformed to a ordinal factor---see ordered()")
        data[,response] <- as.ordered(data[,response])
      }
      nthresholds <- length(levels(data[,response]))-1
    }
    nointercept <- length(grep("-1",gsub(" ", "", as.character(formula))))>0
    if(nointercept){
      constrVars <- getVars(formula, data)
      #constr <- vector("list", length(constrVars))
      #names(constr) <- constrVars
      #for(i in 1:length(constrVars)) constr[[i]] <- matrix(rep(1, nthresholds))
      #return(constr)
    } else{
      constrVars <- getVars(formula, data)
      constrVars <- c("(Intercept)", constrVars)
    }
    
    constr <- vector("list", length(constrVars))
    names(constr) <- constrVars
    
    for(i in 1:length(constrVars)) constr[[i]] <- matrix(rep(1, nthresholds))
    
    if(nthresholds%%2==0){
      matr <- diag(nthresholds/2)
      if(!nointercept) constr[["(Intercept)"]] <- rbind(matr, apply(t(matr*-1),2,rev))
      if(!is.null(restrict)){
        restrict <- restrict[restrict%in%constrVars] # added 23.09.2013
        for(i in restrict) constr[[i]] <- rbind(matr, apply(t(matr*-1),2,rev))
      }
    } else{
      if(nthresholds==1){
        constr <- NULL #constr[["(Intercept)"]] <- matrix(1,ncol=1)
      } else{
        cols <- (nthresholds-1)/2
        matr1 <- diag(cols)
        matr2 <- diag(cols)*-1
        if(!nointercept) constr[["(Intercept)"]] <- rbind(matr1, rep(0,times=cols), apply(t(matr2),2,rev))
        if(!is.null(restrict)){
          restrict <- restrict[restrict%in%constrVars] # added 23.09.2013
          for(i in restrict) constr[[i]] <- 
            rbind(matr1, rep(0,times=cols), apply(t(matr2),2,rev))
        }
      }
    }
    return(constr)
  }
  family <- match.arg(family)
  
  mod<- eval.parent(substitute(vglm(formula, data, family=do.call(family, args=family.control),
                              constraints=getConstr(formula=formula, data=data, restrict=restrict), ...))) 
  mod@call <- match.call()
  return(mod)
}

##' Function for internal usage
##' 
##' @keywords internal
##' 

getObjects <- function(des, prefix="GAMMA"){
  names <- colnames(des)
  obj <- names[grep(prefix, names)][-length(grep(prefix,names))]
  return(obj)
}

##' Ranking based on Estimates
##' 
##' Extracts the estimated parameters and sorts them based on their estimated values
##' 
##' @usage getRank(ordBTL, decreasing=TRUE, prefix=NULL, reference=NULL)
##' 
##' @param ordBTL a fitted model returned by \code{ordBTL}.
##' @param decreasing logical. Should the sort be increasing or decreasing? 
##' @param prefix (optional) a character that is included in the names of the parameters; only the parameters are returned that include this character (\code{prefix=NULL} extracts all estimated parameters).
##' @param reference (optional) a character specifying the reference object.
##' 
##' @return matrix containing the parameter estimates.
##' 
##' @seealso
##' \code{\link{ordBTL}}
##' 
##' @author Giuseppe Casalicchio
##' 
##' @example inst/examples/getRank_ex.R
##' @export 

getRank <- function(ordBTL, decreasing=TRUE, prefix=NULL, reference=NULL){
  if(!is.null(prefix)) 
    coefs <- summaryvglm(ordBTL)@coef3[grep(prefix, row.names(summaryvglm(ordBTL)@coef3)),, drop=FALSE] else
      coefs <- summaryvglm(ordBTL)@coef3
  
  if(!is.null(reference)){
    if(any(grepl(reference, rownames(coefs)))) 
      warning("Isn't '", reference, "' already included?")
    coefs <- rbind(coefs, c(0, NA, NA))
    rownames(coefs)[nrow(coefs)] <- reference
  } 
  rank <- coefs[order(coefs[,"Estimate"], decreasing=decreasing),, drop=FALSE]
  return(rank)
}

