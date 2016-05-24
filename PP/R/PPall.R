#' Estimate Person Parameters
#' 
#' Compute person parameters for the 1,2,3,4-PL model and for the GPCM. Choose between ML, WL, MAP, EAP and robust estimation. Use this function if 4-PL items and GPCM items are mixed for each person.
#'
#' 
#' For a test with both: dichotomous and polytomous items which have been scaled under 1/2/3/4-PL model or the GPCM, use this function to estimate the person ability parameters. You have to define the appropriate model for each item. 
#' 
#'Please note, that  \code{robust} estimation with (Huber ability estimate) polytomous items is still experimental!
#' 
#'@param respm An integer matrix, which contains the examinees responses. A persons x items matrix is expected.
#'@param thres A numeric matrix which contains the threshold parameter for each item. If the first row of the matrix is not set to zero (only zeroes in the first row) - then a row-vector with zeroes is added by default.
#'@param slopes A numeric vector, which contains the slope parameters for each item - one parameter per item is expected. 
#'@param lowerA A numeric vector, which contains the lower asymptote parameters (kind of guessing parameter) for each item. In the case of polytomous items, the value must be 0.
#'@param upperA numeric vector, which contains the upper asymptote parameters for each item. In the case of polytomous items, the value must be 1.
#'@param theta_start A vector which contains a starting value for each person. Currently this is necessary to supply, but soon it will be set automatically if nothing is committed.
#'@param mu A numeric vector of location parameters for each person in case of MAP estimation. If nothing is submitted this is set to 0 for each person for MAP estimation.
#'@param sigma2 A numeric vector of variance parameters for each person in case of MAP or EAP estimation. If nothing is submitted this is set to 1 for each person for MAP estimation.
#'@param type Which maximization should be applied? There are five valid entries possible: "mle", "wle", "map", "eap" and "robust". To choose between the methods, or just to get a deeper understanding the papers mentioned below are quite helpful. The default is \code{"wle"} which is a good choice in many cases.
#'
#'@param model2est A character vector with length equal to the number of submitted items. It defines itemwise the response model under which the item parameter was estimated. There are 2 valid inputs up to now: \code{"GPCM"} and \code{"4PL"}.
#'
#'@param maxsteps The maximum number of steps the NR algorithm will take. Default = 100.
#'@param exac How accurate are the estimates supposed to be? Default is 0.001.
#'@param H In case \code{type = "robust"} a Huber ability estimate is performed, and \code{H} modulates how fast the downweighting takes place (for more Details read Schuster & Yuan 2011).
#'@param ctrl More controls:
#'
#'\itemize{
#' \item \code{killdupli}: Should duplicated response pattern be removed for estimation (estimation is faster)? This is especially resonable in case of a large number of examinees and a small number of items.  Use this option with caution (for map and eap), because persons with different \code{mu} and \code{sigma2} will have different ability estimates despite they responded identically. Default value is \code{FALSE}.
#'
#'\item \code{skipcheck}: Default = FALSE. If TRUE data matrix and arguments are not checked - this saves time e.g. when you use this function for simulations.
#'
#'}
#'
#'@template resulttemplate
#'
#' @seealso \link{PP_gpcm}, \link{PP_4pl}, \link{JKpp}, \link{PV}
#'
#'@export
#'
#'@author Manuel Reif
#'
#'@references Baker, Frank B., and Kim, Seock-Ho (2004). Item Response Theory - Parameter Estimation Techniques. CRC-Press.
#'
#'Barton, M. A., & Lord, F. M. (1981). An Upper Asymptote for the Three-Parameter Logistic Item-Response Model.
#'
#'Magis, D. (2013). A note on the item information function of the four-parameter logistic model. Applied Psychological Measurement, 37(4), 304-315.
#'
#'Muraki, Eiji (1992). A Generalized Partial Credit Model: Application of an EM Algorithm. Applied Psychological Measurement, 16, 159-176.
#'
#'Muraki, Eiji (1993). Information Functions of the Generalized Partial Credit Model. Applied Psychological Measurement, 17, 351-363.
#'
#'Samejima, Fumiko (1993). The bias function of the maximum likelihood estimate of ability for the dichotomous response level. Psychometrika,  58, 195-209.
#'
#'Samejima, Fumiko (1993). An approximation of the bias function of the maximum likelihood estimate of a latent variable for the general case where the item responses are discrete. Psychometrika,  58, 119-138.
#'
#'Schuster, C., & Yuan, K. H. (2011). Robust estimation of latent ability in item response models. Journal of Educational and Behavioral Statistics, 36(6), 720-735.
#'
#'Wang, S. and Wang, T. (2001). Precision of Warm's Weighted Likelihood Estimates for a Polytomous Model in Computerized Adaptive Testing. Applied Psychological Measurement, 25, 317-331.
#'
#'Warm, Thomas A. (1989). Weighted Likelihood Estimation Of Ability In Item Response Theory. Psychometrika, 54, 427-450.
#'
#'Yen, Y.-C., Ho, R.-G., Liao, W.-W., Chen, L.-J., & Kuo, C.-C. (2012). An empirical evaluation of the slip correction in the four parameter logistic models with computerized adaptive testing. Applied Psychological Measurement, 36, 75-87.
#'

#'@example ./R/.examples.R
#'@rdname PPall
#'
PPall <- function(respm, thres, slopes, lowerA, upperA, theta_start=NULL,
                  mu = NULL, sigma2 = NULL, type="wle", model2est,
                  maxsteps=100, exac=0.001,H=1,ctrl=list())
{
### 
call <- match.call()  
attr(call, "date") <- date() 
attr(call,"version") <- packageVersion("PP")
###


## --------- user controls
cont <- list(killdupli=FALSE,skipcheck=FALSE)

user_ctrlI <- match(names(ctrl),names(cont))
if(any(is.na(user_ctrlI)))
  {
    notex <- names(ctrl)[is.na(user_ctrlI)]
    warning("The following options in ctrl do not exist: ", paste(notex,collapse=", "))
    ctrl       <- ctrl[!is.na(user_ctrlI)]
    user_ctrlI <- user_ctrlI[!is.na(user_ctrlI)]
  }

cont[user_ctrlI] <- ctrl



## starting values

if(is.null(theta_start))
  {
    theta_start <- rep(0,nrow(respm))
  }

#---

# if(cont$cdiag) 
#   {
#   cont$killdupli <- FALSE
#   warning("killdupli in 'ctrl' is forced to FALSE!\n")  
#   }



######## take care of the threshold matrix

if(is.matrix(thres))
  {
  #iimm <- nrow(thres) == 1
  
  if(any(thres[1,] != 0))
    {
    thres <- rbind(0,thres)
    }
  
  # compute the maximal score per item
   maxsc <- apply(thres,2,function(x)(length(x) - sum(is.na(x)))-1)
  # are there any items with more than 2 categories?
   #allebigger <- any(maxsc > 1)
  }

## ----------------------


# check for errors and warnings etc  --------------

## --------- check user inputs
if(!cont$skipcheck) # to save time e.g. when running simulations
{
  ## --------- check user inputs
  checkINP(respm, thres, slopes, theta_start, lowerA, upperA, type)
  # -----------------------------------------------------  
}

model2est <- match.arg(model2est,c("GPCM","4PL"),several.ok = TRUE)
if(length(model2est) < ncol(respm)) stop("Check the properties for the submitted model2est vector!\n")



# in case map is chosen and no mu and/or sigma2 is/are submitted.
if( (any(is.null(mu)) | any(is.null(sigma2))))
  {
  if(any(is.null(mu))) 
    {
    mu <- rep(0,nrow(respm))
    if(type == "map") warning("all mu's are set to 0! \n")
    }
  
  if(any(is.null(sigma2))) 
    {
      sigma2 <- rep(1,nrow(respm))
      if(type == "map") warning("all sigma2's are set to 1! \n")
    }
  
  }

## ----------------------



# ----- MLE: Inf and NA ---------------#
if(type=="mle" | type=="robust")
  {
    resPPx <- ansol(respm,maxsc)  
    respm <- respm[!is.na(resPPx[,2]),]
  }



# ----- remove duplicated -------------# 

if(cont$killdupli)
  {
  # this vector contains the indices, to recreate the large matrix in the end of the estimation
  dupvec <- make_dup(respm)
  respm <- respm[dupvec$ndpat,,drop=FALSE]
  }


# ----- estimation procedure -------------# 

    
     if(any(lowerA[model2est == "GPCM"] != 0)) warning("At least for one item with #categories > 1, there is an lower asymptote != 0 submitted. This will be ignored.\n")

     if(any(upperA[model2est == "GPCM"] != 1)) warning("At least for one item with #categories > 1, there is an upper asymptote != 1 submitted. This will be ignored.\n")
       
modelcat <- paste(unique(model2est),collapse=", ")

cat("Estimating: mixed ",modelcat, "... \n")
cat("type =",type,"\n")



if(type %in% c("mle","wle","map"))
  {
    # ----- estimation procedure -------------# 
    
    resPP <-  NR_mixed(respm,DELTA = thres,ALPHA = slopes, LOWA = lowerA, UPPA = upperA, THETA = theta_start,model=model2est, wm=type,maxsteps,exac,mu,sigma2,H)
    
    resPP$resPP[,2] <- sqrt(resPP$resPP[,2])
    
  } else if(type == "robust")
  {
    warning("Robust estimation for GPCM is still very experimental! \n")
    
    resPP <- NR_mixed(respm,DELTA = thres,ALPHA = slopes, LOWA= lowerA, UPPA = upperA, THETA = theta_start,model=model2est, wm=type,maxsteps,exac,mu,sigma2,H=H) 
    
    resPP$resPP[,2] <- sqrt(resPP$resPP[,2])
    
  } else if(type == "eap")
    {
      resPP <- list()
      resPP$resPP <- eap_mixed(respm, thres, slopes, lowerA=lowerA,
                               upperA=upperA, mu = mu, sigma2 = sigma2, model2est=model2est)
      resPP$nsteps <- 0    
    }


  

### result preperation --------------------------

if(cont$killdupli)
  {
   resPP$resPP <- resPP$resPP[dupvec$posvec,]
  }

if(type=="mle" | type=="robust")
  {
   resPPx[!is.na(resPPx[,2]),] <- resPP$resPP
   resPP$resPP <- resPPx
  }

## ---------------------------------------------

colnames(resPP$resPP) <- c("estimate","SE")

ipar <- list(respm=respm,thres=thres,slopes=slopes,lowerA=lowerA,
             upperA=upperA,theta_start=theta_start,mu=mu,sigma2=sigma2,
             cont=cont,model2est=model2est,H=H)


if(cont$killdupli)
  {
    ipar$dupvec <- dupvec 
  }


## ---------------------------------------------
cat("Estimation finished!\n")
rescall <- list(resPP=resPP,call=call,type=type,ipar=ipar)
class(rescall) <- c("gpcm4pl","ppeo")

  


return(rescall)
}













