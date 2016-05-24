#' Estimate Person Parameters for the 4-PL model
#' 
#' Compute Person Parameters for the 1/2/3/4-PL model and choose between five common estimation techniques: ML, WL, MAP, EAP and a robust estimation. All item parameters are treated as fixed.
#'
#' 
#' With this function you can estimate:
#' \itemize{
#' \item \bold{1-PL model} (Rasch model) by submitting: the data matrix, item difficulties and \bold{nothing else}, since the 1-PL model is merely a 4-PL model with: any slope = 1, any lower asymptote = 0 and any upper asymptote = 1!
#' \item \bold{2-PL model} by submitting: the data matrix, item difficulties and slope parameters. Lower and upper asymptotes are automatically set to 0 und 1 respectively.
#' \item \bold{3-PL model} by submitting anything except the upper asymptote parameters
#' \item \bold{4-PL model} ---> submit all parameters ...
#'}
#' 
#' The probability function of the 4-PL model is:
#'\deqn{P(x_{ij} = 1 | \hat \alpha_i, \hat\beta_i, \hat\gamma_i, \hat\delta_i, \theta_j ) = \hat\gamma_i + (\hat\delta_i-\hat\gamma_i) \frac{exp(\hat \alpha_i (\theta_{j} - \hat\beta_{i}))}{\,1 + exp(\hat\alpha_i (\theta_{j} - \hat\beta_{i}))}}
#' 
#' In our case \eqn{\theta} is to be estimated, and the four item parameters are assumed as fixed (usually these are estimates of a former scaling procedure).
#' 
#' The 3-PL model is the same, except that \eqn{\delta_i = 1, \forall i}.
#' 
#' In the 2-PL model \eqn{\delta_i = 1, \gamma_i = 0, \forall i}.
#' 
#' In the 1-PL model \eqn{\delta_i = 1, \gamma_i = 0, \alpha_i = 1, \forall i}.
#' 
#' .
#' 
#' The \bold{robust} estimation method, applies a Huber-type estimator (Schuster & Yuan, 2011), which downweights responses to items which provide little information for the ability estimation. First a residuum is estimated and on this basis, the weight for each observation is computed.
#' 
#' residuum:
#' \deqn{r_i = \alpha_i(\theta - \beta_i)}
#' 
#' weight:
#' 
#' \deqn{w(r_i) =  1 \rightarrow if\, |r_i| \leq H}
#' \deqn{w(r_i) = H/|r| \rightarrow if\, |r_i| > H}
#' 
#'@param respm An integer matrix, which contains the examinees responses. A persons x items matrix is expected.
#'@param thres A numeric vector or a numeric matrix which contains the threshold parameter for each item. If a matrix is submitted, the first row must contain only \bold{zeroes}!
#'@param slopes A numeric vector, which contains the slope parameters for each item - one parameter per item is expected.
#'@param lowerA A numeric vector, which contains the lower asymptote parameters (kind of guessing parameter) for each item.
#'@param upperA numeric vector, which contains the upper asymptote parameters for each item.
#'@param theta_start A vector which contains a starting value for each person. Currently this is necessary to supply, but soon it will be set automatically if nothing is committed.
#'@param mu A numeric vector of location parameters for each person in case of MAP or EAP estimation. If nothing is submitted this is set to 0 for each person for MAP estimation.
#'@param sigma2 A numeric vector of variance parameters for each person in case of MAP or EAP estimation. If nothing is submitted this is set to 1 for each person for MAP estimation.
#'@param type Which maximization should be applied? There are five valid entries possible: "mle", "wle", "map", "eap" and "robust". To choose between the methods, or just to get a deeper understanding the papers mentioned below are quite helpful. The default is \code{"wle"} which is a good choice in many cases. 
#'@param maxsteps The maximum number of steps the NR Algorithm will take. Default = 100.
#'@param exac How accurate are the estimates supposed to be? Default is 0.001.
#'@param H In case \code{type = "robust"} a Huber ability estimate is performed, and \code{H} modulates how fast the downweighting takes place (for more Details read Schuster & Yuan 2011). 
#'@param ctrl more controls:
#'\itemize{
#' \item \code{killdupli}: Should duplicated response pattern be removed for estimation (estimation is faster)? This is especially resonable in case of a large number of examinees and a small number of items.  Use this option with caution (for map and eap), because persons with different \code{mu} and \code{sigma2} will have different ability estimates despite they responded identically. Default value is \code{FALSE}.
#'
#'\item \code{skipcheck}: Default = FALSE. If TRUE data matrix and arguments are not checked - this saves time e.g. when you use this function for simulations.
#'
#'}
#'
#'@template resulttemplate
#'
#' @seealso \link{PPall}, \link{PP_gpcm}, \link{JKpp}, \link{PV}
#'
#' @useDynLib PP
#' @importFrom Rcpp evalCpp
#'
#'
#'@export
#'
#'@author Manuel Reif
#'@references 
#'Baker, Frank B., and Kim, Seock-Ho (2004). Item Response Theory - Parameter Estimation Techniques. CRC-Press.
#'
#' Barton, M. A., & Lord, F. M. (1981). An Upper Asymptote for the Three-Parameter Logistic Item-Response Model.
#'
#' Birnbaum, A. (1968). Some latent trait models and their use in inferring an examinee's ability. In Lord, F.M. & Novick, M.R. (Eds.), Statistical theories of mental test scores. Reading, MA: Addison-Wesley.
#'
#'Magis, D. (2013). A note on the item information function of the four-parameter logistic model. Applied Psychological Measurement, 37(4), 304-315.
#'
#'Samejima, Fumiko (1993). The bias function of the maximum likelihood estimate of ability for the dichotomous response level. Psychometrika,  58, 195-209.
#'
#'Samejima, Fumiko (1993). An approximation of the bias function of the maximum likelihood estimate of a latent variable for the general case where the item responses are discrete. Psychometrika,  58, 119-138.
#'
#'Schuster, C., & Yuan, K. H. (2011). Robust estimation of latent ability in item response models. Journal of Educational and Behavioral Statistics, 36(6), 720-735.
#'
#'Warm, Thomas A. (1989). Weighted Likelihood Estimation Of Ability In Item Response Theory. Psychometrika, 54, 427-450.
#'
#'Yen, Y.-C., Ho, R.-G., Liao, W.-W., Chen, L.-J., & Kuo, C.-C. (2012). An empirical evaluation of the slip correction in the four parameter logistic models with computerized adaptive testing. Applied Psychological Measurement, 36, 75-87.
#'

#'@example ./R/.examples_4pl.R
#'@keywords Person Parameters, 4pl
#'@rdname PP_4pl
#'



PP_4pl <- function(respm, thres, slopes=NULL, lowerA=NULL, upperA=NULL, theta_start=NULL,
                   mu = NULL, sigma2 = NULL, type="wle", maxsteps=40, exac=0.001,H=1,ctrl=list())
{
  
  ### 
  call <- match.call()  
  attr(call, "date") <- date() 
  attr(call,"version") <- packageVersion("PP")
  ###

nitem <- ncol(respm)  
  
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
  
#   if(cont$cdiag) 
#     {
#       cont$killdupli <- FALSE
#       warning("killdupli in 'ctrl' is forced to FALSE!\n")  
#     }
  
  

## --------- threshold matrix
  
## take care of the threshold 'matrix' - for the 4PL model
# a thres-vector is allowed - but for the internal routines
# it had to be reshaped as a matrix

  
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
    
  } else if(is.vector(thres))
      {
        #iimm <- TRUE 
        thres <- rbind(0,thres)
        #allebigger <- FALSE
        maxsc <- apply(thres,2,function(x)(length(x) - sum(is.na(x)))-1)
      }



  
## --------- check user inputs
if(!cont$skipcheck) # to save time e.g. when running simulations
{
  ## --------- check user inputs
  checkINP(respm, thres, slopes, theta_start,lowerA,upperA, type)
  # -----------------------------------------------------  
}

# in case map is chosen and no mu and/or sigma2 is/are submitted.
if( (any(is.null(mu)) | any(is.null(sigma2))))
  {
  if(any(is.null(mu))) 
    {
    mu <- rep(0,nrow(respm))
    if(type %in% c("map","eap")) warning("all mu's are set to 0! \n")
    }
  
  if(any(is.null(sigma2))) 
    {
    sigma2 <- rep(1,nrow(respm))
    if(type %in% c("map","eap")) warning("all sigma2's are set to 1! \n")
    }
  
  }

  ## ----------------------
  
  
  
# ----- CHOOSE MODEL -------------#

if(is.null(lowerA) & is.null(upperA) & is.null(slopes))
{modest <- "1pl"} else if(is.null(lowerA) & is.null(upperA)) 
{modest <- "2pl"} else if(is.null(upperA))
{modest <- "3pl"} else if(is.null(lowerA))
{modest <- "3pl_upperA"} else 
{modest <- "4pl"}


cat("Estimating: ",modest,"model ... \n")
cat("type =",type,"\n")

  
  
  #### conditional controls
  
  if(modest == "4pl")
  {
    if(!all(length(slopes) == c(length(lowerA),ncol(respm),length(upperA),ncol(thres)))) stop("Check length of sumitted vectors!\n")  
  } else if(modest == "3pl")
    {
      if(!all(length(slopes) == c(length(lowerA),ncol(respm),ncol(thres)))) stop("Check length of sumitted vectors!\n")  
    } else if(modest == "2pl")
      {
        if(!all(length(slopes) == c(ncol(respm),ncol(thres)))) stop("Check length of sumitted vectors!\n")  
      } else if(modest == "3pl_upperA")
        {
          if(!all(length(slopes) == c(ncol(respm),ncol(thres),length(upperA)))) stop("Check length of sumitted vectors!\n")  
        }
    
  
  
  # ----- ----- -------------#
  
  
  # add NA and Inf in case of mle estimation and full oder 0 score
  if(type=="mle" | type=="robust")
  {
    resPPx <- ansol(respm,maxsc)  
    respm <- respm[!is.na(resPPx[,2]),]
  }
  
  
  
  ## kill duplicated if desired #################
  
  if(cont$killdupli)
  {
    # this vector contains the indices, to recreate the large matrix in the end of the estimation
    dupvec <- make_dup(respm)
    respm <- respm[dupvec$ndpat,,drop=FALSE]
  }
  
  
  ##############################################
  
  
  
  
  
  # ----- estimation procedure -------------# 
  

# prepare not submitted vectors
  if(modest == "1pl")
    {
      lowerA <- rep(0,nitem)  
      upperA <- rep(1,nitem)
      slopes <- rep(1,nitem)
    } else if(modest == "2pl")
            {
              lowerA <- rep(0,nitem)  
              upperA <- rep(1,nitem)  
            } else if(modest == "3pl")
              {
                upperA <- rep(1,nitem)  
              } else if(modest == "3pl_upperA")
                {
                  lowerA <- rep(0,nitem)    
                }
        
if(type %in% c("mle","wle","map"))
  {
    
  resPP <- NR_4PL(respm,DELTA = thres,ALPHA = slopes, LOWA = lowerA, UPPA = upperA, THETA = theta_start, wm=type,maxsteps,exac,mu,sigma2,H=H)
  
  resPP$resPP[,2] <- sqrt(resPP$resPP[,2])
    
  } else if(type == "robust")
    {
    
    resPP <- NR_4PL(respm,DELTA = thres,ALPHA = slopes, LOWA = lowerA, UPPA = upperA, THETA = theta_start, wm=type,maxsteps,exac,mu,sigma2,H=H)
    
    resPP$resPP[,2] <- sqrt(resPP$resPP[,2])
    
    } else if(type == "eap")
        {
          resPP <- list()
          resPP$resPP <- eap_4pl(respm, thres, slopes, lowerA=lowerA, upperA=upperA,
                         mu = mu, sigma2 = sigma2)
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
             upperA=upperA,theta_start=theta_start,mu=mu,sigma2=sigma2,cont=cont,H=H)


if(cont$killdupli)
  {
  ipar$dupvec <- dupvec 
  }

## ---------------------------------------------
cat("Estimation finished!\n")
rescall <- list(resPP=resPP,call=call,type=type,ipar=ipar)
class(rescall) <- c("fourpl","ppeo")



return(rescall)
  
  

  
  
  
}
