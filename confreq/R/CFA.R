#' @title Configural Frequencies Analysis Main Function 
#' @export CFA
#' @exportClass CFA
#' @description Calculates various coefficients for the Configural Frequencies Analysis (CFA) defining main- and (optionaly) interaction effects. The core principle is to use \code{\link{glm}} in package \code{stats} to calculate the expected counts considering a designmatrix, which is constructed based on an formular definition given in argument \code{form}. 
#' @details This is the main function of the package. It internaly calls several functions of the package \code{\link{confreq}} which are also available as single functions. For clasification of the observed patterns into 'Types' and 'Antitypes' according to Linert  (1971), a S3 summary method for the resulting object of class \code{"CFA"} can be applied - see \code{\link{summary.CFA}}.
#' 
#' @param patternfreq an object of class \code{"Pfreq"}, which is data in pattern frequencies representation - see function \code{\link{dat2fre}}.
#' 
#' @param alpha a numeric giving the alpha level for testing (default set to \code{alpha=.05})
#'    
#' @param form either a character expression which can be coerced into a model formula with the function \code{as.formula} in the package \code{stats}. If this argument is left empty (at default \code{form=NULL}) the (internal) function \code{design_cfg_cfa()} will return a designmatrix coding only main effects and no interactions -- for a designmatrix refering to  three variables (V1, V2, V3) for example, leaving the argument \code{form} empty will be equivalent to assigning the character \code{"~ V1 + V2 + V3"} to the argument (\code{form="~ V1 + V2 + V3"}).
#'  A special Case is to define a null-model or rather a cfa model of order zero. In such a model no (main) effects are considered. This can be achieved bei passing the character expression \code{"null"} to the argument \code{form} -- so: \code{form = "null"} -- not to be confound with the default setting of this argument \code{form=NULL}. Another option is to define your own designmatrix and assign it to this argument (\code{form}) in this case the object assigned to \code{form} musst be of class \code{"matrix"} and must logicaly match to the argument \code{patternfreq}, which is currently not checked! - but simply assumed.   
#' 
#' @param ccor either a logical (TRUE / FALSE) determining wether to apply a continuity correction or not. When set to \code{ccor=TRUE} continuity correction is applied for expected values 5 =< expected =< 10. For \code{ccor=FALSE} no continuity correction is applied. Another option is to set \code{ccor=c(x,y)} where x is the lower and y the upper bound for expected values where continuity correction is applied. So \code{ccor=c(5,10)} is equivalent to \code{ccor=TRUE}.
#' 
#' @param family argument passed to \code{\link{glm.fit}} with default set to \code{poisson()}
#' 
#' @param intercept argument passed to \code{\link{glm.fit}} with default set to \code{FALSE} 
#' 
#' @param method charcter defining the estimation method for expected frequencies with default set to \code{method="log"} to estimate the expected frequencies using \code{\link{glm}}. An other option is to set this argument to \code{method="margins"} which will result in expected frequencies calculated based on the margins of the multidimensional contigency table. Only main effects models are posible in this case and thus the arguments \code{form}, \code{family} and \code{intercept} are ignored.
#' 
#' @param blank unsed only if argument \code{method} is set to \code{method="margins"} otherwise ignored. Should be either (1) character vector defining the pattern (with spaces between variable categories), which will be ignored for calculation of expected frequencies; or (2) a numeric vector defining the position(s) of the pattern in object of class \code{"Pfreq"} (see. argument \code{patternfreq}), which will be ignored for calculation of expected frequencies. At default (\code{blank=NULL}) all possible pattern, as listed in object of class \code{"Pfreq"}, are included for calculation of expected frequencies.  
#'  
#' @param ... additional parameters passed through to other functions.
#' @return an object of class \code{CFA} with results.
#' @references Lienert, G. A. (1971). Die Konfigurationsfrequenzanalyse: I. Ein neuer Weg zu Typen und Syndromen. \emph{Zeitschrift für Klinische Psychologie und Psychotherapie, 19}(2), 99-115. 
#' 
#' @examples #######################################
#' ######### some examples ########
#' data(LienertLSD)
#' LienertLSD
#' res1 <- CFA(LienertLSD)
#' summary(res1)
#' ## testing with (full) interactions
#' res2 <- CFA(LienertLSD,form="~ C + T + A + C:T + C:A + T:A + C:T:A")
#' summary(res2)
#' #' ## testing the null model
#' res3 <- CFA(LienertLSD,form="null")
#' summary(res3)
#' #######################
#' data(suicide)
#' suicide
#' # suicide data is in non tabulated data representation - so it must be tabulated !
#' res4 <- CFA(dat2fre(suicide))  
#' summary(res4)

############### start of function definition ##################
CFA<-function(patternfreq, alpha=.05, form=NULL, ccor=FALSE, family=poisson(), intercept=FALSE, method="log", blank=NULL,...){
  
if(any(class(patternfreq)=="Pfreq") != TRUE){stop("patternfreq must be an object of class 'Pfreq'","\n","see func. dat2fre()", call. = TRUE) }

kategorie <- sapply(lapply(patternfreq[,1:(dim(patternfreq)[2]-1) ],levels),length)
  
pattern <- do.call(paste, patternfreq[,1:(dim(patternfreq)[2]-1) ])
  
observed <- patternfreq[,dim(patternfreq)[2]] 

# condition added 22-06-2015
if(method=="log"){
  if(class(form)!="matrix"){
    if(length(form)==0){form<-paste("~", paste(names(patternfreq)[1:(length(patternfreq)-1)],collapse=" + "))}
    
    designmatrix <- design_cfg_cfa(kat=kategorie , form = form, ...) 
    
    usedform <- form
  }
  
  if(class(form)=="matrix"){designmatrix <- form ; usedform <- "designmatrix" } # !!! no further checks !!!
  
  # expected <- expected_cfa(des=designmatrix, observed=observed, family=family, intercept=intercept, ...) # padded out 24.10.2014
  glmfitres <- glm.fit(x=designmatrix, y=observed ,family=family, intercept = intercept) #, ... added 24.10.2014
  expected <- glmfitres$fitted.value # added 24.10.2014
  #aic <- glmfitres$aic # added 24.10.2014
  class(glmfitres) <- c("glm", "lm" )# this is a trick!!! necessary for the following three lines of code - added 24.10.2014
  loglik <- logLik(glmfitres)# cheked against code below OK! added 24.10.2014
  aic <- AIC(glmfitres)# cheked against code below OK! added 24.10.2014
  bic <- BIC(glmfitres)# cheked against code below OK! added 24.10.2014
  # glmres <- glm(formula = paste("Freq",form) , data=patternfreq, , family = family) # added 24.10.2014
  # expected <- fitted(glmres) # OK
  # # resid(glmres)
  # # predict(glmres)
  # loglik <- logLik(glmres)
  # aic <- AIC(glmres)
  # bic <- BIC(glmres)
}

if(method=="margins"){
  expected <- expected_margin_cfa(Pfreq = patternfreq, blank = blank) # added 22-06-2015
  loglik <- NA # added 22-06-2015
  aic <- NA # added 22-06-2015
  bic <- NA # added 22-06-2015
  # nur für df
  form<-paste("~", paste(names(patternfreq)[1:(length(patternfreq)-1)],collapse=" + "))
  designmatrix <- design_cfg_cfa(kat=kategorie , form = form, ...) 
  df <- df_des_cfa(designmatrix)
  designmatrix <- NA # added 22-06-2015
  usedform <- "margins" # added 22-06-2015
}

erg <- data.frame(pat.=pattern, obs.=observed,exp.=expected,do.call(cbind,chi_local_test_cfa(observed,expected)),ex.bin.test=binomial_test_cfa(observed,expected), z_tests_cfa(observed,expected,ccor=ccor),p.stir=stirling_cfa(observed=observed, expected=expected,cum=TRUE,verb=FALSE),density.stir=stirling_cfa(observed=observed, expected=expected,cum=FALSE,verb=FALSE))

chi.square <- sum(erg$Chi)
  
if(method=="log"){df <- df_des_cfa(designmatrix)} ## added 22-06-2015
  
chi.square.p <- (1-pchisq(chi.square,df))
  
lr.chi <- lr(observed,expected) ## added 20. October 2014 JHH
  
lr.p <- (1-pchisq(lr.chi,df)) ## added 20. October 2014 JHH
  
bonferroni <- alpha/length(expected)
  
result <- list( local.test = erg, bonferroni.alpha=bonferroni, global.test = list(pearson = list(Chi=chi.square,df=df,pChi=chi.square.p,alpha=alpha), likelihood.ratio = list(Chi=lr.chi,df=df,pChi=lr.p,alpha=alpha), infocrit=list(loglik=loglik, AIC=aic, BIC=bic)),designmatrix=designmatrix, variables=kategorie, used.formula=usedform) 
  
class(result)<-c("CFA","list")
 
return(result)
}
# End of function ------ 