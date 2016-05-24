#' @title Item Fit Indices
#' @export pairwise.item.fit
#' @exportClass pairwise_item_fit
#' @description function for calculating item fit indices. The procedures for calculating the fit indices are based on the formulas given in Wright & Masters, (1982, P. 100), with further clarification given in \code{http://www.rasch.org/rmt/rmt34e.htm}.
#' @details contrary to many IRT software using Ml based item parameter estimation, \code{pairwise} will not exclude persons, showing perfect response vectors (e.g. c(0,0,0) for dataset with three variables), prior to the scaling. Therefor the fit statistics computed with \code{pairwise} may deviate somewhat from the fit statistics produced by IRT software using Ml based item parameter estimation (e.g. R-package \code{eRm}), depending on the amount of persons with perfect response vectors in the data.
#' @param pers_obj an object of class \code{"pers"} as a result from function \code{\link{pers}}
#' @param na_treat value to be assigned to residual cells which have missing data in the original response matrix. default is set to \code{na_treat=NA} to ignore these cells in further calculations. An option is to set these residuals to 0 using \code{na_treat=0}, which implys that they are imputed as 'fitting data', i.e., zero residuals. This can attenuate contrasts (see. http://www.rasch.org/rmt/rmt142m.htm).
#' @return an object of class \code{c("pairwise_item_fit", "data.frame")} containing item fit indices.
#' @references Wright, B. D., & Masters, G. N. (1982). \emph{Rating Scale Analysis.} Chicago: MESA Press.
#' @examples ########
#' data(sim200x3)
#' result <- pers(pair(sim200x3))
#' pairwise.item.fit(pers_obj=result) # item fit statistic
####################################################
####################################################


pairwise.item.fit <- function(pers_obj,na_treat=NA){
  # needs internal functions pvx,  pvx.matrix and expscore
  obj <- expscore(pers_obj,na_treat=na_treat) # calls internal function  
  emp_resp <- pers_obj$pair$resp
  Eni <- obj$Eni # expected scores (Expected Mean of ...) gegencheck eRm OK
  Wni <- obj$Wni # Variance of ... gegencheck eRm OK
  Cni <- obj$Cni # Kurtosis of ... gegencheck eRm OK
  Yni <- obj$Yni # score residual ... gegencheck eRm OK
  Zni <- obj$Zni # standardised residual ... gegencheck eRm (st.res in itemfit.ppar) OK
  Y2ni <- obj$Y2ni 
  Z2ni <- obj$Z2ni #standardised residual squared ... gegencheck eRm (sq.res in itemfit.ppar) OK
#-----------------------------------------------------------------
# Nna_v <- colSums(!is.na(Z2ni))
Nna_v <- colSums(!is.na(emp_resp))

Chi <- colSums(Z2ni,na.rm=TRUE) # ... gegencheck eRm (ifit in itemfit.ppar) OK
df <- Nna_v-1 # sowieso besser als eRm, da wird das -1 vergessen  
pChi <- 1-pchisq(Chi, df) # p-value  
#loc.chi.square.p<-1-pchisq(loc.chi.square,loc.df)  
#-----------------------------------------------------------------

## Variance Uq2i of -> Unweighted Mean Square Ui () -------
Uq2i  <-  (colSums( (Cni / (Wni^2)), na.rm = TRUE) / (Nna_v)^2 - (1/Nna_v)   ) # ... gegencheck eRm (qsq.outfitMSQ in itemfit.ppar) OK 
Uqi <- sqrt(Uq2i)

## Unweighted Mean Square Ui (OUTFIT.MEANSQ)-------
# so macht es eRm als alternative (dritte stelle hintem komma versch.):   i.outfitMSQ <- Chi/df
#colSums(!is.na(Z2ni))
Ui <- colSums(Z2ni, na.rm = TRUE)/Nna_v   # nicht N wegen missings!
Uikorr <- Ui+1-mean(Ui) # korr. outfit --> siehe oRM.pdf seite 8 oben und oRM.R Zeile 115 und Ordinales_Rasch_Modell.pdf seite 68

## Standardised (Un)weighted Mean Square Ti (OUTFIT.ZSTD)-------
UTi <- ( ( (Ui^(1/3)) -1) * (3/Uqi) ) + (Uqi/3) # ... gegencheck eRm (i.outfitZ in itemfit.ppar) formel stimmt - werte leicht unterschiedlich - in eRm werden perfecte resp. vorher rausgeworfen OK
UTikorr <- ( ( (Uikorr^(1/3)) -1) * (3/Uqi) ) + (Uqi/3)

#-----------------------------------------------------------------

## Variance Vq2i of -> Weighted Mean Square Vi (INFIT) -------
Vq2i  <- colSums( (Cni - (Wni^2)), na.rm = TRUE) / (colSums(Wni, na.rm = TRUE)^2) # ... gegencheck eRm (qsq.infitMSQ in itemfit.ppar) OK
Vqi <- sqrt(Vq2i)

## Weighted Mean Square Vi (INFIT.MEANSQ)-------
# Vi <- colSums(Z2ni*Wni, na.rm = TRUE)/colSums(Wni, na.rm = TRUE) # eRm style -> identisch
Vi <- colSums(Y2ni, na.rm = TRUE) / colSums(Wni, na.rm = TRUE) # ... gegencheck eRm (i.infitMSQ in itemfit.ppar) OK
Vikorr <- Vi+1-mean(Vi) # korr. outfit --> siehe oRM.pdf seite 8 oben und oRM.R Zeile 115 und Ordinales_Rasch_Modell.pdf seite 68

## Standardised Weighted Mean Square Ti (INFIT.ZSTD)-------
VTi <- ( (Vi^(1/3)-1) * (3/Vqi) ) + (Vqi/3) # ... gegencheck eRm (i.infitZ in itemfit.ppar) unterschiede! aber formel stimmt OK
VTikorr <- ( (Vikorr^(1/3)-1) * (3/Vqi) ) + (Vqi/3)

#-----------------------------------------------------------------
erg <- as.data.frame(list(Chi=round(Chi,4), df=df, p=round(pChi,4), OUTFIT.MSQ=round(Ui,4) , OUTFIT.ZSTD=round(UTi,4) ,INFIT.MSQ=round(Vi,4), INFIT.ZSTD=round(VTi,4) , OUTFIT.MSQ.REL=round(Uikorr,4), OUTFIT.ZSTD.REL=round(UTikorr,4) ,INFIT.MSQ.REL=round(Vikorr,4), INFIT.ZSTD.REL=round(VTikorr,4) ))

class(erg) <- c("pairwise.item.fit","data.frame")
return( erg )
}
