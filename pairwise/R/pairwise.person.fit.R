#' @title Person Fit Indices
#' @export pairwise.person.fit
#' @exportClass pairwise_person_fit
#' @description function for calculating person fit indices. The procedures for calculating the fit indices are based on the formulas given in Wright & Masters, (1982, P. 100), with further clarification given in \code{http://www.rasch.org/rmt/rmt34e.htm}.
#' @details contrary to many IRT software using ML based item parameter estimation, \code{pairwise} will not exclude persons, showing perfect response vectors (e.g. c(0,0,0) for dataset with three variables), prior to scaling. Therefor the fit statistics computed with \code{pairwise} may deviate somewhat from the fit statistics produced by IRT software using ML based item parameter estimation (e.g. R-package \code{eRm}), depending on the amount of persons with perfect response vectors in the data.
#' @param pers_obj an object of class \code{"pers"} as a result from function \code{\link{pers}}.
#' @param na_treat value to be assigned to residual cells which have missing data in the original response matrix. default is set to \code{na_treat=NA} to ignore these cells in further calculations. An option is to set these residuals to 0 using \code{na_treat=0}, which implys that they are imputed as 'fitting data', i.e., zero residuals. This can attenuate contrasts (see. http://www.rasch.org/rmt/rmt142m.htm).
#' @return an object of class \code{c("pairwise_person_fit", "data.frame")} containing person fit indices
#' @examples ########
#' data(sim200x3)
#' result <- pers(pair(sim200x3))
#' pairwise.person.fit(pers_obj=result) # item fit statistic
####################################################
####################################################


pairwise.person.fit <- function(pers_obj,na_treat=NA){
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
P_Nna_v <- rowSums(!is.na(emp_resp))

P_Chi <- rowSums(Z2ni,na.rm=TRUE) # ... gegencheck eRm (pfit in personfit.ppar) OK
P_df <- P_Nna_v-1 # OK  
P_pChi <- 1-pchisq(P_Chi, P_df) # p-value   
#-----------------------------------------------------------------
## Variance Uq2i of -> Unweighted Mean Square Ui () -------
P_Uq2i  <-  (rowSums( (Cni / (Wni^2)), na.rm = TRUE) / (P_Nna_v)^2 - (1/P_Nna_v)   ) # ... gegencheck eRm (qsq.outfitMSQ in personfit.ppar) ~ OK 
P_Uqi <- sqrt(P_Uq2i)

## Unweighted Mean Square Ui (OUTFIT.MEANSQ)-------
# so macht es eRm als alternative (dritte stelle hintem komma versch.):   i.outfitMSQ <- Chi/df
P_Ui <- rowSums(Z2ni, na.rm = TRUE)/P_Nna_v   # nicht m wegen missings!
P_Uikorr <- P_Ui+1-mean(P_Ui, na.rm = TRUE) # analog zu korr. outfit --> siehe oRM.pdf seite 8 oben und oRM.R Zeile 115 und Ordinales_Rasch_Modell.pdf seite 68
## Standardised (Un)weighted Mean Square Ti (OUTFIT.ZSTD)-------
P_UTi <- ( ( (P_Ui^(1/3)) -1) * (3/P_Uqi) ) + (P_Uqi/3) # ... gegencheck eRm (p.outfitZ in itemfit.ppar) formel stimmt - werte leicht unterschiedlich - in eRm werden perfecte resp. vorher rausgeworfen OK
P_UTikorr <- ( ( (P_Uikorr^(1/3)) -1) * (3/P_Uqi) ) + (P_Uqi/3) 

#-----------------------------------------------------------------

## Variance Vq2i of -> Weighted Mean Square Vi (INFIT) -------
P_Vq2i  <- rowSums( (Cni - (Wni^2)), na.rm = TRUE) / (rowSums(Wni, na.rm = TRUE)^2) # ... gegencheck eRm (qsq.infitMSQ in itemfit.ppar) OK
P_Vqi <- sqrt(P_Vq2i)

## Weighted Mean Square Vi (INFIT.MEANSQ)-------
P_Vi <- rowSums(Y2ni, na.rm = TRUE) / rowSums(Wni, na.rm = TRUE) # ... gegencheck eRm (p.infitMSQ in personfit.ppar) OK
P_Vikorr <- P_Vi+1-mean(P_Vi, na.rm = TRUE) # analog zu korr. outfit --> siehe oRM.pdf seite 8 oben und oRM.R Zeile 115 und Ordinales_Rasch_Modell.pdf seite 68

## Standardised Weighted Mean Square Ti (INFIT.ZSTD)-------
P_VTi <- ( (P_Vi^(1/3)-1) * (3/P_Vqi) ) + (P_Vqi/3) # ... gegencheck eRm (p.infitZ in personfit.ppar) unterschiede! aber formel stimmt OK
P_VTikorr <- ( (P_Vikorr^(1/3)-1) * (3/P_Vqi) ) + (P_Vqi/3)

#-----------------------------------------------------------------
erg <- as.data.frame(list(Chi=round(P_Chi,4), df=P_df, p=round(P_pChi,4), OUTFIT.MSQ=round(P_Ui,4) , OUTFIT.ZSTD=round(P_UTi,4) ,INFIT.MSQ=round(P_Vi,4), INFIT.ZSTD=round(P_VTi,4), OUTFIT.MSQ.REL=round(P_Uikorr,4), OUTFIT.ZSTD.REL=round(P_UTikorr,4), INFIT.MSQ.REL=round(P_Vikorr,4), INFIT.ZSTD.REL=round(P_VTikorr,4)    ))

class(erg) <- c("pairwise.person.fit","data.frame")
return( erg )
}
