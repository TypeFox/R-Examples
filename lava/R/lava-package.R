
##' Estimation and simulation of latent variable models
##'
##' Framwork for estimating parameters and simulate data from Latent Variable
##' Models.
##'
##' @name lava-package
##' @importFrom graphics plot lines points abline points text layout
##'     par plot.new plot.window title rect locator segments image
##'     mtext box axis polygon matplot contour contour.default
##'     identify
##' @importFrom grDevices xy.coords col2rgb rgb colors rainbow
##'     topo.colors gray.colors palette colorRampPalette heat.colors
##' @importFrom utils stack combn read.csv getTxtProgressBar
##'     setTxtProgressBar txtProgressBar tail modifyList
##'     getFromNamespace packageVersion write.table methods data
##' @importFrom stats density deriv effects lm family simulate vcov
##'     var cov cor coef model.frame model.weights as.formula
##'     model.matrix rnorm rchisq runif rlnorm pnorm qnorm na.omit AIC
##'     terms logLik qt pt update update.formula confint approxfun
##'     pchisq confint.default formula fft uniroot rbinom predict sd
##'     addmargins residuals dnorm quantile qf cov2cor qchisq
##'     get_all_vars p.adjust rpois rgamma printCoefmat rt glm nlminb
##' @importFrom methods new as
##' @aliases lava-package lava
##' @docType package
##' @author Klaus K. Holst Maintainer: <k.k.holst@@biostat.ku.dk>
##' @keywords package
##' @examples
##'
##' lava()
##'
NULL

##' Longitudinal Bone Mineral Density Data
##'
##' Bone Mineral Density Data consisting of 112 girls randomized to receive
##' calcium og placebo. Longitudinal measurements of bone mineral density
##' (g/cm^2) measured approximately every 6th month in 3 years.
##'
##'
##' @name calcium
##' @docType data
##' @format A data.frame containing 560 (incomplete) observations. The 'person'
##' column defines the individual girls of the study with measurements at
##' visiting times 'visit', and age in years 'age' at the time of visit. The
##' bone mineral density variable is 'bmd' (g/cm^2).
##' @source Vonesh & Chinchilli (1997), Table 5.4.1 on page 228.
##' @keywords datasets
NULL

##' Longitudinal Bone Mineral Density Data (Wide format)
##'
##' Bone Mineral Density Data consisting of 112 girls randomized to receive
##' calcium og placebo. Longitudinal measurements of bone mineral density
##' (g/cm^2) measured approximately every 6th month in 3 years.
##' @name bmd
##' @docType data
##' @source Vonesh & Chinchilli (1997), Table 5.4.1 on page 228.
##' @format data.frame
##' @keywords datasets
##' @seealso calcium
NULL

##' Simulated data
##'
##' Simulated data
##' @name brisa
##' @docType data
##' @format data.frame
##' @source Simulated
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name bmidata
##' @docType data
##' @format data.frame
##' @keywords datasets
NULL

##' Hubble data
##'
##' Velocity (v) and distance (D) measures of 36 Type Ia super-novae from the Hubble
##' Space Telescope
##' @name hubble
##' @docType data
##' @format data.frame
##' @source Freedman, W. L., et al. 2001, AstroPhysicalJournal, 553, 47.
##' @keywords datasets
NULL

##' Hubble data
##'
##' @name hubble2
##' @seealso hubble
##' @docType data
##' @format data.frame
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name indoorenv
##' @docType data
##' @format data.frame
##' @source Simulated
##' @keywords datasets
NULL

##' Missing data example
##'
##' Simulated data generated from model
##' \deqn{E(Y_i\mid X) = X, \quad cov(Y_1,Y_2\mid X)=0.5}
##'
##' The list contains four data sets
##' 1) Complete data
##' 2) MCAR
##' 3) MAR
##' 4) MNAR (missing mechanism depends on variable V correlated with Y1,Y2)
##' @examples
##' data(missingdata)
##' e0 <- estimate(lvm(c(y1,y2)~b*x,y1~~y2),missingdata[[1]]) ## No missing
##' e1 <- estimate(lvm(c(y1,y2)~b*x,y1~~y2),missingdata[[2]]) ## CC (MCAR)
##' e2 <- estimate(lvm(c(y1,y2)~b*x,y1~~y2),missingdata[[2]],missing=TRUE) ## MCAR
##' e3 <- estimate(lvm(c(y1,y2)~b*x,y1~~y2),missingdata[[3]]) ## CC (MAR)
##' e4 <- estimate(lvm(c(y1,y2)~b*x,y1~~y2),missingdata[[3]],missing=TRUE) ## MAR
##' @name missingdata
##' @docType data
##' @format list of data.frames
##' @source Simulated
##' @keywords datasets
NULL

##' Example data (nonlinear model)
##'
##' @name nldata
##' @docType data
##' @format data.frame
##' @source Simulated
##' @keywords datasets
NULL

##' Example SEM data (nonlinear)
##'
##' Simulated data
##' @name nsem
##' @docType data
##' @format data.frame
##' @source Simulated
##' @keywords datasets
NULL

##' Example SEM data
##'
##' Simulated data
##' @name semdata
##' @docType data
##' @source Simulated
##' @format data.frame
##' @keywords datasets
NULL

##' Serotonin data
##'
##' This simulated data mimics a PET imaging study where the 5-HT2A
##' receptor and serotonin transporter (SERT) binding potential has
##' been quantified into 8 different regions. The 5-HT2A
##' cortical regions are considered high-binding regions
## 'which are a priori known to yield quite similar and highly correlated
##' measurements.  These measurements can be regarded as proxy measures of
##' the extra-cellular levels of serotonin in the brain
##' \tabular{rll}{
##'         day    \tab numeric \tab Scan day of the year \cr
##'         age    \tab numeric \tab Age at baseline scan \cr
##'         mem    \tab numeric \tab Memory performance score \cr
##'         depr   \tab numeric \tab Depression (mild) status 500 days after baseline \cr
##'         gene1  \tab numeric \tab Gene marker 1 (HTR2A) \cr
##'         gene2  \tab numeric \tab Gene marker 2 (HTTTLPR) \cr
##'         cau \tab numeric \tab SERT binding, Caudate Nucleus \cr
##'         th  \tab numeric \tab SERT binding, Thalamus \cr
##'         put \tab numeric \tab SERT binding, Putamen \cr
##'         mid \tab numeric \tab SERT binding, Midbrain \cr
##'         aci \tab numeric \tab 5-HT2A binding, Anterior cingulate gyrus \cr
##'         pci  \tab numeric \tab 5-HT2A binding, Posterior cingulate gyrus \cr
##'         sfc \tab numeric \tab 5-HT2A binding, Superior frontal cortex \cr
##'         par \tab numeric \tab 5-HT2A binding, Parietal cortex \cr
##' }
##' @name serotonin
##' @docType data
##' @format data.frame
##' @source Simulated
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @seealso serotonin
##' @name serotonin2
##' @docType data
##' @format data.frame
##' @source Simulated
##' @keywords datasets
NULL

##' Twin menarche data
##'
##' Simulated data
##' \tabular{rll}{
##'         id    \tab numeric \tab Twin-pair id \cr
##'         zyg    \tab character \tab Zygosity (MZ or DZ) \cr
##'         twinnum    \tab numeric \tab Twin number (1 or 2) \cr
##'         agemena    \tab numeric \tab Age at menarche (or censoring) \cr
##'         status    \tab logical \tab Censoring status (observed:=T,censored:=F) \cr
##'         bw    \tab numeric  \tab Birth weight \cr
##'         msmoke    \tab numeric \tab Did mother smoke? (yes:=1,no:=0) \cr
##' }
##' @name twindata
##' @docType data
##' @format data.frame
##' @keywords datasets
##' @source Simulated
NULL


##' For internal use
##'
##' @title For internal use
##' @name startvalues
##' @rdname internal
##' @author Klaus K. Holst
##' @keywords utilities
##' @export
##' @aliases
##' startvalues0 startvalues1 startvalues2 startvalues3
##' starter.multigroup
##' addattr modelPar modelVar matrices pars pars.lvm
##' pars.lvmfit pars.glm score.glm procdata.lvmfit modelPar modelVar
##' matrices reorderdata graph2lvm igraph.lvm subgraph finalize
##' index.lvm index.lvmfit index reindex index<-
##' survival survival<- ordinal ordinal<-
##' rmvn dmvn NR logit expit tigol
##' randomslope randomslope<- lisrel variances offdiags describecoef
##' parlabels rsq stdcoef CoefMat CoefMat.multigroupfit deriv updatelvm
##' checkmultigroup profci estimate.MAR missingModel Inverse
##' gaussian_logLik.lvm addhook gethook multigroup Weight fixsome
##' parfix parfix<- merge IV parameter contr parsedesign Specials decomp.specials
##' getoutcome index index<-
NULL
