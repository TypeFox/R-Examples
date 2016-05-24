gkcov<-function(x,y,gk.sigmamu=taulc,...){
#
# Compute robust covariance using the Gnanadesikan-Kettenring
# estimator.
# (cf. Marrona & Zomar, 2002, Technometrics
#
val<-.25*(gk.sigmamu(x+y,...)-gk.sigmamu(x-y,...))
val
}