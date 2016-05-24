sample.caco<-function(object,p.cases=1,caco.ratio=1, verbose=TRUE){
#
# simulate a case-control sampling from a cohort of families (object)
#
# cases and controls are sampled from probands. Their complete families are returned
#
# p.cases    = proportion of cases to be sampled among those available if p.cases>1 sample exactly p.cases
# caco.ratio = number of controls per case
#
if (!inherits(object,"kin.cohort.sample"))
  stop("Please provide a family data.frame created by kc.simul")

probands <- object$rel==0
cancer   <- object$cancer==1

# selection of cases
#
cases <- object$famid[probands & cancer]
ncases <- length(cases)

if (p.cases < 1)
   ncases <- round(ncases*p.cases)
else if (p.cases >1) {
   if (ncases > p.cases) {
     ncases <- p.cases
     warning(paste ("A max of ",ncases, "will be used"), call.=FALSE)
   }
}

cases <- sample(cases, size=ncases)

# selection of controls
#
controls <- object$famid[probands & !cancer]
maxcontrols<-length(controls)
ncontrols <- round(ncases*caco.ratio)
if (ncontrols>maxcontrols){
  ncontrols <- maxcontrols
  warning(paste ("A max of ",maxcontrols, "will be used"), call.=FALSE)
}

# random sample (no matching)
controls <- sample(controls, size=ncontrols)

if(verbose)
 cat(length(cases), "cases and", length(controls), " controls have been selected\n")

object[object$famid %in% c(cases,controls),]

}
