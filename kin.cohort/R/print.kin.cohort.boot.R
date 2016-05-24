`print.kin.cohort.boot` <-
function(x, cumrisk = TRUE, hazard = FALSE, HR=TRUE, conf = 0.95, digits=5, show=TRUE, logrank=TRUE, ...){
  if (!inherits(x,"kin.cohort.boot")) stop("Object should be class kin.cohort.boot")

  boot.ci<-function (t, conf )
  {
    # percentile bootstrap confidence intervals
    # function borrowed from R package boot: 
    # S original by Angelo Canty <cantya@mcmaster.ca>. R port by Brian Ripley <ripley@stats.ox.ac.uk>
    #
    alpha <- (1 + c(-conf, conf))/2
    t <- t[is.finite(t)]
    r <- length(t)
    if(r>0){
      rk <- (r + 1) * alpha
      k <- trunc(rk)
      inds <- 1:length(k)
      out <- inds
      kvs <- k[k > 0 & k < r]
      tstar <- sort(t, partial = sort(union(c(1, r), c(kvs, kvs + 1))))
      ints <- (k == rk)
      if (any(ints))
          out[inds[ints]] <- tstar[k[inds[ints]]]
      out[k == 0] <- tstar[1]
      out[k == r] <- tstar[r]
      not <- function(v) xor(rep(TRUE, length(v)), v)
      temp <- inds[not(ints) & k != 0 & k != r]
      temp1 <- qnorm(alpha[temp])
      temp2 <- qnorm(k[temp]/(r + 1))
      temp3 <- qnorm((k[temp] + 1)/(r + 1))
      tk <- tstar[k[temp]]
      tk1 <- tstar[k[temp] + 1]
      out[temp] <- tk + (temp1 - temp2)/(temp3 - temp2) * (tk1 - tk)
    } else out<-c(NA,NA)
    out
  }

e<- x$estimate
k<- length(e$knots)
if (show){
   cat("Kin-cohort analysis\n")
   if(inherits(x,"wacholder"))  cat("\nWacholder Moments method\n")
   if(inherits(x,"chatterjee")) cat("\nMarginal Likelihood method\n")
}

ncateg<-ifelse(is.list(x$cumrisk),2,4)

#
# cumulative risk
#

if (ncateg==2){

  g0.lab<-colnames(e$cumrisk)[1:2]

  pco<-x$cumrisk$Noncarrier
  pca<-x$cumrisk$Carrier
  pcr<-x$cumrisk$"C.R.R."

  lim.pco<-t(sapply(1:k,function(x)boot.ci(pco[,x], conf)))
  lim.pca<-t(sapply(1:k,function(x)boot.ci(pca[,x], conf)))
  lim.pcr<-t(sapply(1:k,function(x)boot.ci(pcr[,x], conf)))
  
  
  med.pco<-apply(pco, 2, median)
  med.pca<-apply(pca, 2, median)
  med.pcr<-apply(pcr, 2, median)
  
  o.cr<-cbind(e$cumrisk[,1],med.pco,lim.pco,
              e$cumrisk[,2],med.pca,lim.pca,
              e$cumrisk[,3],med.pcr,lim.pcr)
  
  colnames(o.cr)<-c(g0.lab[1],"med boot","(95%","C.I.)",
                    g0.lab[2],"med boot","(95%","C.I.)",
                    "Cum.Risk Ratio","med boot","(95%","C.I.)")
  rownames(o.cr)<-rownames(e$cumrisk)
if(cumrisk & show){
  cat("\nCumulative risk\n")
  print(round(o.cr,digits) )
}
out<-list(cumrisk=o.cr, knots=e$knots, ncateg=ncateg)

#
# hazard
#
if(inherits(x,"chatterjee")){
  
	hco<-x$hazard$Noncarrier
	hca<-x$hazard$Carrier
	hr<-x$hazard$"H.R."
	
	lim.hco<-t(sapply(1:k,function(x)boot.ci(hco[,x], conf)))
	lim.hca<-t(sapply(1:k,function(x)boot.ci(hca[,x], conf)))
	lim.hr <-t(sapply(1:k,function(x)boot.ci(hr[,x],  conf)))
	
	med.hco<-apply(hco, 2, median)
	med.hca<-apply(hca, 2, median)
	med.hr<- apply(hr,  2, median)
	
	o.hz<-cbind(e$hazard[,1],med.hco,lim.hco,
	           e$hazard[,2],med.hca,lim.hca,
	           e$hazard[,3],med.hr, lim.hr)
	
	
	colnames(o.hz)<-c(g0.lab[1],"med boot","(95%","C.I.)",
	                 g0.lab[2],"med boot","(95%","C.I.)",
	                 "Hazard Ratio","med boot","(95%","C.I.)")
	rownames(o.hz)<-rownames(e$hazard)

	if(hazard & show){
	  cat("\nHazard\n")
	  print(round(o.hz,digits) )
	}

   out<-list(cumrisk=o.cr, hazard=o.hz, knots=e$knots, ncateg=ncateg)
}

} else { # 2 genes 4 categories
########################################################################
#
  cr<-x$cumrisk
  dcr<-dim(cr) # ages x groups x B
  crb<-array(NA, c(dcr[1],4, dcr[2])) # last dimension= median, lo, hi

  for (i in 1:dcr[2]) # groups (7)
     crb[,,i]<-t(sapply(1:dcr[1],function(x)c(e$cumrisk[x,i],median(cr[x,i,]),boot.ci(cr[x,i,], conf))))
  dimnames(crb)<-list(ages=rownames(x$estimate$cumrisk), c("Estimate","med boot","(95%","C.I.)"), group=colnames(x$estimate$cumrisk) )

if(cumrisk & show){
  cat("\nCumulative risk\n")
  print(round(crb,digits) )
}
out<-list(cumrisk=crb, knots=e$knots, ncateg=ncateg)

#
# hazard
#
if(inherits(x,"chatterjee")){

  hz<-x$hazard
  dhz<-dim(hz)
  hzb<-array(NA, c(dhz[1],4, dhz[2])) # last dimension= median, lo, hi

  for (i in 1:dhz[2]) # groups (7)
     hzb[,,i]<-t(sapply(1:dhz[1],function(x)c(e$cumrisk[x,i],median(hz[x,i,]),boot.ci(hz[x,i,], conf))))
  dimnames(hzb)<-list(ages=rownames(x$estimate$hazard), c("Estimate","med boot","(95%","C.I.)"), group=colnames(x$estimate$hazard) )

	if(hazard & show){
	  cat("\nHazard\n")
	  print(round(hzb,digits) )
	}
   out<-list(cumrisk=crb, hazard=hzb, knots=e$knots, ncateg=ncateg)
}
}

# HR

if(HR & inherits(x,"chatterjee")){
  lim.logHR<-t(sapply(1:ncol(x$logHR),function(i)boot.ci(x$logHR[,i], conf)))
  med.logHR<-apply(x$logHR, 2, median)

  logHR<-cbind(e$logHR, med.logHR, lim.logHR)
  colnames(logHR)<-c("HR","med boot","(95%","C.I.)")
  rownames(logHR)<-names(e$logHR)
  if(show){
	  cat("\nAverage Hazard Ratio\n")
	  print(round(exp(logHR),digits))
  }
}

if(logrank & show & !is.null(e$logrank)){
   cat("\nLogrank test p-value\n")
   print(signif(e$logrank, digits))
}

class(out)<-c(class(x), "print.kin.cohort.boot")
invisible(out)
}

