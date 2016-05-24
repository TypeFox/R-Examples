seg.lm.fit.boot<-function(y, XREG, Z, PSI, w, offs, opz, n.boot=10, size.boot=NULL, jt=FALSE,
    nonParam=TRUE, random=FALSE){
#random se TRUE prende valori random quando e' errore: comunque devi modificare qualcosa (magari con it.max)
#     per fare restituire la dev in corrispondenza del punto psi-random
#nonParm. se TRUE implemneta il case resampling. Quello semiparam dipende dal non-errore di
#----------------------------------
#  sum.of.squares<-function(obj.seg){
#      #computes the "correct" SumOfSquares from a segmented" fit
#      b<-obj.seg$obj$coef
#      X<-qr.X(obj.seg$obj$qr) #X<-model.matrix(obj.seg)
#      X<-X[,!is.na(b)]
#      b<-b[!is.na(b)]
#      rev.b<-rev(b)
#      rev.b[1:length(obj.seg$psi)]<-0
#      b<-rev(rev.b)
#      new.fitted<-drop(X%*%b)
#      new.res<- obj.seg$obj$residuals + obj.seg$obj$fitted - new.fitted
#      ss<-sum(new.res^2)
#      ss
#      }
#--------
extract.psi<-function(lista){
#serve per estrarre il miglior psi..
    	dev.values<-lista[[1]]
    	psi.values<-lista[[2]]
    	dev.ok<-min(dev.values)
    	id.dev.ok<-which.min(dev.values)
    	if(is.list(psi.values))  psi.values<-matrix(unlist(psi.values),
    		nrow=length(dev.values), byrow=TRUE)
    	if(!is.matrix(psi.values)) psi.values<-matrix(psi.values)
    	psi.ok<-psi.values[id.dev.ok,]
    	r<-list(SumSquares.no.gap=dev.ok, psi=psi.ok)
    	r
	}
#-------------
      visualBoot<-opz$visualBoot
      opz.boot<-opz
      opz.boot$pow=c(1.1,1.2)
      opz1<-opz
      opz1$it.max <-1
      n<-length(y)
      o0<-try(seg.lm.fit(y, XREG, Z, PSI, w, offs, opz, return.all.sol=FALSE), silent=TRUE)
      rangeZ <- apply(Z, 2, range) #serve sempre
      if(!is.list(o0)) {
          o0<- seg.lm.fit(y, XREG, Z, PSI, w, offs, opz, return.all.sol=TRUE)
          o0<-extract.psi(o0)
          if(!nonParam) {warning("using nonparametric boot");nonParam<-TRUE}
          }
      if(is.list(o0)){
        est.psi00<-est.psi0<-o0$psi
        ss00<-o0$SumSquares.no.gap
        if(!nonParam) fitted.ok<-fitted(o0)
        } else {
          if(!nonParam) stop("the first fit failed and I cannot extract fitted values for the semipar boot")
          if(random) {
            est.psi00<-est.psi0<-apply(rangeZ,2,function(r)runif(1,r[1],r[2]))
            PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
            o0<-try(seg.lm.fit(y, XREG, Z, PSI1, w, offs, opz1), silent=TRUE)
            ss00<-o0$SumSquares.no.gap
          } else {
          est.psi00<-est.psi0<-apply(PSI,2,mean)
          ss00<-opz$dev0
        }
        }

      all.est.psi.boot<-all.selected.psi<-all.est.psi<-matrix(, nrow=n.boot, ncol=length(est.psi0))
      all.ss<-all.selected.ss<-rep(NA, n.boot)
      if(is.null(size.boot)) size.boot<-n

#      na<- ,,apply(...,2,function(x)mean(is.na(x)))

      Z.orig<-Z
      if(visualBoot) cat(0, " ", formatC(opz$dev0, 3, format = "f"),"", "(No breakpoint(s))", "\n")
      count.random<-0
      for(k in seq(n.boot)){
          PSI <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
          if(jt) Z<-apply(Z.orig,2,jitter)
          if(nonParam){
              id<-sample(n, size=size.boot, replace=TRUE)
              o.boot<-try(seg.lm.fit(y[id], XREG[id,,drop=FALSE], Z[id,,drop=FALSE], PSI[id,,drop=FALSE],
                w[id], offs[id], opz.boot), silent=TRUE)
          } else {
              yy<-fitted.ok+sample(residuals(o0),size=n, replace=TRUE)
              o.boot<-try(seg.lm.fit(yy, XREG, Z.orig, PSI, weights, offs, opz.boot), silent=TRUE)
          }
          if(is.list(o.boot)){
            all.est.psi.boot[k,]<-est.psi.boot<-o.boot$psi
            } else {
            est.psi.boot<-apply(rangeZ,2,function(r)runif(1,r[1],r[2]))
            }
            PSI <- matrix(rep(est.psi.boot, rep(nrow(Z), length(est.psi.boot))), ncol = length(est.psi.boot))
            opz$h<-max(opz$h*.9, .2)
            opz$it.max<-opz$it.max+1
            o<-try(seg.lm.fit(y, XREG, Z.orig, PSI, w, offs, opz, return.all.sol=TRUE), silent=TRUE)
            if(!is.list(o) && random){
                est.psi0<-apply(rangeZ,2,function(r)runif(1,r[1],r[2]))
                PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
                o<-try(seg.lm.fit(y, XREG, Z, PSI1, w, offs, opz1), silent=TRUE)
                count.random<-count.random+1
              }
            if(is.list(o)){
              if(!"coefficients"%in%names(o$obj)) o<-extract.psi(o)
              all.est.psi[k,]<-o$psi
              all.ss[k]<-o$SumSquares.no.gap
              if(o$SumSquares.no.gap<=ifelse(is.list(o0), o0$SumSquares.no.gap, 10^12)) o0<-o
              est.psi0<-o0$psi
              all.selected.psi[k,] <- est.psi0
              all.selected.ss[k]<-o0$SumSquares.no.gap #min(c(o$SumSquares.no.gap, o0$SumSquares.no.gap))
              }
            if(visualBoot) {
              flush.console()
              spp <- if (k < 10) "" else NULL
              cat(k, spp, "", formatC(o0$SumSquares.no.gap, 3, format = "f"), "\n")
              }
            } #end n.boot
      all.selected.psi<-rbind(est.psi00,all.selected.psi)
      all.selected.ss<-c(ss00, all.selected.ss)

      ris<-list(all.selected.psi=drop(all.selected.psi),all.selected.ss=all.selected.ss, all.psi=all.est.psi, all.ss=all.ss)

      if(is.null(o0$obj)){
          PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
          o0<-try(seg.lm.fit(y, XREG, Z, PSI1, w, offs, opz1), silent=TRUE)
      }
      if(!is.list(o0)) return(0)
      o0$boot.restart<-ris
      return(o0)
      }