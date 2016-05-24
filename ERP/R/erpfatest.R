erpfatest <-
function(dta,design,design0=NULL,method="BH",nbf=NULL,nbfmax=15,alpha=0.05,pi0=1,
                                  wantplot=FALSE,s0=NULL,min.err=1e-03,maxiter=5,verbose=FALSE) {  
                                  
   ifa = function(Psi, B) {
     if (class(B) == "numeric") 
         B = matrix(B, ncol = 1)
     q = ncol(B)
     Phi = rep(0, length(Psi))
     Phi[abs(Psi) > 1e-05] = 1/Psi[abs(Psi) > 1e-05]
     PhiB = tcrossprod(Phi, rep(1, q))
     PhiB = PhiB * B
     G = diag(q) + t(B) %*% PhiB
     GinvtPhiB = tcrossprod(solve(G), PhiB)
     Phib2 = tcrossprod(PhiB, t(GinvtPhiB))
     iS = diag(Phi) - Phib2
     PhiB2 = crossprod(PhiB, B)
     GinvtPhiB2 = crossprod(solve(G), PhiB2)
     Phib2 = tcrossprod(PhiB, t(GinvtPhiB2))
     iSB = PhiB - Phib2
     return(list(iS = iS, iSB = iSB))
   }

   emfa = function (data, nbf, min.err = 1e-06,verbose=FALSE) {

      m = ncol(data)
      n = nrow(data)
      mdta = t(rep(1,n))%*%data/n
      vdta = (t(rep(1,n))%*%data^2/n)-mdta^2
      sddta = sqrt(n/(n-1))*sqrt(vdta)
      cdta = data-rep(1,n)%*%mdta
      if (nbf == 0) {
         B = NULL
         Psi = rep(1, m)
         Factors = NULL
      }
      if (nbf > 0) {
         svddta = svd(cdta/sqrt(n-1),nv=n)
         evalues = (svddta$d[1:nbf])^2
         evectors = svddta$v[,1:nbf]

         if (nbf > 1) 
            B = evectors[, 1:nbf] %*% diag(sqrt(evalues[1:nbf]))
         if (nbf == 1) 
            B = matrix(evectors, nrow = m, ncol = 1) * sqrt(evalues[1])
         Psi = sddta^2 - (B^2%*%rep(1,nbf))[,1]
         crit = 1
         while (crit > min.err) {
            iS = ifa(Psi,B)
            xiSB = cdta%*%iS$iSB
            Cyz = t(cdta)%*%xiSB/(n-1)
            Czz = t(iS$iSB)%*%Cyz+diag(nbf)-t(B)%*%iS$iSB
            Bnew = Cyz%*%solve(Czz)
            Psinew = sddta^2 - (Bnew^2%*%rep(1,nbf))[,1]
            crit = mean((Psi - Psinew)^2)
            B = Bnew
            Psi = Psinew
         }
         sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
         G = solve(diag(nbf) + sB %*% t(sB))
         sB = scale(t(B), center = FALSE, scale = Psi)
         Factors = cdta%*%t(sB)%*%t(G)
      }
      res = list(B = B, Psi = Psi,Factors=Factors)
      return(res)
   }

   bivprob = function (rho,lower,upper=-lower,mean=0) {
      nu = 0
      low = rep(as.double((lower - mean)),2)
      upp = rep(as.double((upper - mean)),2)
      if (any(lower == upper)) return(0)
      infin = c(2, 2)
      infin = as.integer(infin)
      low = replace(low, low == -Inf, 0)
      upp = replace(upp, upp == Inf, 0)
      rho = as.double(rho)
      prob = as.double(0)
      a = lapply(rho,function(r,low,upp) biv.nt.prob(df=Inf,lower=low,upper=upp,mean=rep(0,2),S=matrix(c(1,r,r,1),2,2)),
         low=low,upp=upp)
      return(unlist(a))
   }

   Dt = function(rho) {
      threshold=0.05
      ut = qnorm(1 - threshold/2)
      delta = unlist(lapply(rho,bivprob,lower=-ut)) - (1 - threshold)^2
      dt = delta/(threshold * (1 - threshold))
      return(dt)
   }

   VarInflation = function(dta,Blist,maxnbfactors,dig,verbose) {
      m = ncol(dta)
      n = nrow(dta)
      vecrho = round(seq(10^(-dig),1,10^(-dig)),digits=dig)
      vecdt = unlist(lapply(vecrho,Dt))
      sampled = sample(1:m,min(200,m))
      sampsize = length(sampled)
      cordata = t(dta[,sampled])%*%dta[,sampled]/(n-1)
      sdt = rep(0,maxnbfactors+1)
      names(sdt) = paste(0:maxnbfactors,"factors")
      for (i in 1:(maxnbfactors+1)) {
         if (verbose) print(paste("Calculating Variance Inflation criterion for the model with",i-1,"factors"))
         B = matrix(Blist[[i]][sampled,],nrow=sampsize)
         sdb = sqrt(1-apply(B^2,1,sum))
         matrho = cordata - B%*%t(B)
         matrho = sweep(matrho,2,FUN="/",STATS=sdb)
         matrho = sweep(matrho,1,FUN="/",STATS=sdb)
         rho = matrho[col(matrho)>row(matrho)]
         rho[abs(rho)>=1] = 1
         veccor = sort(round(abs(rho),digits=dig))
         duplic = duplicated(veccor)
         vduplic = sort(unique(veccor[duplic]))
         vunic = setdiff(unique(veccor),vduplic)
         dtunic = vecdt[is.element(vecrho,vunic)]
         dtduplic = vecdt[is.element(vecrho,vduplic)]
         vmatch = match(vecrho,veccor,0)
         nboccur = diff(c(vmatch[vmatch>0],length(veccor)+1))
         nboccur = nboccur[nboccur>1]
         sdt[i] = 2*(m-1)*(sum(dtunic)+crossprod(nboccur,dtduplic))/(sampsize*(sampsize-1))  
      }
      return(sdt) 
   }

   nbfactors = function(dta,maxnbfactors=15,diagnostic.plot=FALSE,min.err=1e-03,verbose=FALSE){
      dig = 2
      m = ncol(dta)
      falist = vector(length=maxnbfactors+1,"list")
      falist[[1]] = list(B=matrix(0,ncol=1,nrow=m))
      falist[-1] = lapply(1:maxnbfactors,emfa,data=dta,min.err=min.err,verbose=verbose)
      Blist = lapply(falist,function(fa,m) matrix(fa$B,nrow=m),m=m)
      sdt = VarInflation(dta,Blist,maxnbfactors,dig,verbose)
      if (diagnostic.plot) {
         dev.new()
         plot(0:maxnbfactors,sdt,ylab="Variance Inflation Criterion",xlab="Number of factors",bty="l",
            lwd=1.25,type="b",pch=16,cex.lab=1.25,cex=1.25,cex.axis=1.25)
      }
      if (which.min(sdt)==1) nbf = 0
      if (which.min(sdt)>1) { 
         jumps = -diff(sdt)/sdt[-length(sdt)]
         opt = max((1:maxnbfactors)[jumps>0.05]) 
         nbf = opt
         #gain = sdt[1]-sdt[opt+1]
         #if (any(abs(sdt[1:(opt+1)]-sdt[opt+1])<0.10*gain)) nbf = min((0:length(sdt))[abs(sdt[1:(opt+1)]-sdt[opt+1])<0.10*gain])  
      }
      list(criterion=sdt,optimalnbfactors=nbf)
   }

   method = match.arg(method,choices=c("BH","holm","hochberg","hommel","bonferroni","BY","fdr","none"))
   if (is.null(design0)) design0 = matrix(1,nrow=nrow(dta),ncol=1)
   erpdta = as.matrix(dta)
   design = as.matrix(design)
   design0 = as.matrix(design0)
   if (typeof(erpdta)!="double") stop("ERPs should be of type double")
   if (nrow(erpdta)!=nrow(design)) stop("dta and design should have the same number of rows")
   if (nrow(erpdta)!=nrow(design0)) stop("dta and design0 should have the same number of rows")
   if (ncol(design)<=ncol(design0)) stop("design0 should have fewer columns than design")
   idsignal = NULL
   for (j in 1:ncol(design)) {
      cj = apply(design0,2,function(x,y) all(x==y),y=design[,j])
      if (all(!cj)) idsignal = c(idsignal,j)
   }
   if (length(idsignal)!=(ncol(design)-ncol(design0))) stop("the null model design0 should be nested into the non-null model design")
   if (typeof(alpha)!="double") stop("alpha should be of type double")
   if ((alpha<=0)|(alpha>=1)) stop("alpha should be in ]0,1[, typically 0.05")
   if (typeof(pi0)!="double") stop("pi0 should be of type double")
   if ((pi0<=0)|(pi0>1)) stop("pi0 should be in ]0,1]")

   if (length(s0)==1) stop("s0 should be either NULL, or of length larger than 2")
   frames = 1:ncol(erpdta)
   if (is.null(s0)) fs0i = integer(0)
   if (length(s0)>2) fs0i = s0   
   if (length(s0)==2) fs0i = which((frames<=s0[1]*diff(range(frames)))|(frames>=s0[2]*diff(range(frames))))

   n = nrow(erpdta)                  # sets the value for n
   T = ncol(erpdta)                  # sets the value for T 

   pdesign = solve(t(design)%*%design)%*%t(design)
   Proj = design%*%pdesign
   Proj0 = design0%*%solve(t(design0)%*%design0)%*%t(design0)
   rdf1 = nrow(design)-ncol(design)
   rdf0 = nrow(design0)-ncol(design0)

   beta = (pdesign%*%erpdta)[idsignal,]
   if (length(idsignal)==1) beta = matrix(beta,nrow=1)

   if (sum(is.element(nbfmax,0:rdf1))!=1) {
      warning(paste("nbfmax should be an integer in [0,",rdf1,"]",sep=""))
      nbfmax = rdf1
   }

   res1 = erpdta-Proj%*%erpdta
   scer1 = as.vector(t(rep(1,n))%*%res1^2)
   res0 = erpdta-Proj0%*%erpdta
   scer0 = as.vector(t(rep(1,n))%*%res0^2)
   fstat = ((scer0-scer1)/(rdf0-rdf1))/(scer1/rdf1) 
   pval = pf(fstat,df1=rdf0-rdf1,df2=rdf1,lower.tail=FALSE)
   
   if (is.null(nbf)) {
      nbfact = nbfactors(scale(res1),maxnbfactors=nbfmax,diagnostic.plot=wantplot,min.err=min.err,verbose=verbose)  
      nbf = nbfact$optimalnbfactors 
   }

   if (verbose) print(paste("AFA with", nbf,"factors"))

   if ((nbf>0)&(maxiter>0)) {

      if (length(s0)>=2) {
      
         err = 1
         betaSig0 = beta
         betaSig = beta

         R = diag(ncol(design))
         R = R[-(1:ncol(design0)),]
         if ((ncol(design)-ncol(design0))==1) R = matrix(R,nrow=1)
         coef1 = pdesign%*%erpdta+solve(t(design)%*%design)%*%t(R)%*%solve(R%*%solve(t(design)%*%design)%*%t(R))%*%(betaSig-R%*%pdesign%*%erpdta)

         k=0

         while ((err>min.err)&(k<=10)) {
            k=k+1
            res = erpdta - design%*%coef1
            fares = emfa(scale(res,center=TRUE,scale=FALSE),nbf=nbf,min.err=min.err,verbose=verbose)
            Psi = fares$Psi
            B = fares$B
            B0 = B[fs0i,]
            if (nbf==1) B0 = matrix(B0,ncol=1)
            iSxx = ifa(Psi[fs0i],B0)
            Sxy = (B0%*%t(B[-fs0i,]))
            betaSigma0 = iSxx$iS%*%Sxy
            beta0 = beta
            beta0[,-fs0i] = beta[,fs0i]%*%betaSigma0 
            betaSig = beta-beta0  
            for (i in 1:nrow(betaSig)) betaSig[i,] = ksmooth(frames,betaSig[i,],bandwidth=0.01*diff(range(frames)),x.points=frames)$y
            err = sqrt(mean((betaSig0-betaSig)^2)) 
            betaSig0 = betaSig 
            coef1 = pdesign%*%erpdta+solve(t(design)%*%design)%*%t(R)%*%solve(R%*%solve(t(design)%*%design)%*%t(R))%*%(betaSig-R%*%pdesign%*%erpdta)
         }
   
         res = erpdta - design%*%coef1   
         fa = emfa(scale(res,center=TRUE,scale=TRUE),nbf=nbf,min.err=min.err,verbose=verbose) 
         Psi = fa$Psi 
         B = fa$B 
         sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
         G = solve(diag(nbf) + sB %*% t(sB))
         sB = scale(t(B), center = FALSE, scale = Psi)
         Factors = scale(res,center=TRUE,scale=TRUE)%*%t(sB)%*%t(G)

         designz0 = cbind(design[,-idsignal],Factors)
         Projz0 = designz0%*%solve(t(designz0)%*%designz0)%*%t(designz0)
         designz1 = cbind(design[,-idsignal],Factors,design[,idsignal])
         idsignalz1 = idsignal+nbf
         pdesignz1 = solve(t(designz1)%*%designz1)%*%t(designz1)
         rdf0 = nrow(designz0)-ncol(designz0)
         rdf1 = nrow(designz1)-ncol(designz1)

         R = diag(ncol(designz1))
         R = R[-(1:ncol(designz0)),]
         if ((ncol(designz1)-ncol(designz0))==1) R = matrix(R,nrow=1)
         coef1 = pdesignz1%*%erpdta+solve(t(designz1)%*%designz1)%*%t(R)%*%solve(R%*%solve(t(designz1)%*%designz1)%*%t(R))%*%(betaSig-R%*%pdesignz1%*%erpdta)
         subjcoef = coef1[1:ncol(designz0),]
         if (ncol(designz0)==1) subjcoef = matrix(subjcoef,nrow=1)
         scer0 = apply((erpdta-designz0%*%subjcoef)^2,2,sum)

         R = diag(ncol(designz1))
         R = R[(1:ncol(designz0)),]
         if (ncol(designz0)==1) R = matrix(R,nrow=1)

         coef1 = pdesignz1%*%erpdta+solve(t(designz1)%*%designz1)%*%t(R)%*%solve(R%*%solve(t(designz1)%*%designz1)%*%t(R))%*%(subjcoef-R%*%pdesignz1%*%erpdta)
         scer1 = apply((erpdta-designz1%*%coef1)^2,2,sum)

         fstat = ((scer0-scer1)/(rdf0-rdf1))/(scer1/rdf1) 
         pval = pf(fstat,df1=rdf0-rdf1,df2=rdf1,lower.tail=FALSE)
   
      }
      
      correctedpval = pi0*p.adjust(pval,method=method)
      fs0 = sort(unique(c(fs0i,which(pval>0.2)))) 
      fs1 = integer(0)

      betaSig0 = beta
      betaSig = beta
      iter = 0 

      while((!setequal(fs0,fs1))&(iter<maxiter)) {

         err = 1
         iter = iter + 1
         if (verbose) print(paste(iter,"/",maxiter," iterations",sep=""))

         betaSig0 = beta
         betaSig = beta

         R = diag(ncol(design))
         R = R[-(1:ncol(design0)),]
         if ((ncol(design)-ncol(design0))==1) R = matrix(R,nrow=1)
         coef1 = pdesign%*%erpdta+solve(t(design)%*%design)%*%t(R)%*%solve(R%*%solve(t(design)%*%design)%*%t(R))%*%(betaSig-R%*%pdesign%*%erpdta)

         k=0

         while ((err>min.err)&(k<=10)) {
            k=k+1
            if (length(fs0)<(ncol(erpdta)-1)) {
               res = erpdta - design%*%coef1
               fares = emfa(scale(res,center=TRUE,scale=FALSE),nbf=nbf,min.err=min.err,verbose=verbose)
               Psi = fares$Psi
               B = fares$B
               B0 = B[fs0,]
               if (nbf==1) B0 = matrix(B0,ncol=1)
               iSxx = ifa(Psi[fs0],B0)
               Sxy = (B0%*%t(B[-fs0,]))
               betaSigma0 = iSxx$iS%*%Sxy
               beta0 = beta
               beta0[,-fs0] = beta[,fs0]%*%betaSigma0 
               betaSig = beta-beta0  
               for (i in 1:nrow(betaSig)) betaSig[i,] = ksmooth(frames,betaSig[i,],bandwidth=0.01*diff(range(frames)),x.points=frames)$y
               err = sqrt(mean((betaSig0-betaSig)^2))
               betaSig0 = betaSig 
               coef1 = pdesign%*%erpdta+solve(t(design)%*%design)%*%t(R)%*%solve(R%*%solve(t(design)%*%design)%*%t(R))%*%(betaSig-R%*%pdesign%*%erpdta)
            }
            if (length(fs0)>=(ncol(erpdta)-1)) {
               betaSig = matrix(0,nrow=nrow(betaSig),ncol=ncol(betaSig))  
               for (i in 1:nrow(betaSig)) betaSig[i,] = ksmooth(frames,betaSig[i,],bandwidth=0.05*diff(range(frames)),x.points=frames)$y
               err = sqrt(mean((betaSig0-betaSig)^2)) 
               betaSig0 = betaSig 
               coef1 = pdesign%*%erpdta+solve(t(design)%*%design)%*%t(R)%*%solve(R%*%solve(t(design)%*%design)%*%t(R))%*%(betaSig-R%*%pdesign%*%erpdta)
            }
         }

         res = erpdta - design%*%coef1   
         fa = emfa(scale(res,center=TRUE,scale=TRUE),nbf=nbf,min.err=min.err,verbose=verbose) 
         Psi = fa$Psi 
         B = fa$B 
         sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
         G = solve(diag(nbf) + sB %*% t(sB))
         sB = scale(t(B), center = FALSE, scale = Psi)
         Factors = scale(res,center=TRUE,scale=TRUE)%*%t(sB)%*%t(G)

         designz0 = cbind(design[,-idsignal],Factors)
         Projz0 = designz0%*%solve(t(designz0)%*%designz0)%*%t(designz0)
         designz1 = cbind(design[,-idsignal],Factors,design[,idsignal])
         idsignalz1 = idsignal+nbf
         pdesignz1 = solve(t(designz1)%*%designz1)%*%t(designz1)
         rdf0 = nrow(designz0)-ncol(designz0)
         rdf1 = nrow(designz1)-ncol(designz1)

         R = diag(ncol(designz1))
         R = R[-(1:ncol(designz0)),]
         if ((ncol(designz1)-ncol(designz0))==1) R = matrix(R,nrow=1)
         coef1 = pdesignz1%*%erpdta+solve(t(designz1)%*%designz1)%*%t(R)%*%solve(R%*%solve(t(designz1)%*%designz1)%*%t(R))%*%(betaSig-R%*%pdesignz1%*%erpdta)
         subjcoef = coef1[1:ncol(designz0),]
         if (ncol(designz0)==1) subjcoef = matrix(subjcoef,nrow=1)
         scer0 = apply((erpdta-designz0%*%subjcoef)^2,2,sum)

         R = diag(ncol(designz1))
         R = R[(1:ncol(designz0)),]
         if (ncol(designz0)==1) R = matrix(R,nrow=1)

         coef1 = pdesignz1%*%erpdta+solve(t(designz1)%*%designz1)%*%t(R)%*%solve(R%*%solve(t(designz1)%*%designz1)%*%t(R))%*%(subjcoef-R%*%pdesignz1%*%erpdta)
         scer1 = apply((erpdta-designz1%*%coef1)^2,2,sum)

         fstat = ((scer0-scer1)/(rdf0-rdf1))/(scer1/rdf1) 
         pval = pf(fstat,df1=rdf0-rdf1,df2=rdf1,lower.tail=FALSE)
         correctedpval = pi0*p.adjust(pval,method=method)
         fs1 = fs0
         fs0 = sort(unique(c(fs0i,which(pval>0.2)))) 

      }

      beta = (pdesignz1%*%erpdta)[idsignalz1,]
      if (length(idsignal)==1) beta = matrix(beta,nrow=1)
 
   }

   if (is.null(pi0)) pi0 = pval.estimate.eta0(pval,diagnostic.plot=FALSE)
   correctedpval = pi0*p.adjust(pval,method=method)
   significant = which(correctedpval<=alpha)
 
   list(pval=pval,correctedpval=correctedpval,significant=significant,pi0=pi0,test=fstat,df1=rdf1,df0=rdf0,nbf=nbf,signal=beta,r2=(1-1/(1+fstat*((rdf0-rdf1)/rdf1))))

}
