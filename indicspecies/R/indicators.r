indicators <- function (X, cluster, group, func="IndVal", max.order = 5, max.indicators=NULL, At = 0, Bt=0, sqrtIVt =0, nboot=0, alpha=0.05, XC = TRUE, enableFixed = FALSE, verbose=FALSE) {
	                 
  func <- match.arg(func, c("IndVal", "IndVal.g"))                                                                                                             
  if(sum(is.na(cluster))>0) stop("Cannot deal with NA values. Remove and run again.")
  if(sum(is.na(X))>0) stop("Cannot deal with NA values. Remove and run again.")
  nsites = nrow(X)

  cluster= as.factor(cluster)
  group <- as.character(group)
  group <- match.arg(group, levels(cluster))                                                                                                             
  if(verbose) cat(paste("Target site group: ",group,"\n", sep=""))
  group.vec <- cluster==group
  group.vec[is.na(group.vec)] = FALSE

  #Get species names
  spplist = names(X)
  if(verbose) cat(paste("Number of candidate species: ",length(spplist),"\n", sep=""))
  if(length(spplist)==1) stop("At least two species are necessary.")
  
  #Select rows that contain the species or the group
  ng = sum(group.vec)
  if(verbose) cat(paste("Number of sites:",nsites,"\n"))
  if(verbose) cat(paste("Size of the site group:",ng,"\n"))
  
  #Study frequency
  freq = colSums(ifelse(X[group.vec,]>0,1,0))/sum(group.vec)
  if(enableFixed) {
    fixedSpecies = (freq==1)
  } else {
    fixedSpecies = rep(FALSE, length(spplist))
  }
  numFixed = sum(fixedSpecies)
  if(verbose & enableFixed) cat(paste("Number of fixed species:",numFixed,"\n"))
  fixedPos = which(fixedSpecies)
  if(verbose & numFixed>0) print(fixedPos)
  

  evalCombs<-function(spvec, dvec, verbose=FALSE) {
    comblist<-vector("list",0)
    sc.ab<-apply(X[,spvec, drop=FALSE],1,min)
    if(sum(sc.ab)>0) {
      scg = sc.ab[group.vec]
      if(func=="IndVal.g") {
        mg = (sum(scg)/ng)
        A = mg/sum(tapply(sc.ab,cluster, "mean"))
      } else {
        A = sum(scg)/sum(sc.ab)
      }
      B = sum(scg>0)/ng
      sqrtIV = sqrt(A*B)
      if(A>=At & B>=Bt & sqrtIV>=sqrtIVt) {  
        comblist<-c(comblist,list(spvec))
      }
      if(B>Bt && (length(spvec)-numFixed)<max.order && length(dvec)>0) {##If B is not too small and order is allowed we can explore further
        for(j in 1:length(dvec)) {
          if(verbose) cat(paste("Starting species ",dvec[j],"..."))          
          comblist<-c(comblist,evalCombs(c(spvec,dvec[j]), dvec[-(1:j)], verbose=FALSE))
          if(verbose) cat(paste(" accepted combinations:",length(comblist),"\n"))
        }      
      }
    } 
    return(comblist)
  }
  k = length(spplist)-numFixed #Number of species to combine
  if(verbose & enableFixed) cat(paste("Number of species to combine: ",k,"\n", sep=""))
  veck = (1:length(spplist))[!fixedSpecies] #Vector of species indices to combine
  
  comblistDef<-vector("list",0)
  if(length(fixedPos)>0) {
    comblistDef<-c(comblistDef,evalCombs(fixedPos, veck,verbose=TRUE))
  } else {
    for(i in 1:length(veck)) {
      cat(paste("Starting species ",veck[i],"..."))
      comblistDef<-c(comblistDef,evalCombs(c(veck[i],fixedPos), veck[-(1:i)], verbose=FALSE))
      cat(paste(" accepted combinations:",length(comblistDef),"\n"))
    }    
  }
  
  #Create structures to store data
  nc = length(comblistDef)
  if(verbose) cat(paste("Number of valid combinations: ",nc,"\n", sep=""))
  if(nc==0) return()
  trim =FALSE
  if(!is.null(max.indicators)) {
    if(nc>max.indicators) {
      if(verbose) cat(paste("Maximum number of valid combinations exceeded.\n", sep=""))    
      trim = TRUE
    }
  }
  Astat = numeric(nc)
  Bstat = numeric(nc)
  if(trim) {
    for(i in 1:nc) {
      spvec = as.numeric(comblistDef[[i]])
      sc.ab <-apply(X[,spvec, drop=FALSE],1,min)
      scg = sc.ab[group.vec]
      if(func=="IndVal.g") {
        mg = (sum(scg)/ng)
        Astat[i] = mg/sum(tapply(sc.ab,cluster, "mean"))
      } else {
        Astat[i] = sum(scg)/sum(sc.ab)
      }
      Bstat[i] = sum(scg>0)/ng
    }
    sqrtIVstat = sqrt(Astat*Bstat)    
    sel = order(sqrtIVstat,decreasing = TRUE)[1:max.indicators]
    Astat = Astat[sel]
    Bstat = Bstat[sel]
    sqrtIVstat = sqrtIVstat[sel]
    nc = max.indicators
    Cvalid<-as.data.frame(matrix(0,nrow=nc,ncol=length(spplist)))
    names(Cvalid)<-spplist
    XC = data.frame(matrix(0, nrow=nrow(X), ncol=nc))
    for(i in 1:nc) {
      spvec = as.numeric(comblistDef[[sel[i]]])
      XC[,i] <-apply(X[,spvec, drop=FALSE],1,min)
      Cvalid[i,spvec] = 1
    }
  } else {
    Cvalid<-as.data.frame(matrix(0,nrow=nc,ncol=length(spplist)))
    XC = data.frame(matrix(0, nrow=nrow(X), ncol=nc))
    names(Cvalid)<-spplist
    for(i in 1:nc) {
      spvec = as.numeric(comblistDef[[i]])
      sc.ab <-apply(X[,spvec, drop=FALSE],1,min)
      XC[,i]<-sc.ab
      scg = sc.ab[group.vec]
      if(func=="IndVal.g") {
        mg = (sum(scg)/ng)
        Astat[i] = mg/sum(tapply(sc.ab,cluster, "mean"))
      } else {
        Astat[i] = sum(scg)/sum(sc.ab)
      }
      Bstat[i] = sum(scg>0)/ng
      Cvalid[i,spvec] = 1
    }
    sqrtIVstat = sqrt(Astat*Bstat)    
  }
  
  #Remove species that do not appear in any valid combination
  selSpp = colSums(Cvalid)>0
  if(verbose) cat(paste("Number of remaining species:",sum(selSpp),"\n"))
  Cvalid = Cvalid[,selSpp]
  nspp <- sum(selSpp)
  
  
  #Calculate bootstrap confidence intervals for sensitivity and ppp of valid combinations
  if(nboot>0) {
  	  if(nboot<100) nboot=100 #Minimum of 100 bootstrap replicates
	  if(verbose) {
  			cat(paste("Calculating bootstrap confidence intervals"))
  	  }
	  dmbA = matrix(NA,nrow=nboot,ncol=nc)
	  dmbB = matrix(NA,nrow=nboot,ncol=nc)
	  dmbIV = matrix(NA,nrow=nboot,ncol=nc)
	  for(b in 1:nboot) {
	  	  if(b%%round(nboot/10)==0 && verbose) cat(".")
		  bi = sample(nsites,replace=TRUE)
		  ngb = sum(group.vec[bi])
		  XCB = as.matrix(XC[bi,])
		  XCBg = as.matrix(XCB[group.vec[bi],])
	  	  if(func=="IndVal.g") {
	  	  	kk <- colSums(apply(XCB,MARGIN=2,FUN=tapply,cluster[bi],"mean"))
		  	dmbA[b,] = colMeans(XCBg)/kk
  		  } else {
  			dmbA[b,] = colSums(XCBg)/colSums(XCB)
  		  }
		  dmbB[b,] = colSums(as.matrix(ifelse(XCBg>0,1,0)))/ngb
		  dmbIV[b,] = sqrt(dmbA[b,]*dmbB[b,])
	  }
	  if(verbose) cat(paste("\n"))
	  dmlowerA = rep(0,nc)
	  dmupperA = rep(0,nc)
	  dmlowerB = rep(0,nc)
	  dmupperB = rep(0,nc)
	  dmlowerIV = rep(0,nc)
	  dmupperIV = rep(0,nc)
	  for(i in 1:nc) {	
			sdmb = sort(dmbA[,i])			
			dmlowerA[i]=sdmb[(alpha/2.0)*nboot]
			dmupperA[i]=sdmb[(1-(alpha/2.0))*nboot]
			sdmb = sort(dmbB[,i])
			dmlowerB[i]=sdmb[(alpha/2.0)*nboot]
			dmupperB[i]=sdmb[(1-(alpha/2.0))*nboot]
			sdmb = sort(dmbIV[,i])
			dmlowerIV[i]= sdmb[(alpha/2.0)*nboot]
			dmupperIV[i]= sdmb[(1-(alpha/2.0))*nboot]
	  }
	  sA = as.data.frame(cbind(Astat,dmlowerA,dmupperA))
  	  names(sA) = c("stat", "lowerCI", "upperCI")
  	  row.names(sA) = row.names(Cvalid)
	  sB = as.data.frame(cbind(Bstat,dmlowerB,dmupperB))
  	  names(sB) = c("stat", "lowerCI", "upperCI")
  	  row.names(sB) = row.names(Cvalid)
	  sIV = as.data.frame(cbind(sqrtIVstat,dmlowerIV,dmupperIV))
  	  names(sIV) = c("stat", "lowerCI", "upperCI")
  	  row.names(sIV) = row.names(Cvalid)
  } else{
  	  sA = Astat
  	  sB = Bstat
  	  sIV = sqrtIVstat
  }
  
  result = list(group.vec =group.vec, candidates = spplist, finalsplist= spplist[selSpp], C=Cvalid, XC=XC, A=sA, B=sB, sqrtIV=sIV)
  class(result) = "indicators"
  return(result)
}


