###################################################################################
## rss
## Objective function
## this function can be minimized to find the parameter values that
## minimize the RSS for a given model and given data
## In our experience this is best done with the package DEoptim. 
## for pareto task use package mco

# dataSAS: SAS part of S(Q)
# paramSASQ: if SAS(Q) should be calculated parametrically
# gammaR: SAS part of G(r)
#            = "param" => parametric envelope is calculated
# avRes: averages results over avRes simulations
# wG: weight for G(r)
# wSAS: weight for SAS(Q)
# pareto: if TRRUE then fuction value is a vector of residuals 
# kind: "exact", "fast", "fastHist", "fast_av" (for polydisperse particles)
# SASscale: "normal", "log", weighted"
# scaleSq: calculated PDF or S(Q) function should be multiplied by scaleSq
# scaleSASq: calculated SAS(Q)  should be multiplied by scaleSASq
#            if gammaR is given it is divided by scaleSASq

rss <- function(par, dataG=NA, dataS=NA, dataSAS=NA,
                part=NA, type="neutron", simPar=NA,  
                PDF.fixed=list(), TotalScatt.fixed=list(), 
                parscale=NA, skel=NA, con=TRUE, oneDW=FALSE, 
                punish=FALSE, 
                gammaR=NA, rvector=NA, fixed=NA,
                wG=0.03, wSAS=0.03, avRes=1, pareto=FALSE){
  if(!is.na(parscale[1])) 
    par2 <- relist( par / parscale, skel)
  else 
    par2 <- relist(par, skel)
  if(is.list(fixed)) 
    par2[names(fixed)] <- fixed

  if(con && !(is.null(par2$r))) {
    if(length(c(par2$rcore)) == 1) ## test for non-param
      if(par2$rcore > par2$r)
        par2$rcore <- par2$r     
  }
  else if(punish) { ## only makes sense now for core shell
    if(length(c(par2$rcore,par2$r)) == 2)
      if(par2$rcore>par$r)
        return(10e15)
  }

  r <- if(is.null(par2$r)) NA else par2$r
  rsigma <- if(is.null(par2$rsigma)) NA else par2$rsigma
  bkg <- if(is.null(par2$bkg)) 0 else par2$bkg  
  shell <- if(is.null(par2$shell)) NA else par2$shell  
  rcore <- if(is.null(par2$rcore)) NA else par2$rcore 
  scaleSq <- if(is.null(par2$scaleSq)) 1 else par2$scaleSq
  scaleSASq <- if(is.null(par2$scaleSASq)) 1 else par2$scaleSASq
  pStack <- if(is.null(par2$pStack)) 0 else par2$pStack
  pDimer <- if(is.null(par2$pDimer)) 0 else par2$pDimer  
  delta <- if(is.null(par2$delta)) 0 else par2$delta
  box <- if(is.null(par2$box)) NA else par2$box
  ellipse <- if(is.null(par2$ellipse)) NA else par2$ellipse
 if(delta==0 && !is.na(PDF.fixed$delta))
    delta <- PDF.fixed$delta
  if(delta==0 && !is.na(TotalScatt.fixed$delta))
    delta <- TotalScatt.fixed$delta
  
  if(!is.null(par2$sigma))
    sig <- par2$sigma
  else{
	sig <- 0
	for(i in 1:length(simPar$atoms))
	  sig <- append(sig, simPar$atoms[[i]]$sigma)
	if(length(simPar$atomsShell)>0){
	  for(i in 1:length(simPar$atomsShell))
	    sig <- append(sig, simPar$atomsShell[[i]]$sigma)
	}
	sig <- sig[-1]
  }
  
  if(is.na(r))
    r<- simPar$r
  if(is.na(rsigma))
    rsigma<- simPar$rsigma
	
  rIn <- r
	
  if(length(par2$latticep)>0)
    latticep <- par2$latticep
  else 
    latticep <- simPar$latticep
  if(length(latticep)==3){ #works for wurtzite structure
    simPar$atoms[[2]]$base[1,3] <- latticep[3]
    simPar$atoms[[2]]$base[2,3] <- 0.5-latticep[3]
  }

  if(length(par2$latticepShell)>0)
    latticepShell <- par2$latticepShell
  else 
    latticepShell <- simPar$latticepShell
  
  if(length(latticepShell)==3){ #works for wurtzite structure
    simPar$atomsShell[[2]]$base[1,3] <- latticepShell[3]
    simPar$atomsShell[[2]]$base[2,3] <- 0.5-latticepShell[3]
  }
  
  if(!is.na(shell) && !is.na(rsigma)){
    if(!is.na(rcore)){
	  if(simPar$distr=="lognormal")
        core <- rlnorm(avRes, meanlog = log(rcore), sdlog = log(rsigma))
      else 
	    core <- exp(rnorm(avRes, log(rcore), log(rsigma)))
		
	  core <- sort(core)
      r <- core + shell
	  rcore <- core
	}
    if(!is.na(r)){
	  if(simPar$distr=="lognormal")
        r <- rlnorm(avRes, meanlog = log(r), sdlog = log(rsigma))
      else 
	    r <- exp(rnorm(avRes, log(r), log(rsigma)))

      r <- sort(r)
	  rcore <- r-shell
	}
  }
   
  if(is.na(shell) && !is.na(rsigma)){
  	if(simPar$distr=="lognormal")
      r <- rlnorm(avRes, meanlog = log(r), sdlog = log(rsigma))
    else 
	  r <- exp(rnorm(avRes, log(r), log(rsigma)))
		
    r <- sort(r)
  }
  kind <- TotalScatt.fixed$kind
  dr <- TotalScatt.fixed$dr
  del <- TotalScatt.fixed$del
  eps <- TotalScatt.fixed$eps
  SASscale <- TotalScatt.fixed$SASscale
  paramSASQ <- TotalScatt.fixed$paramSASQ
    
  gQ_SAS_av <- gQ_av <- mod_av <- gQ_SAS <- gQ <-  mod <- 0  
  simulate_particle <- TRUE
  
  ####################
  #
  # !MAIN CYCLE OVER AVRES!
  #
  ####################
  for(j in 1:avRes){
    p <- 0
    if(kind !="fast_av" && !is.na(rsigma)){
      r_tmp <- r[j]
	  rcore_tmp <- rcore[j]
	}
	else{
	  r_tmp <- r
      rcore_tmp <- rcore
	} 	 				
	
    if(!is.na(part[1]) && j==1)   # if on the first step we already have particle
      simulate_particle <- FALSE  # we don't have to simulate new ones...
      
    if(simulate_particle)
      part <- simPart(atoms=simPar$atoms, sym=simPar$sym, latticep=latticep, r=r_tmp,
	                  atomsShell=simPar$atomsShell, symShell=simPar$symShell, 
					  latticepShell=latticepShell, rcore=rcore_tmp, 
				      rcenter=simPar$rcenter, box=box, ellipse=ellipse, 
					  move=simPar$move, rotShell=simPar$rotShell, rcenterShell=simPar$rcenterShell, 
				      pDimer=pDimer,  pStack=pStack)
					  
    if(oneDW && j==1)
      sig <- rep(sig, attributes(part)$nAtomTypes)
#========================================================	 
	if(!is.na(dataS[1])){
      cat("Calculating Scattering function... \n")
      gQ <- calcTotalScatt(part, type=type, kind=kind, 
	                        dQ=TotalScatt.fixed$dQ,
                            minQ=TotalScatt.fixed$minQ,
                            maxQ=TotalScatt.fixed$maxQ,	  
	                        scatterLength=TotalScatt.fixed$scatterLength, 
                            scatterFactor=TotalScatt.fixed$scatterFactor,
                            sigma=sig, 
                            delta=delta, n=TotalScatt.fixed$n, 
						    dr=dr, del=del, eps=eps)$gQ	
    }
					
#========================================================
	if(!is.na(dataG[1])){
	  cat("Calculating  PDF... \n")	
      mod <- calcPDF(part,type=type,
                   dr=PDF.fixed$dr,
                   minR=PDF.fixed$minR,
                   maxR=PDF.fixed$maxR,
                   scatterLength=PDF.fixed$scatterLength,
                   scatterFactor=PDF.fixed$scatterFactor,
                  )
      mod <- broadPDF(mod, sigma=sig, delta=delta, n=PDF.fixed$n, nAtomTypes = attributes(part)$nAtomTypes)$gr
	  mod[which(is.na(mod))] <- 0
      if(any(is.na(mod)))
        mod <- rep(10e15, length(rvector))
      mod[which(is.infinite(mod))] <- 0  
	}
	
    TotalScattSAS=list(dQ=TotalScatt.fixed$dQ_SAS, minQ=TotalScatt.fixed$minQ_SAS, maxQ=TotalScatt.fixed$maxQ_SAS, scatterLength=TotalScatt.fixed$scatterLength)
#========================================================	
	if(!is.na(dataSAS[1]) && (!paramSASQ)){
      cat("Calculating SAS ... \n")     
	  gQ_SAS <- calcTotalScatt(part,type=type, kind=kind,
                          dQ=TotalScattSAS$dQ,
                          minQ=TotalScattSAS$minQ,
                          maxQ=TotalScattSAS$maxQ,
	                      scatterLength=TotalScatt.fixed$scatterLength, 
                          scatterFactor=TotalScatt.fixed$scatterFactor,
                          sigma=sig, 
                          delta=delta, n=TotalScatt.fixed$n, 
						  dr=dr, del=del, eps=eps)$gQ	
	}


#=============================================================
    gQ_SAS_av <- gQ_SAS_av + gQ_SAS
    gQ_av <- gQ_av + gQ
    mod_av <- mod_av + mod
  }	
  gQ_SAS <- gQ_SAS_av/avRes
  gQ <- gQ_av/avRes
  mod <- mod_av/avRes	

  if(!is.na(dataG[1])){
    if(!any(is.na(rvector)))
      mod <- 4*pi*rvector*mod
  		 
    if(!any(is.na(gammaR))){
	  if(gammaR[1]=="param"){
		if(!is.na(shell))
		  sasvector <- GrSASCS(rvector, Rcore=rcore, Rpart=r, latticep=latticep, latticepShell=latticepShell, 
		                       N1=PDF.fixed$sym, N2=PDF.fixed$sym, sym=simPar$sym, symShell=simPar$symShell)
		else
		  sasvector <- GrSAS(rvector, Rcore=rcore, Rpart=r, latticep=latticep, latticepShell=latticepShell, 
		                       N1=PDF.fixed$sym, N2=PDF.fixed$sym, sym=simPar$sym, symShell=simPar$symShell)
		
		mod <- mod  - sasvector
	  }
	  else
	    mod <- mod  - gammaR/scaleSASq

	} 
	if(PDF.fixed$termRip && !any(is.na(mod)))
	  mod <- termRip(mod, Qmax=PDF.fixed$Qmax, dr=PDF.fixed$dr, maxR=PDF.fixed$maxR, maxRTermRip=PDF.fixed$maxRTermRip)
    if(!is.na(PDF.fixed$Qdamp))
	  mod <- mod*exp(-(PDF.fixed$Qdamp*rvector)^2/2)
  }	
		
  if(!is.na(dataS[1]) && TotalScatt.fixed$convolution && !any(is.na(gQ)))
	gQ <- gaussConvol(SQ=gQ, Q=seq(TotalScatt.fixed$minQ, TotalScatt.fixed$maxQ, TotalScatt.fixed$dQ), Qdamp=TotalScatt.fixed$Qdamp)
 
  if(!is.na(dataSAS[1]) && paramSASQ){
    cat("Calculating  parametric SAS(Q)... \n")
	Q=seq(TotalScatt.fixed$minQ_SAS, TotalScatt.fixed$maxQ_SAS, TotalScatt.fixed$dQ_SAS)
	if(!is.na(rsigma))
	  gQ_SAS <- IqSASP(Q=Q, shell=shell, Rpart=rIn, latticep=latticep, latticepShell=latticepShell, 
	                   scatterLength=c(TotalScatt.fixed$f1, TotalScatt.fixed$f2), pDimer=pDimer, 
			           N1=TotalScatt.fixed$N1, N2=TotalScatt.fixed$N2,
			           sym=simPar$sym, symShell=simPar$symShell, rsigma=rsigma) 
    else
	  gQ_SAS <- IqSAS(Q=Q, Rcore=rcore, Rpart=r, latticep=latticep, latticepShell=latticepShell,
	                  scatterLength=c(TotalScatt.fixed$f1, TotalScatt.fixed$f2), 
                      N1=TotalScatt.fixed$N1, N2=TotalScatt.fixed$N2,  
					  pDimer=pDimer, sym=simPar$sym, symShell=simPar$symShell)				   
  }
##
  mod <- scaleSq*mod
  if(!is.na(dataG[1]))	
    rssG <- sqrt(sum((mod-dataG)^2)/( sum(dataG^2))) * wG
  else 
    rssG <- 0  
##	
  gQ <- scaleSq*gQ
  if(!is.na(dataS[1]))		
    rssS <- sqrt(sum((gQ-dataS)^2)/(sum(dataS^2))) 
  else 
    rssS <- 0
##		
  if(!is.na(dataSAS[1])){
    dataSAS <- dataSAS - bkg
    gQ_SAS <- scaleSASq*gQ_SAS
	
    if(SASscale=="normal")
	  rssSAS <- sqrt(sum((gQ_SAS-dataSAS)^2)/(sum(dataSAS^2))) * wSAS	
	if(SASscale=="log")
	  rssSAS <- sqrt(sum((log(abs(gQ_SAS))-log(abs(dataSAS)))^2)/sum( (log(abs(dataSAS)))^2 )) * wSAS	
	if(SASscale=="weighted")
      rssSAS <- sqrt(sum((gQ_SAS/dataSAS-1)^2)/(length(dataSAS))) * wSAS		
  }
  else 
    rssSAS <- 0
	
  if(is.na(rssS)) rssS <- 1e15
  if(is.na(rssSAS)) rssSAS <- 1e15
  if(is.na(rssG)) rssG <- 1e15

  include <- logical(3)
  include[1] <- !is.na(dataSAS[1])
  include[2] <- !is.na(dataG[1])
  include[3] <- !is.na(dataS[1])
  
  res <- numeric(3)
  if(include[1])
    res[1] <- rssSAS	
  else 
    res[1] <- 0

  if(include[2])
    res[2] <- rssG
  else 
    res[2] <- 0

  if(include[3])
    res[3] <- rssS
  else 
    res[3] <- 0  
	
  if(pareto)
    res[which(include==TRUE)]
  else
    sqrt(res[1]^2 + res[2]^2 + res[3]^2)
}


