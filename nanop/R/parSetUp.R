simPartPar <- function(sym = "fcc", symShell=NA,  latticep = 4.08, latticepShell=NA, 
                       atoms=list(), atomsShell=list(), distr = "lognormal",
                       r=NA, rsigma=NA, rcore=NA, box=NA, ellipse=NA, shell=NA,
                       rcenter=FALSE, move=TRUE, rotShell=FALSE, rcenterShell=FALSE){

  nAtomTypes <- length(atoms) + length(atomsShell)
  
  list(sym=sym, symShell=symShell, nAtomTypes = nAtomTypes, distr=distr,
       latticep=latticep, r=r,rsigma=rsigma, atoms=atoms, atomsShell=atomsShell, 
       latticepShell = latticepShell, box=box, ellipse=ellipse, shell=shell, 
       rcore = rcore, rcenter=rcenter, rotShell=rotShell, 
	   move=move, rcenterShell=rcenterShell)

}
	   
###########################################################################
	   
PDFPar <- function(dr=.01, minR=1, maxR=20, 
                   termRip=FALSE, Qmax=30, maxRTermRip=20,
                   scatterLength=NA, scatterFactor=NA,
				   n=2, delta=NA, N1=4, N2=4,
				   Qdamp=NA)
  list(dr=dr, minR=minR, maxR=maxR, 
       termRip=termRip, Qmax=Qmax, N1=N1, N2=N2,
	   Qdamp=Qdamp, n=n, delta=delta, maxRTermRip=maxRTermRip,
       scatterLength=scatterLength, scatterFactor=scatterFactor)

	   
	   
###########################################################################	   
	   
TotalScattPar <- function(dQ=.01, minQ=0.771, maxQ=35, minQ_SAS=0.001, maxQ_SAS=0.771, dQ_SAS=0.005,
                     dr = 0.001, del = 0.01, eps=1e-3,  kind="fastHist", SASscale="normal",        
					 convolution = FALSE, Qdamp=0.0457, delta=NA, n=2, paramSASQ=FALSE,  
                     scatterFactor=NA, scatterLength=NA, N1=4, N2=4, f1=5, f2=5)			 
  list(dQ=dQ, dQ_SAS=dQ_SAS, minQ=minQ,maxQ=maxQ, minQ_SAS=minQ_SAS, maxQ_SAS=maxQ_SAS, kind=kind, SASscale=SASscale, 
        dr=dr,  del=del, eps=eps, scatterLength=scatterLength, scatterFactor=scatterFactor, 
		convolution=convolution, Qdamp=Qdamp, delta=delta, n=n, paramSASQ=paramSASQ,
		N1=N1, N2=N2, f1=f1, f2=f2)
  
