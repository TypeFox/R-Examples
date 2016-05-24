#'## Lindeman Trophic Aggregations
#'## Singh P. July 2014.
#'## Source Ulanowicz and Kay 1991. Environmental Software 6:131-142.
#'## ----------------------------------------------------------------------
enaTroAgg <- function (x){
  if (class(x) != "network") {
    stop("x is not a network class object")
  }

                                        # Initials

  liv <- x %v% "living"     ##Living vector
  nl = sum(liv)             ##No. of living nodes
  N <- length(liv)
                                        #Living vector check
  liv2<-rep(FALSE,N)
  liv2[1:nl]<-rep(TRUE,nl)
  if(identical(liv,liv2) == FALSE) {
      stop('Non-living nodes must be at the end of the list.')
  }

  flow <- as.matrix(x, attrname = 'flow')
  XCHNGE<-flow
  Feeding_Cycles   <- cycliv(x)
  XCHNGE[1:nl,1:nl] <- Feeding_Cycles$ResidualFlows
  Ti <- x %v% "input"       ##AINPUT
  T <- Ti + apply(XCHNGE, 2, sum)

  exp<-x %v% "export"; exp[is.na(exp)] <- 0
  res<-x %v% "respiration"; res[is.na(res)] <- 0
  # ---------------------------------------------------------

                                        # Determining In-Migration and Obligate Producers
  BINPUT <- x %v% 'input'
  if(identical(XCHNGE,flow[1:nl,1:nl])) {print('Cycle free feeding transfers')}
  NMIG <- NPRM <- 0 ###NMIG - In-migration of Heterotrophs; NPRM - no. of obligate primary producers
  CANON <- rep(0,N)
  ###A compartment is an obligate producer iff it receives sustenance from no other compartment
  for(NP in 1:nl) {
  	if(Ti[NP]<=0){next}
  	for(i in 1:N) {
            ##Otherwise the input represents in-migration
  		if(XCHNGE[i,NP]<=0){next}
  		NMIG=NMIG+1
                ##Record the Migratory Input Temporary in CANON for use below
  		CANON[NP]=Ti[NP]
  		break
  	}
  	NPRM = NPRM+1
  }
  if(NPRM<=0){
  	warning("No Unambiguous primary producers found! All inputs assumed to be primary production!")
  	NMIG <- 0
  	CANON<- rep(0,N)
  	}
  mig.input<-rep(0,nl)
  if(NMIG>0) {
  	###Migratory Inputs. To be treated as Non-Primary Inflows
  	mig.input <- CANON[1:nl]
  }
  #---------------------------------------------------------

                                        #Recreate the matrix of feeding coefficients without the migratory inputs
  TL <- rep(0,N)
  FEED <- flow*0
  for(i in 1:N) {
  	##Store Throughputs in Vector TL
  	TL[i]<-0
  	for(j in 1:N) {TL[i]<-TL[i]+XCHNGE[j,i]}
  	TL[i]<-TL[i]+Ti[i]-CANON[i]
  	if(TL[i]<=0){TL[i]=1}
  	for(j in 1:N) {FEED[j,i]<-XCHNGE[j,i]/TL[i]}
  	BINPUT[i]<-(Ti[i]-CANON[i])/TL[i]
  }
  #---------------------------------------------------------

                                        #Create Lindeman Transformation Matrix
  CANON<- rep(0,N)
  A<-flow*0
  A[1,1:nl]<-BINPUT[1:nl]
  for(k in 2:nl) {
  	KM1=k-1
  	for(l in 1:nl) {
  		if((k<=2)&&(nl<N)) {
  			for(K2 in (nl+1):N) {A[k,l]<-A[k,l]+FEED[K2,l]}
  		}
  		for(j in 1:nl) {
  			A[k,l]<-A[k,l]+(A[KM1,j]*FEED[j,l])
  		}
  	}
  }
  if(nl<N) {for (i in (nl+1):N) {A[N,i]<-1}}
  ## A is the required Lindeman transformation matrix
  rownames(A) <- 1:N

                                        # 2. Effective Trophic Levels
  MF = matrix(1:N, nrow=N, ncol=N, byrow = 'FALSE')
  etl = rep(1,N)
  etl[1:nl] = apply((MF*A)[1:nl,1:nl],2,sum)


  ci <- Ti
  ci <- A %*% Ti
  ci <- as.vector(ci)


                                        # 3. Canonical Exports
  cel = exp
  ce = A %*% exp
  ce1 = as.vector(ce)

                                        # 4. Canonical Respirations
  crl = res
  cr = A%*%res
  cr1=as.vector(cr)

                                        # 5. Grazing Chain
  gc <- rep(0,nl)
  gc[1] <- sum(A[1,]*T)
  AT = A %*% flow %*% t(A)
  gc[2:nl]=apply(AT[1:(nl-1),1:nl,drop=FALSE],1,sum)

                                        # 6. Returns to Detrital Pool
  rtd <- AT[1:nl,N]
  rtd <- as.vector(rtd)

                                        # 7. Detrivory
  dtry <- sum(AT[N,1:nl])

                                        # 8. Input to Detrital Pool
  U <- A %*% Ti
  dinp<-0
  if(nl<N) {  dinp <- sum(U[(nl+1):N]) }

                                        # 9. Circulation within Detrital Pool
  dcir <- AT[N,N]

                                        # 10. Lindeman Spine
  ls = gc

  ls[1] = sum(rtd[2:nl]) + gc[1] + dinp
  ls[2] = gc[2]+dtry
                                        # 11. Trophic Efficiencies
  te=ls
  for(i in 1:nl){
    if(te[i]<=0){break}
    te[i]=te[i+1]/te[i]
  }
  te[is.na(te)] <- 0

                                        # Output Listing
  Detrivory<-dtry; DetritalInput<-dinp; DetritalCirc<-dcir
  ns <- cbind(Detrivory, DetritalInput, DetritalCirc, Feeding_Cycles$ns)
  if(NMIG>0) {
  	out <- list(Feeding_Cycles=Feeding_Cycles[1:(length(Feeding_Cycles)-1)], A = A[1:nl,1:nl], ETL = etl, M.Flow = mig.input, CI = ci, CE = ce1, CR = cr1, GC = gc, RDP = rtd, LS = ls,TE = te, ns=ns)
  }
  else{
  	out <- list(Feeding_Cycles=Feeding_Cycles[1:(length(Feeding_Cycles)-1)], A = A[1:nl,1:nl], ETL = etl, CE = ce1, CR = cr1, GC = gc, RDP = rtd, LS = ls,TE = te, ns=ns)

  }

  return(out)

  }#End of Function troAgg




