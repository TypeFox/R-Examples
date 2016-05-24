## pipe.design version 0.4
## MJS 01/04/2016

if(getRversion() >= "2.15.1") globalVariables(c("y","z","ci"))


## N - sample size
## S - Number of trial simulations (default 1)
## c - cohort size
## theta - TTL
## pi - true p(DLT) matrix. If not specified only next recommended dose is returned.
## prior.med - prior median p(DLT) matrix
## prior.ss - prior sample size matrix
## strategy - dose escalation strategy. Options are 
	## "ss" - Sample size strategy. Choose from the admissible set of doses based on inverse of the sample size
	## "ss-random" - Random sample size strategy. Choose from the admissible set of doses based on randomisation, weighted by the inverse of the sample size
## admis - Choice of admissible doses around the MTC. Choices are
	## "adjacent" - All doses that are adjacent the MTC are admissible (even corner doses)
	## "closest" - All doses that are closest to the MTC
## constraint - dose-skipping constraint (defaults to no constraint). Options are
	## "neighbouring" - any dose combination can be chosen up to one dose level above OR BELOW current drug A and drug B levels (i.e. neighbouring current dose)
	## "no.dose.skip" - any dose combination can be chosen up to one dose level above ANY previously experimented drug A and drug B levels
	## "neighbouring-nodiag" - As "neighbouring" but prohibiting diagonal escalation (diagonal de-escalation is still allowed)
	## "no.dose.skip-nodiag" - As "no.dose.skip" but prohibiting diagonal escalation (diagonal de-escalation is still allowed)
## epsilon - should there be a safety constraint imposed if the posterior probability(dose > MTC)>epsilon. This averages over posterior distribution of MTC contours
	## Defaults to NULL
	## Any number is the posterior tail probability for which any dose that exceeds this is not experimented on.
## mode - Testing mode. Options are
	## "sim" - no testing (used for simulations) [the default]
	## "nodlt" - Every patient is DLT free
	## "alldlt"  - Every patient suffers a DLT
## data (Optional) A named data frame giving information about dose and toxicity from previously recruited patients. If missing, then it is assumed that no data have thus far been collected. Contains the following variables:
    ##patient Recruited patient numbers, 1,...,n
    ##doseA Dose levels of Drug A 
    ##doseB Dose levels of Drug B
    ##tox An indicator variable for each patient (1=toxicity, 0=no toxicity)
## a - Matrix of values for prior Beta(a,b) distributions [Alternative to specifying prior.med and prior.ss]
## b - Matrix of values for prior Beta(a,b) distributions [Alternative to specifying prior.med and prior.ss]
## alternate - deescalate if above the MTC, escalate if below - defaults to FALSE
## uppertox.constraint - should an upper toxicity safety constraint be imposed?
	## Defaults to NULL
	## Any number is the upper toxicity level from which no dose should be experimented or recommended. 
	## N.B. This constraint is NOT documented in the Mander and Sweeting 2015 paper and its performance is often substandard to using the safety constraint that weights over all possible MTCs, due to rigidity. It should therefore be used with caution. 
## stop - an alternative stopping rule, whereby we stop if posterior prob. of being > TTL at lowest dose is > stop
	## Defaults to null
	## Any number >0 and <1
	## N.B. Not documented in Mander and Sweeting 2015.
## non.admissible - Dose combinations that are never to be experimented on 
## seed - Seed to set so results can be reproduced
## contour.select - TO COME IN FUTURE VERSION
## reweight - TO COME IN FUTURE VERSION
## R - TO COME IN FUTURE VERSION
## P - TO COME IN FUTURE VERSION

pipe.design<-function(N=dim(data)[1]+1,S=1,c,theta,pi=NULL,prior.med=NULL,prior.ss=NULL,strategy,admis,constraint=NULL,epsilon=NULL,mode="sim",data=matrix(nrow=0,ncol=0),a=NULL,b=NULL,alternate=FALSE,uppertox.constraint=NULL,stop=NULL,non.admissible=NULL,seed=NULL){

  ## Contour select option that will be used in future version
  contour.select="mode"  
  ## Reweighting options that will be used in future version
  reweight=FALSE
  R=0
  P=0

  if(!is.null(seed)){
    set.seed(seed)
  }
  
	# Check that either prior.med and prior.ss are specified, or a and b
	if(is.null(prior.med) & is.null(a)){
		stop("Either `prior.med' and `prior.ss' or `a' and `b' must be specified")
	}
	# Dimensions of two-agent design
	if(!is.null(prior.med)){
		I=dim(prior.med)[1]
		J=dim(prior.med)[2]
	} else {
		I=dim(a)[1]
		J=dim(a)[2]
	}
	
	# Do some checks of inputs
	if(!(strategy %in% c("ss","ss-random"))) stop("strategy must be one of `ss' or `ss-random'")
## ,"psmooth" TO COME IN FUTURE VERSION  or `psmooth'"
	if(!(admis %in% c("adjacent","closest"))) stop("admis must be one of `adjacent' or `closest'")
  # Following check also allows constraint to take value "none" for ShinyPIPE application
  if(!is.null(constraint) & !(constraint %in% c("none","neighbouring","no.dose.skip","neighbouring-nodiag","no.dose.skip-nodiag"))) stop("constraint must be one of `neighbouring', `no.dose.skip', `neighbouring-nodiag' or `no.dose.skip-nodiag'")
	if(!(mode %in% c("sim","nodlt","alldlt"))) stop("mode must be one of `sim', `nodlt', or `alldlt'")
	if(mode=="sim" & is.null(pi)){
		# cat("`pi' must be specified to conduct a simulation study, only next recommended dose will be given \n")
		N<-dim(data)[1]+c
	}
	if(!is.null(uppertox.constraint)){
		if((theta>uppertox.constraint | uppertox.constraint>1)) stop("uppertox.constraint must be a number between theta and 1")
	}
	if(!is.null(epsilon)){
		if(epsilon<0 | epsilon>1) stop("epsilon must be a number between 0 and 1")
	}
	if(ceiling(N/c)!=floor(N/c)) stop("Total sample size `N' must be divisible by cohort size `c'")
	if (nrow(data)>0) {
        if (any(!(c("patient", "doseA", "doseB","tox") %in% names(data)))) 
            stop("data must have variables named 'patient', 'doseA', 'doseB' and 'tox'")
        data <- data[order(data$patient), ]
        if (any(data$patient != 1:dim(data)[1])) 
            stop("'patient' variable in data must be an ascending vector of positive integers")
        if (any(!(data$tox %in% c(0, 1)))) 
            stop("'tox' variable in data must be a vector of zeros (no toxicity) and ones (toxicity)")
        if (any(!(data$doseA %in% 1:I))) 
            stop(paste("'doseA' variable in data must contain the dose levels (1 to ",I,")", sep = ""))
        if (any(!(data$doseB %in% 1:J))) 
            stop(paste("'doseB' variable in data must contain the dose levels (1 to ",J,")", sep = ""))
	}            

  ## Check that non.admissible, if given, is the correct dimension
  if(!is.null(non.admissible)){
    if(dim(non.admissible)[1]!=I | dim(non.admissible)[2]!=J)
        stop(paste("Expecting non.admissible to be an ",I,"x",J," matrix",sep=""))
  }
  ## Check that non.admissible is populated by TRUE and FALSE
  if(!is.null(non.admissible)){
    if(!is.logical(non.admissible))
      stop("non.admissible must be a matrix of logicals")
  }
  
	## No. of dose combinations
	k=I*J

	## Monotonic matrices
	matrices<-monotonic.matrices(I,J)

	hsegments <- lapply(matrices , function(m) {
	  mp <- rbind( rep(0,J) , m , rep(1,J) )
	  mp[-1,] - mp[-(I+2),]
	})
	
	vsegments <- lapply(matrices , function(m) {
	  mp <- cbind( rep(0,I) , m , rep(1,I) )
	  mp[,-1] - mp[,-(J+2)]
	})
	
	# Number of doses being recommended for next cohort 
	# Currently PIPE only recommends 1 dose-combination for the next cohort
	doses<-1

	# Set up list of No. DLTs and Number Patients, recommended doses and number RPII doses per simulation
	r.sim<-n.sim<-list()
	rec.i.sim<-rec.j.sim<-array(NA,dim=c(S,N/c,doses))
	rec<-matrix(0,ncol=J,nrow=I)
	n.rpII<-vector()
	## Find a and b parameters from beta distribution with median = prior.med and prior strength given by prior.ss
	if(is.null(a) & is.null(b)){
		prior <- beta.med(prior.med,prior.ss)
		a<-prior$a
		b<-prior$b
	}

	# Set up more lists
	mat.list=uppermat.list=uppermat2.list=n.list=r.list=pi.theta.list=rpII.list=list()
	means <- cdfs <- list()
	mat.lik.list <- list()
	h.lik.list <- v.lik.list <- list()
	dom.list <- admis.list <- list() 

	
	
	for(s in 1:S){ # Loop over simulations
		# Initialise No DLTs and No. patients for each dose combination
		r=matrix(0,nrow=I,ncol=J)
		n=matrix(0,nrow=I,ncol=J)
				
		## Prior probability that each dose is less than or equal to theta
		p<-pbeta(theta,a,b)
		# If uppertox.constraint is specified find probability that each dose is less than or equal to uppertox
		if(!is.null(uppertox.constraint)){
			pconstraint<-pbeta(uppertox.constraint,a,b)
		} else {
			pconstraint<-NULL
		}
		# Initialise recommended doses for dimension i and j
		rec.i=rec.j=matrix(nrow=0,ncol=doses)

		## Where does the first cohort get dosed?
		create<-mtc.create(matrices,p,constraint,pconstraint,epsilon,admis,rec.i,rec.j,n,
		                   contour.select=contour.select,hsegments=hsegments,vsegments=vsegments,S,non.admissible,reweight,R,P)
		mat<-create$mat			
		
#		if(strategy=="ss" | strategy=="ss-random"){   ### NEEDED ANY MORE?
			## CHOOSE NEXT DOSE AS ONE LEAST EXPERIMENTED ON
			pi.theta<-1/(a+b)
#		} 
				
		# Using admissible doses, and prior sample sizes choose the next dose
		nxt<-mtc(create$dominant,create$admissible,strategy,rec.i,rec.j,pi.theta,mat,p,alternate,psmooth=create$mat.lik)

		
		# Store results if only one trial run
		if(S==1){
		  ## Posterior distributions based upon each contour
		  thetas <- seq(0.05,0.95,0.05)
		  ps<-lapply(thetas , pbeta , shape1=a , shape2=b)
		  
		  mtc.nums <- lapply(ps , function(p) unlist(lapply(matrices,function(l){prod((1-p)[l==1])*prod(p[l==0])})) )
		  mtc.liks <- lapply(mtc.nums , function(mtc.num) mtc.num/sum(mtc.num) )
		  
		  ## Posterior p(>MTC) for each dose combination
		  mat.liks <- lapply(mtc.liks , function(mtc.lik) sapply(1:length(mtc.lik),function(l){matrices[[l]]*mtc.lik[[l]]},simplify=F) )
		  cdfs[[1]] <- lapply(mat.liks , function(mat.lik) Reduce('+',mat.lik) )
		  means[[1]] <- 0.95 - 0.05*Reduce("+" , cdfs[[1]])
		  
			# This is most likely (modal) monotonic contour
			mat.list[[1]]=create$mat.mode
			mat.lik.list[[1]] <- create$mat.lik

			h.lik.list[[1]] <- create$h.lik
			v.lik.list[[1]] <- create$v.lik

			# This is monotonic matrix corresponding to upper toxicity constrain
			uppermat.list[[1]]=create$matupper
			# This is monotonic matrix corresponding to weighted Posterior p(>MTC) for each dose combination
			uppermat2.list[[1]]=create$matupper2
			dom.list[[1]] <- create$dominant
			admis.list[[1]] <- create$admissible
			n.list[[1]]=n
			r.list[[1]]=r
			pi.theta.list[[1]]=pi.theta
			
			
		}
	
		rec.i=nxt$rec.i
		rec.j=nxt$rec.j


		for(m in 1:(N/c)){ # Loop over cohorts
			if(doses==1){		
				if(!is.null(data) & nrow(data)>=m*c){ # If data is already specified then add it	
					if(length(unique(data$doseA[(m*c-c+1):(m*c)]))>1 | length(unique(data$doseB[(m*c-c+1):(m*c)]))>1){
						stop("Data given does not have all patients in the same cohort on the same dose combination")				
					} 
					for(pt in (m*c-c+1):(m*c)){
						r[data$doseA[pt],data$doseB[pt]]<-r[data$doseA[pt],data$doseB[pt]]+data$tox[pt]
						n[data$doseA[pt],data$doseB[pt]]<-n[data$doseA[pt],data$doseB[pt]]+1
					}
					rec.i[m,1]<-data$doseA[m*c]
					rec.j[m,1]<-data$doseB[m*c]
				} else if(mode=="alldlt" | mode=="nodlt") { # If nodlt or alldlt specified, add 0 or c respectively, else generate from binomial with true p(DLT)=pi
					r[rec.i[m,1],rec.j[m,1]]<-r[rec.i[m,1],rec.j[m,1]]+ifelse(mode=="nodlt",0,c)
					n[rec.i[m,1],rec.j[m,1]]<-n[rec.i[m,1],rec.j[m,1]]+c
				} else if(!is.null(pi)) {
					r[rec.i[m,1],rec.j[m,1]]<-r[rec.i[m,1],rec.j[m,1]]+rbinom(1,c,pi[rec.i[m,1],rec.j[m,1]])
					n[rec.i[m,1],rec.j[m,1]]<-n[rec.i[m,1],rec.j[m,1]]+c
				} else {
					break
				}
			} else { # This section is not currently used and corresponds to if more than one dose recommended after each cohort
				r[rec.i[m,1],rec.j[m,1]]<-r[rec.i[m,1],rec.j[m,1]]+ifelse(mode=="nodlt",0,ifelse(mode=="alldlt",c/2,rbinom(1,c/2,pi[rec.i[m,1],rec.j[m,1]])))
				n[rec.i[m,1],rec.j[m,1]]<-n[rec.i[m,1],rec.j[m,1]]+c/2
				r[rec.i[m,2],rec.j[m,2]]<-r[rec.i[m,2],rec.j[m,2]]+ifelse(mode=="nodlt",0,ifelse(mode=="alldlt",c/2,rbinom(1,c/2,pi[rec.i[m,2],rec.j[m,2]])))
				n[rec.i[m,2],rec.j[m,2]]<-n[rec.i[m,2],rec.j[m,2]]+c/2
			}
		
			# Calculate new posterior probailities
			p<-pbeta(theta,a+r,b+n-r)
			if(!is.null(uppertox.constraint)){
				pconstraint<-pbeta(uppertox.constraint,a+r,b+n-r)
			} else {
				pconstraint<-NULL
			}
			create<-mtc.create(matrices,p,constraint,pconstraint,epsilon,admis,rec.i,rec.j,n,contour.select=contour.select,hsegments=hsegments,vsegments=vsegments,S,non.admissible,reweight,R,P)
			mat<-create$mat		
			## CHOOSE NEXT DOSE AS ONE LEAST EXPERIMENTED ON
			## OR WEIGHTED RANDOMISATION BASED ON INVERSE SAMPLE SIZE (INCLUDING PRIOR SS)
			pi.theta<-1/(a+b+n)
			
			
  		if(S==1){ # For a single trial (i.e. no simulation) store information after each cohort
			  
			   ## Posterior distributions based upon each contour
			  thetas <- seq(0.05,0.95,0.05)
			  ps<-lapply(thetas , pbeta , shape1=a+r , shape2=b+n-r)
			  mtc.nums <- lapply(ps , function(p) unlist(lapply(matrices,function(l){prod((1-p)[l==1])*prod(p[l==0])})) )
			  mtc.liks <- lapply(mtc.nums , function(mtc.num) mtc.num/sum(mtc.num) )
			  
			  ## Posterior p(>MTC) for each dose combination
			  mat.liks <- lapply(mtc.liks , function(mtc.lik) sapply(1:length(mtc.lik),function(l){matrices[[l]]*mtc.lik[[l]]},simplify=F) )
			  cdfs[[m+1]] <- lapply(mat.liks , function(mat.lik) Reduce('+',mat.lik) )
			  means[[m+1]] <- 0.95 - 0.05*Reduce("+" , cdfs[[m+1]])
			  
			  ## Most likely (modal) contour
			  mat.list[[m+1]]=create$mat.mode
			  mat.lik.list[[m+1]] <- create$mat.lik
			
				h.lik.list[[m+1]] <- create$h.lik
				v.lik.list[[m+1]] <- create$v.lik     

				uppermat.list[[m+1]]=create$matupper
				uppermat2.list[[m+1]]=create$matupper2
				dom.list[[m+1]] <- create$dominant
				admis.list[[m+1]] <- create$admissible
				n.list[[m+1]]=n
				r.list[[m+1]]=r
				pi.theta.list[[m+1]]=pi.theta
				
			}
			
			## IF NO DOSES ARE ADMISSIBLE THEN STOP THE TRIAL FOR SAFETY
			if(all(!create$admissible)){
				rec.i<-rbind(rec.i,matrix(0,nrow=N/c+1-m,ncol=doses))
				rec.j<-rbind(rec.j,matrix(0,nrow=N/c+1-m,ncol=doses))
				break
			}
			if(!is.null(stop)){
				## If lowest dose combination has posterior probability of being greater than the TTL of > stop then stop the trial
				if(1-p[1,1]>stop){
					rec.i<-rbind(rec.i,matrix(0,nrow=N/c+1-m,ncol=doses))
					rec.j<-rbind(rec.j,matrix(0,nrow=N/c+1-m,ncol=doses))
					break
				}
			}
		
			# Find next dose
			nxt<-mtc(create$dominant,create$admissible,strategy,rec.i,rec.j,pi.theta,mat,p,alternate,psmooth=create$mat.lik)
			rec.i=nxt$rec.i
			rec.j=nxt$rec.j
		}
		r.sim[[s]]<-r
		n.sim[[s]]<-n
		## Recommendations made for each cohort throughout the trial
		rec.i.sim[s,,]<-rec.i[-(m+1),]
		rec.j.sim[s,,]<-rec.j[-(m+1),]
		
		## Recommended PII dose combinations
		## THESE ARE DOSE COMBINATIONS THAT HAVE BEEN EXPERIMENTED ON, ARE CLOSEST TO ESTIMATED MTC_THETA
		## AND ARE LOWER THAN UPPER CONSTRAINT CONTOUR, OR p(dose>MTC)< epsilon
		## CAN ONLY RECOMMEND PII DOSES IF TRIAL IS NOT STOPPED EARLY
		if(any(create$admissible)){
			create.rpII<-mtc.create(matrices,p,constraint="none",pconstraint=pconstraint,epsilon,admis="closest",rec.i,rec.j,n,
			                        contour.select=contour.select,hsegments=hsegments,vsegments=vsegments,S,non.admissible,reweight=reweight,R=R,P=P)
			rpIIs<-create.rpII$dominant & create.rpII$mat==0 & n!=0 & create$matupper==0 & create$matupper2==0 
			rpII.i<-row(mat)[rpIIs]
			rpII.j<-col(mat)[rpIIs]
			for(i in 1:length(rpII.i)){
				rec[rpII.i[i],rpII.j[i]]<-rec[rpII.i[i],rpII.j[i]]+1
			}
			n.rpII[s]<-length(rpII.i)
		} else {
		## If trial has stopped early
			rpII.i<-rpII.j<-numeric(0)
			n.rpII[s]<-0
		}
		rpII<-rbind(rpII.i,rpII.j)
		rownames(rpII)<-c("rpII.A","rpII.B")
		rpII.list[[s]]<-rpII

		if(S>1) cat(s,"\n")
	}

	exp<-Reduce('+',n.sim)/sum(Reduce('+',n.sim))
	no.not.treated<-N*S-sum(Reduce('+',n.sim))
	dlts<-sapply(1:S,function(l){sum(r.sim[[l]])/sum(n.sim[[l]])})
	results<-list(r.sim=r.sim,n.sim=n.sim,rec.i.sim=rec.i.sim,rec.j.sim=rec.j.sim,exp=exp,rec=rec/sum(rec),dlts=dlts,mat.list=mat.list,uppermat.list=uppermat.list,uppermat2.list=uppermat2.list,r.list=r.list,n.list=n.list,n.rpII=n.rpII,
	              no.not.treated=no.not.treated,pi=pi,theta=theta,a=a, b=b,
	              pi.theta.list=pi.theta.list, cdfs=cdfs, means=means, mat.lik.list=mat.lik.list , 
	              h.lik.list=h.lik.list , v.lik.list=v.lik.list , dom.list=dom.list , admis.list=admis.list,
	              rpII.list=rpII.list)
	if(S>1){
		class(results)<-"pipe.sim"
	} else {
		class(results)<-"pipe"
	}
	return(results)
}


## Function that returns all monotonic matrices of dimension IxJ
monotonic.matrices<-function(I,J){
	comb.col<-combinations(2, J, c(0,1), repeats.allowed=TRUE)
	n.com1<-dim(comb.col)[1]
	comb.row<-combinations(n.com1,I,repeats.allowed=TRUE)
	n.com2<-dim(comb.row)[1]
	matrices<-sapply(1:n.com2,function(i){comb.col[comb.row[i,],]},simplify=F)
	return(matrices)
}

## Obtain doses closest to MTC
closest<-function(mat){
	I<-nrow(mat)
	J<-ncol(mat)
	dominantu<-mat==1 & rbind(0,mat[-I,]) %in% c(0,2) & cbind(0,mat[,-J]) %in% c(0,2)
	dominantl<-mat==0 & rbind(mat[-1,],1) %in% c(1,2) & cbind(mat[,-1],1) %in% c(1,2)
	dominant<-dominantl | dominantu
	dominant
}


## Create admissible dose matrix
mtc.create<-function(matrices,p,constraint,pconstraint,epsilon,admis,rec.i,rec.j,n,contour.select=contour.select,hsegments=hsegments,vsegments=vsegments,S,non.admissible,reweight,R,P){
	m<-dim(rec.i)[1]

	## Assess all possible MTC contours to find most likely
	mtc.num<-unlist(lapply(matrices,function(l){prod((1-p)[l==1])*prod(p[l==0])}))
	mtc.lik<-mtc.num/sum(mtc.num)
	
  ## Re-weight contours
  if(reweight){
    I<-dim(matrices[[1]])[1]
	  J<-dim(matrices[[1]])[2]
    Q<-unlist(lapply(matrices,sum)) ## Using Q instead of S as S == Total no. of sims
	  ## Weighting  controlled by R and P (F has been renamed R - as F represents FALSE in R language)
	  newweight<- 1 + R*abs(2/(I*J)*Q - 1 )^P
	  mtc.lik<-newweight*mtc.lik/sum(newweight*mtc.lik)
  }

	## Should we choose contour based on posterior median or mode?
	if(contour.select=="median" | S==1){	 ## Only evaluate this code if median contour is specified or if not a simulation (i.e. S==1)
	  h.lik <- Reduce( "+" , lapply(1:length(mtc.lik) , function(i) mtc.lik[[i]]*hsegments[[i]]) )
	  v.lik <- Reduce( "+" , lapply(1:length(mtc.lik) , function(i) mtc.lik[[i]]*vsegments[[i]]) )     
	  nv <- ncol(h.lik)
	  nh <- nrow(v.lik)
	  h.med <- apply(h.lik , 2 , w.median , x=0:nh)
	  v.med <- apply(v.lik , 1 , w.median , x=0:nv)	 
	  
	  mat.hmed <- do.call(cbind , lapply(h.med , function(i) rep(0:1,c(i,nh-i))))
	  mat.vmed <- do.call(rbind , lapply(v.med , function(i) rep(0:1,c(i,nv-i))))
	  
	  mtc.hmed <- which(sapply( matrices , function(m) sum(abs(m-mat.hmed))==0))
	  mtc.vmed <- which(sapply( matrices , function(m) sum(abs(m-mat.vmed))==0))
	  
	  if(mtc.lik[mtc.hmed]>=mtc.lik[mtc.vmed]) mat.med <- mat.hmed
	  else mat.med <- mat.vmed  
	  
	  mtc.mode<-which.max(mtc.lik)
	  mat.mode<-matrices[[mtc.mode]]
	  
	} else {
	  mtc.mode<-which.max(mtc.lik)
	  mat.mode<-matrices[[mtc.mode]]
	  h.lik<-v.lik<-mat.med<-NULL  
  }

	
	if(contour.select=="median") mat <- mat.med
	else mat <- mat.mode 

	I<-dim(mat)[1]
	J<-dim(mat)[2]

	## Find upper toxicity constraint contour if uppertox.constraint is used
	if(!is.null(pconstraint)){
		upper.lik<-which.max(unlist(lapply(matrices,function(l){prod((1-pconstraint)[l==1])*prod(pconstraint[l==0])})))
		matupper<-matrices[[upper.lik]]
	} else {
		matupper<-matrix(0,nrow=I,ncol=J)
	}

	## Find doses that do not satisfy constraint if weightedMTC.constraint is used
	## Posterior p(>MTC) for each dose combination
	mat.lik <- Reduce('+', sapply(1:length(mtc.lik),function(l){matrices[[l]]*mtc.lik[[l]]},simplify=F) )
	if(!is.null(epsilon)){
	  weight.pMTC<-mat.lik
		matupper2<-weight.pMTC>=epsilon
	} else {
	  weight.pMTC<-NULL
		matupper2<-matrix(0,nrow=I,ncol=J)
	}
	
	if(grepl("neighbouring",constraint)){
		## IF NEIGHBOURING CONSTRAINT AND MORE THAN ONE DOSE RECOMMENDATION PER COHORT THEN USE UNION OF BOTH ADMISSIBLE REGIONS
		if(dim(rec.i)[2]>1){
			admissible1<- row(mat)<=max(rec.i[m,1],0)+1 & col(mat)<=max(rec.j[m,1],0)+1 & row(mat)>=max(rec.i[m,1],0)-1 & col(mat)>=max(rec.j[m,1],0)-1
			admissible2<- row(mat)<=max(rec.i[m,2],0)+1 & col(mat)<=max(rec.j[m,2],0)+1 & row(mat)>=max(rec.i[m,2],0)-1 & col(mat)>=max(rec.j[m,2],0)-1
			admissible<- admissible1 | admissible2
		} else {
			admissible<- row(mat)<=max(rec.i[m,1],0)+1 & col(mat)<=max(rec.j[m,1],0)+1 & row(mat)>=max(rec.i[m,1],0)-1 & col(mat)>=max(rec.j[m,1],0)-1
		}
	} else if(grepl("no.dose.skip",constraint)){
		## IF NO.DOSE.SKIP CONSTRAINT AND MORE THAN ONE DOSE RECOMMENDATION PER COHORT THEN USE UNION OF BOTH ADMISSIBLE REGIONS
		if(dim(rec.i)[2]>1){
			admissible1<- row(mat)<=max(rec.i[,1],0)+1 & col(mat)<=max(rec.j[,1],0)+1
			admissible2<- row(mat)<=max(rec.i[,2],0)+1 & col(mat)<=max(rec.j[,2],0)+1
			admissible<- admissible1 | admissible2
		} else {
			admissible<- row(mat)<=max(rec.i[,1],0)+1 & col(mat)<=max(rec.j[,1],0)+1 
		}
	} else {
		## IF NO NEIGHBOURING CONSTRAINT THEN ALL DOSE COMBINATIONS ARE ADMISSIBLE	
		admissible<- matrix(TRUE,nrow=I,ncol=J)		
	}
	
	## IF A NODIAG CONSTRAINT IS ADDITIONALLY SPECIFIED
	if(grepl("nodiag",constraint)){
		if(dim(rec.i)[1]>=1){
			if(rec.i[m,1]!=I & rec.j[m,1]!=J){
				admissible[rec.i[m,1]+1,rec.j[m,1]+1]<-FALSE
			}
		}
	}

	## IF a NON-ADMISSIBLE RANGE IS ADDITIONALLY SELECTED
	if(!is.null(non.admissible)){
	  admissible<-admissible & !non.admissible
	}
	
	
	if(!is.null(pconstraint) | !is.null(epsilon)){
		admissible<-admissible & matupper==0 & matupper2==0
		## IF THERE ARE NO ADMISSIBLE DOSES LEFT (THAT IS ALL NEIGHBOURING DOSES ARE NOW UNSAFE)
		## CHOOSE CLOSEST DOSE THAT IS SAFE (TO FIRST COHORT DOSE)
		if(all(!admissible)){
			test<-abs(rec.i[m,1]-row(mat))+abs(rec.j[m,1]-col(mat))
			admissible<- test==min(c(test[matupper==0 & matupper2==0],-Inf)) & matupper==0 & matupper2==0
		}
	}
	
	separate<-FALSE
	if(admis=="adjacent"){
  		## ANY DOSE COMBINATION ADJACENT TO THE MTC IS ADMISSIBLE
		if(I<2 | J<2) stop("Admissible doses can only be calculated when both drugs have more than one level")
		admat<-mat
		dominantu<-admat==1 & (rbind(0,admat[-I,])==0 | cbind(0,admat[,-J])==0 | rbind(0,cbind(0,admat[,-J])[-I,])==0)
		dominantl<-admat==0 & (rbind(admat[-1,],1)==1 | cbind(admat[,-1],1)==1 | rbind(cbind(admat[,-1],1)[-1,],1)==1)
		dominant<-dominantl | dominantu
		## If dominant and admissible regions are separate choose the "closest" dose in the admissible region
		if(!any(dominant & admissible)){
			separate<-TRUE
		}
	} 
	if(admis=="closest" | separate==TRUE){
		## ONLY DOSE COMBINATIONS CLOSEST TO THE MTC ARE ADMISSIBLE
		admat <- mat
		## SET ALL DOSES OUTSIDE ADMISSIBLE RANGE TO 2 AND ALLOW ANY TOUCHING CONSTRAINT TO BE DOMINANT (IF SATISFY OTHER CRITERIA)
		admat[admissible==FALSE]=2
							
		dominant<-closest(admat)
	}
	return(list(dominant=dominant,admissible=admissible,mat=mat,mat.mode=mat.mode,mat.med=mat.med,matupper=matupper,matupper2=matupper2,weight.pMTC=weight.pMTC,mat.lik=mat.lik,h.lik=h.lik,v.lik=v.lik))
}

## Choose next dose from admissible doses that use either smallest sample size or weighted randomisation of sample size
mtc<-function(dominant,admissible,strategy,rec.i,rec.j,pi.theta,mat,p,alternate,psmooth){

	m<-dim(rec.i)[1]
	I<-dim(pi.theta)[1]
	J<-dim(pi.theta)[2]
	k<-I*J

	if(alternate==T){
	  ## If more than one admissible & dominant dose then go below MTC if last dose is above, or vice-versa
	  ## If this is the first cohort then go below MTC
	  if(sum(dominant & admissible)>1){
	    if(m==0){
	      if(any(mat[dominant & admissible]==0)) dominant[mat==1]<-FALSE
	    } else {
	      if(mat[rec.i[[m]],rec.j[[m]]]==1){
	        if(any(mat[dominant & admissible]==0)) dominant[mat==1]<-FALSE
	      } else{
	        if(any(mat[dominant & admissible]==1)) dominant[mat==0]<-FALSE
	      }		
	    }
	  }
	}

	## Strategy "ss": Select the dominant dose with smallest sample size
	## If there are more than two dose combinations that are equal choose from them at random
	if(strategy=="ss"){
	  test<- pi.theta==max(pi.theta[dominant & admissible]) & dominant & admissible
	  ## If still more than one dose comb. then choose one at random
	  chosen=ifelse(sum(test)>1,sample(sum(test),1),1)
	  rec.i<-rbind(rec.i,row(pi.theta)[test][chosen])
	  rec.j<-rbind(rec.j,col(pi.theta)[test][chosen])
	} 
	
	## weighted_mtc strategy not currently implemented
	else if(strategy=="ss-random"  | strategy=="weighted_mtc"){
	  ## Strategy "ss-random": Select the dominant dose with probability weighted by inverse sample size or weighted by the (weighted) probability of being a closest dose to the true MTC
	  pi.theta[!(dominant & admissible)]=0
	  chosen=sample(k,1,prob=pi.theta)
	  rec.i<-rbind(rec.i,row(pi.theta)[chosen])
	  rec.j<-rbind(rec.j,col(pi.theta)[chosen])
	} 
	
	## p strategy not currently implemented
	else if(strategy=="p"){
	  ## Strategy "p": Select the dominant dose with posterior probability of being less than p closest to 0.5 (i.e. most uncertain)
	  p[!(dominant & admissible)]=2
	  test<- abs(p-0.5)==min(abs(p-0.5)) & dominant & admissible
	  ## If still more than one dose comb. then choose one at random
	  chosen=ifelse(sum(test)>1,sample(sum(test),1),1)
	  rec.i<-rbind(rec.i,row(pi.theta)[test][chosen])
	  rec.j<-rbind(rec.j,col(pi.theta)[test][chosen])
	} 
	
	## Strategy "psmooth": Select the dominant dose with posterior probability using the smoothed probabilities of being less than p that is closest to 0.5 (i.e. most uncertain)	## 
	else if(strategy=="psmooth"){
	  psmooth[!(dominant & admissible)]=2
	  test<- abs(psmooth-0.5)==min(abs(psmooth-0.5)) & dominant & admissible
	  ## If more than one dose comb. choose one at random
	  chosen=ifelse(sum(test)>1,sample(sum(test),1),1)
	  rec.i<-rbind(rec.i,row(pi.theta)[test][chosen])
	  rec.j<-rbind(rec.j,col(pi.theta)[test][chosen])
	}
	
	## IS THIS STRATEGY NEEDED?
	else if(strategy=="ssp") {  
	  test<- pi.theta==max(pi.theta[dominant & admissible]) & dominant & admissible
	  p[!(test)]=2  
	  test <- abs(p-0.5)==min(abs(p-0.5)) 
	  chosen=ifelse(sum(test)>1,sample(sum(test),1),1)
	  rec.i<-rbind(rec.i,row(pi.theta)[test][chosen])
	  rec.j<-rbind(rec.j,col(pi.theta)[test][chosen])
	}
	
	return(list(rec.i=rec.i,rec.j=rec.j))	
}


## Obtain a and b parameters for beta prior from median and sample size using numerical optimisation
beta.med<-function(prior.med,prior.ss){
  if(any(prior.med==0 | prior.med==1)) stop("`prior.med' must be greater than 0 and less than 1")
	betaprior1 = function(K, x, p) {
        m.lo = 0
        m.hi = 1
        flag = 0
        while (flag == 0) {
            m0 = (m.lo + m.hi)/2
            p0 = pbeta(x, K * m0, K * (1 - m0))
            if (p0 < p) 
                m.hi = m0
            else m.lo = m0
            if (abs(p0 - p) < 1e-04) 
             flag = 1
   		}
		return(m0)
	}
	a.med<-b.med<-matrix(NA,nrow=dim(prior.med)[1],ncol=dim(prior.med)[2])
	for(i in 1:dim(prior.med)[1]){
		for(j in 1:dim(prior.med)[2]){
			a.med[i,j]<-prior.ss[i,j]*betaprior1(prior.ss[i,j],prior.med[i,j],0.5)
			b.med[i,j]<-prior.ss[i,j]-a.med[i,j]
		}
	}
	return(list(a=a.med,b=b.med))
}


## Function to plot dose-escalation steps for a PIPE design
## Arguments:
## x: a PIPE object as obtained from pipe
## type: What data should be shown on the figure. Options are
## 			"b" (default): Show both the number of DLTs (numerator) and number of patients recruited (denominator) at each dose combination
## 			"r": Only show number of DLTs (numerator) at each dose combination
## 			"n": Only show number of patients recruited (denominator) at each dose combination
## pi: True p(DLT) matrix
## theta: Target toxicity level
## epsilon.line: Plot safety constraint line formed by a weighted average of all possible MTCs
## uppertox.constraint.line: Plot safety constraint line based on most likely MTC at upper toxicity constraint
plot.pipe<-function(x,type="b",pi=x$pi,theta=x$theta,epsilon.line=TRUE,uppertox.constraint.line=FALSE,add.empirical.data=FALSE,...){
	mat<-x$mat.list
	c<-length(mat)
	I<-dim(mat[[1]])[1]
	J<-dim(mat[[1]])[2]
	cohort.size=sum(x$n.sim[[1]])/(c-1)
	xlevels=rep(1:(I+1)-0.5,c)
	ylevels=c(sapply(1:c,function(k){c(apply(mat[[k]],1,function(i){min(which(i==1),J+1)-0.5}),0.5)}))
	cohort=rep(1:c,each=(I+1))
	if(!is.null(pi)){
		mat.true<-pi>theta
		x.true<-1:(I+1)-0.5
		y.true<-c(apply(mat.true,1,function(i){min(which(i==1),J+1)-0.5}),0.5)
	}
	## Estimated MTC
	df<-data.frame(x=xlevels,y=ylevels,cohort=cohort)
	ncol<-2 ## Minimum number of coloured tiles to be shown in plot
	## Data to be shown
	if(type=="n"){
		df2<-data.frame(x=rep(rep(1:I,J),c-1),y=rep(rep(1:J,each=I),c-1),z=factor(unlist(x$n.list[-1]),levels=unique(sort(unlist(x$n.list[-1]),decreasing=T))),cohort=rep(1:(c-1),each=I*J))
		ncol<-length(levels(df2$z))
	} else if(type=="r"){
		df2b<-data.frame(x=rep(rep(1:I,J),c-1),y=rep(rep(1:J,each=I),c-1),z=unlist(x$r.list[-1]),cohort=rep(1:(c-1),each=I*J))
		df2b<-df2b[df2b$z!=0,]
		df2b$z<-factor(df2b$z,levels=unique(sort(df2b$z,decreasing=T)))
	} else if(type=="b"){
		df2<-data.frame(x=rep(rep(1:I,J),c-1),y=rep(rep(1:J,each=I),c-1),z=factor(unlist(x$n.list[-1]),levels=unique(sort(unlist(x$n.list[-1]),decreasing=T))),cohort=rep(1:(c-1),each=I*J))
		ncol<-length(levels(df2$z))
		df2b<-data.frame(x=rep(rep(1:I,J),c-1),y=rep(rep(1:J,each=I),c-1),z=unlist(x$r.list[-1]),cohort=rep(1:(c-1),each=I*J))
		df2b<-df2b[df2b$z!=0,]
		df2b$z<-factor(df2b$z,levels=unique(sort(df2b$z,decreasing=T)))
	}
	## True MTC 
	if(!is.null(pi)){
		df3<-data.frame(x=x.true,y=y.true,cohort=rep(c,I+1))
	}
	## Recommended PII doses
	df4<-data.frame(x=rep(1:I,J),y=rep(1:J,each=I),z=factor(c(cohort.size*as.numeric(x$rec!=0)),levels=c(cohort.size,0)),cohort=rep(c,I*J))
	
	## Empirical tox. probs
	if(add.empirical.data){
		a.post<-x$a+x$r.sim[[1]]
		b.post<-x$b+(x$n.sim[[1]]-x$r.sim[[1]])
		emp.medians<-round(qbeta(0.5,a.post,b.post),2)
		emp.lower2.5<-round(qbeta(0.025,a.post,b.post),2)
		emp.upper2.5<-round(qbeta(0.975,a.post,b.post),2)
		df6<-data.frame(x=rep(1:I,J),y=rep(1:J,each=I),median=100*c(emp.medians),lower=100*c(emp.lower2.5),upper=100*c(emp.upper2.5),cohort=rep(c,I*J))
		df6$ci<-paste("(",df6$lower,", ",df6$upper,")",sep="")
	}
	v1<-ggplot()+
		geom_step(aes(x=x,y=y),data=df,size=2)+facet_wrap(~cohort)+
		xlab("Drug A level")+ylab("Drug B level")
	v2<-ggplot()+
		geom_step(aes(x=x,y=y),data=df[df$cohort==c,],size=2)+
		xlab("Drug A level")+ylab("Drug B level")
	if(!is.null(pi)){
		v1<-v1+geom_step(aes(x=x,y=y),data=df3,size=1.5,colour="green",linetype=4)
		v2<-v2+geom_step(aes(x=x,y=y),data=df3,size=1.5,colour="green",linetype=4)
	}
	if(type=="n" | type=="b"){
		v1<-v1+geom_tile(aes(x=x,y=y,fill = z),alpha=0.5,data=df2)
	}
	if(type=="r" | type=="b"){
		if(dim(df2b)[1]>0){
			v1<-v1+geom_point(aes(x=x,y=y,shape = z),size=2,data=df2b)+scale_shape(name="Number DLTs")	
		}
	}
	v1<-v1+geom_tile(aes(x=x,y=y,fill = z),alpha=0.5,data=df4)+scale_fill_manual(name="Number pts.",values=c(rainbow(ncol-1),"#FFFFFFFF"))
	if(length(x$uppermat.list)>0 & uppertox.constraint.line){
		uppermat<-x$uppermat.list
		x.high<-rep(1:(I+1)-0.5,c)
		y.high<-c(sapply(1:c,function(k){c(apply(uppermat[[k]],1,function(i){min(which(i==1),J+1)-0.5}),0.5)}))
		df5<-data.frame(x=x.high,y=y.high,cohort=cohort)
		v1<-v1+geom_step(aes(x=x,y=y),data=df5,size=1,colour="red")
	}
	if(length(x$uppermat2.list)>0 & epsilon.line){
		uppermat2<-x$uppermat2.list
		x.high<-rep(1:(I+1)-0.5,c)
		y.high<-c(sapply(1:c,function(k){c(apply(uppermat2[[k]],1,function(i){min(which(i==1),J+1)-0.5}),0.5)}))
		df5<-data.frame(x=x.high,y=y.high,cohort=cohort)
		v1<-v1+geom_step(aes(x=x,y=y),data=df5,size=1,colour="red4",linetype=4)
	}
	print(v1)
	if(add.empirical.data){
		v2<-v2+geom_text(aes(x=x,y=y+0.1,label=median),data=df6)+geom_text(aes(x=x,y=y-0.1,label=ci),data=df6)
		dev.new()
		print(v2)
	}
	#if(length(x$plot.density)>0){
	#	for(m in 1:length(x$plot.density)){
	#			print(x$plot.density[[m]])
	#	}
	#}
}

## Function to plot dose-escalation steps for a PIPE design
## Arguments:
## x: a PIPE object as obtained from pipe
## pi: True p(DLT) matrix
## theta: Target toxicity level
## plot: What operating characteristics should be plotted? Options are:
##		"exp": Experimentation proportions heat map
##		"rec": Recommendation proportions heat map
plot.pipe.sim<-function(x,pi=x$pi,theta=x$theta,plot="both",...){
	exp<-x$exp
	rec<-x$rec
	I<-dim(x$n.sim[[1]])[1]
	J<-dim(x$n.sim[[1]])[2]

	mat.true<-pi>theta
	x.true<-1:(I+1)-0.5
	y.true<-c(apply(mat.true,1,function(i){min(which(i==1),J+1)-0.5}),0.5)
	
	s<-length(x$n.sim)
	
	if(plot=="exp"){
		# Experimentation proportions plot
		df<-data.frame(x=rep(1:I,J),y=rep(1:J,each=I),z=100*c(exp))
		df2<-data.frame(x=x.true,y=y.true)
		v1<-ggplot()+geom_tile(aes(x=x,y=y,fill = z),data=df)+scale_fill_gradient(name="Experimentation percentages",low="white",high="red")+
			geom_step(aes(x=x,y=y),data=df2,size=1.5,colour="green",linetype=4)+
			xlab("Drug A level")+ylab("Drug B level")
		print(v1)
	}

	if(plot=="rec"){
		# Recommendation proportions
		df<-data.frame(x=rep(1:I,J),y=rep(1:J,each=I),z=100*c(rec))
		df2<-data.frame(x=x.true,y=y.true)
		v1<-ggplot()+geom_tile(aes(x=x,y=y,fill = z),data=df)+scale_fill_gradient(name="Recommendation percentages",low="white",high="red")+
			geom_step(aes(x=x,y=y),data=df2,size=1.5,colour="green",linetype=4)+
			xlab("Drug A level")+ylab("Drug B level")
		print(v1)
	}

	if(plot=="both"){
		# Experimentation and Recommendation proportions plot
		df.exp<-data.frame(x=rep(1:I,J),y=rep(1:J,each=I),z=100*c(exp),type="Experimentation")
		df.rec<-data.frame(x=rep(1:I,J),y=rep(1:J,each=I),z=100*c(rec),type="Recommendation")
		df<-rbind(df.exp,df.rec)
		df2<-data.frame(x=x.true,y=y.true)
		v1<-ggplot()+geom_tile(aes(x=x,y=y,fill = z),data=df)+facet_grid(~type)+scale_fill_gradient(name="Percent",low="white",high="red")+
			geom_step(aes(x=x,y=y),data=df2,size=1.5,colour="green",linetype=4)+
			xlab("Drug A level")+ylab("Drug B level")
		print(v1)
	}
}


print.pipe<-function(x,...){
	I=dim(x$r.sim[[1]])[1]
	J=dim(x$r.sim[[1]])[2]
	n<-x$n.sim[[1]]
	r<-x$r.sim[[1]]
	mat<-x$mat.list[[length(x$mat.list)]]
	matupper<-x$uppermat.list[[length(x$uppermat.list)]]
	matupper2<-x$uppermat2.list[[length(x$uppermat2.list)]]
	tab1<-t(n)[J:1,]
	rownames(tab1)<-paste("Level",J:1)
	colnames(tab1)<-paste("Level",1:I)
	names(dimnames(tab1))<-c("Drug B","Drug A")
	tab2<-t(r)[J:1,]
	rownames(tab2)<-paste("Level",J:1)
	colnames(tab2)<-paste("Level",1:I)
	names(dimnames(tab2))<-c("Drug B","Drug A")
	tab3<-t(mat)[J:1,]
	rownames(tab3)<-paste("Level",J:1)
	colnames(tab3)<-paste("Level",1:I)
	names(dimnames(tab3))<-c("Drug B","Drug A")
	tab4<-t((matupper | matupper2)*1)[J:1,]
	rownames(tab4)<-paste("Level",J:1)
	colnames(tab4)<-paste("Level",1:I)
	names(dimnames(tab4))<-c("Drug B","Drug A")
	cat("\n Number of patients dosed:\n")
	print(tab1)
	cat("\n Toxicities observed:\n")
	print(tab2)
	if(length(x$r.list)==length(x$rec.i.sim)){
		cat("\n Next recommended dose level: \n Dose A: ",x$rec.i.sim[1,length(x$rec.i.sim),1],"\n Dose B: ",x$rec.j.sim[1,length(x$rec.j.sim),1],"\n")
	}
	cat("\n MTC:\n")
	print(tab3)
	cat("\n Upper toxicity constraint (1 indicates doses not allowed):\n")
	print(tab4)
}

###### Function to print experimentation and recommendation percentages from a simulated PIPE design
## Arguments:
## x: a PIPE object as obtained from pipe
## pi: True p(DLT) matrix
## cut.points: cut points on the true DLT range to present the operating characteristics
## digits: Number of decimal places to be used for reporting
## print: (defaul=TRUE). Should the output be printed?
print.pipe.sim<-function(x,pi=x$pi,cut.points=c(0,15,25,35,45.1,100)/100,digits=1,print=TRUE,...){
	exp<-x$exp
	rec<-x$rec

	cuts<-cut(pi,cut.points,right=F)
	# Experimentation proportions table
	exp.table<-sapply(levels(cuts),function(i){sum(exp[cuts==i])})

	# Recommendation proportions table
	rec.table<-sapply(levels(cuts),function(i){sum(rec[cuts==i])})

	if(print){
		cat("\n Experimentation percentages by true toxicity: \n")
		print(round(100*exp.table,digits))

		cat("\n Recommendation percentages by true toxicity: \n")
		print(round(100*rec.table,digits))		
	}
	return(list(exp.table=exp.table,rec.table=rec.table))
}



plothists <- function(ds) {
  
  temp <- ds$cdfs[[length(ds$cdfs)]]
  
  temp2 <- lapply(0:length(temp) , function(i) {
    if(i==0) 1-temp[[1]]
    else if(i==length(temp)) temp[[length(temp)]]
    else temp[[i]]-temp[[i+1]]  
  })
  
  temp3 <- array(NA , dim=c(1+length(temp),nrow(temp2[[1]]),ncol(temp2[[1]])) )
  for(i in 1:nrow(temp2[[1]])) {
    for(j in 1:ncol(temp2[[1]])) {
      temp3[,i,j] <- sapply(temp2 , function(m) m[i,j])
    }
  }
  
  #  apply(temp3 , c(2,3) , sum)
  
  par(mfrow=dim(temp3)[3:2] , mar=rep(0.5,4),oma=c(4,4,1,1)) ## Plot drug A along x-axis and drug B along y-axis
  for(j in dim(temp3)[3]:1) { ## Levels of Drug B decreasing
    for(i in 1:dim(temp3)[2]) { ## Levels of Drug A increasing
      barplot(temp3[,i,j] , ylim=c(0,0.2) , border=NA , space=0 , width=rep(0.05,21) , axes=F , col=rep(c("GREEN","RED"),c(ds$theta/0.05 , (1-ds$theta)/0.05)+1e-15 ))
      box()
      if(length(ds$cdfs)==1) lines(seq(0,1,0.01) , 0.2*dbeta(seq(0,1,0.01) , ds$a[i,j] , ds$b[i,j]))
      else lines(seq(0,1,0.01) , 
                 0.2*dbeta(seq(0,1,0.01) , ds$a[i,j]+ds$r.sim[[1]][i,j] , ds$b[i,j]+ds$n.sim[[1]][i,j]-ds$r.sim[[1]][i,j]))
      if(j==1) mtext(i,side=1,line=1)
      if(i==1) mtext(j,side=2,line=1)
    }
  }
  mtext("Drug A level", side=1, outer=TRUE,line=2)
  mtext("Drug B level", side=2, outer=TRUE,line=2)
  par(mfrow=c(1,1))
  
}


plotsegs <- function(object) {
  
  h.lik<-object$h.lik[[length(object$h.lik)]]
  v.lik<-object$v.lik[[length(object$v.lik)]]  
  modal<-object$mat.list[[length(object$mat.list)]] 
  dom<-object$dom.list[[length(object$dom.list)]] 
  admis<-object$admis.list[[length(object$admis.list)]]
  n<-object$n.list[[length(object$n.list)]] 
  r<-object$r.list[[length(object$r.list)]] 
  mat.lik<-object$mat.lik.list[[length(object$mat.lik.list)]] 
  means<-object$means[[length(object$means)]]  
  rec.i<-object$rec.i.sim 
  rec.j<-object$rec.j.sim

  opar <- par(mar=rep(0.5,4))
  nv <- ncol(h.lik)
  nh <- nrow(v.lik)
  plot(c(0 , nh) , c(0 , nv), type="n" , axes=F , xlab="" , ylab="")  
  box()
  #abline(h=0:(-nh) , col="grey")  
  #abline(v=0:nv , col="grey")
  
  ## Previous design point
  if(!is.null(rec.i)) {
    allx <- rec.i[1,,1]
    ally <- rec.j[1,,1]
    if(length(allx)>1) {
      x <- allx[length(allx)-1]
      y <- ally[length(ally)-1]
      ## polygon(c(y-1,y,y,y-1) , c(-x,-x,1-x,1-x) , col="grey90")
      polygon(c(x,x,x-1,x-1),c(y-1,y,y,y-1), col="grey90") 
    }  
    x <- allx[length(allx)]
    y <- ally[length(ally)]
    ##polygon(c(y-1,y,y,y-1) , c(-x,-x,1-x,1-x) , col="grey")
    polygon(c(x,x,x-1,x-1),c(y-1,y,y,y-1), col="grey")    
  }
  
  for(h in 0:nh) {
    for(v in 0:nv) {
      ##if(v<nv) segments(v , -h , v+1 , -h , lwd=25*h.lik[h+1,v+1])
      if(v<nv) segments(h,   v , h , v+1 , lwd=25*h.lik[h+1,v+1])
      ##if(h<nh) segments(v , -h , v , -(h+1) , lwd=25*v.lik[h+1,v+1])
      if(h<nh) segments(h,   v , h+1 , v , lwd=25*v.lik[h+1,v+1])
    }
  }
  
  h.med <- apply(h.lik , 2 , w.median , x=0:nh)
  v.med <- apply(v.lik , 1 , w.median , x=0:nv)
  ##segments(0:(nv-1) , -h.med , 1:nv , -h.med , col="RED" , lwd=2)
  segments(h.med,0:(nv-1) , h.med , 1:nv , col="RED" , lwd=2)  
  ##segments(v.med , -(0:(nh-1)) , v.med , -(1:nh) , col="RED" , lwd=2)
  segments((0:(nh-1)), v.med , (1:nh), v.med , col="RED" , lwd=2)  
  if(!is.null(modal)) {
    modalplus <- cbind(0 , rbind(0,modal,1) , 1)
    h.mod <- apply( modal==0 , 2 , function(x) { w<-which(x);if(length(w)==0) w<-0; max(w) } ) 
    v.mod <- apply( modal==0 , 1 , function(x) { w<-which(x);if(length(w)==0) w<-0; max(w) } )     
  ##  segments(0:(nv-1) , -h.mod , 1:nv , -h.mod , col=c("BLUE","PURPLE")[1+(h.mod==h.med)] , lwd=2)
    segments(h.mod, 0:(nv-1) , h.mod , 1:nv, col=c("BLUE","PURPLE")[1+(h.mod==h.med)] , lwd=2)  
  ##  segments(v.mod , -(0:(nh-1)) , v.mod , -(1:nh) , col=c("BLUE","PURPLE")[1+(v.mod==v.med)] , lwd=2)
    segments((0:(nh-1)), v.mod , (1:nh) , v.mod , col=c("BLUE","PURPLE")[1+(v.mod==v.med)] , lwd=2)  
  }
  
  if(!is.null(dom)) {
    for(i in 1:nrow(dom)) {
      for(j in 1:ncol(dom)) {
        ## polygon(j - c(0.7,0.5,0.5,0.7,0.7) ,  c(0.6,0.6,0.4,0.4,0.6) - i , col=c("RED","GREEN")[1+as.numeric(dom[i,j])])
        polygon(i-c(0.7,0.5,0.5,0.7,0.7), j - c(0.6,0.6,0.4,0.4,0.6) ,   col=c("RED","GREEN")[1+as.numeric(dom[i,j])])  
      }
    }
  }
  if(!is.null(admis)) {
    for(i in 1:nrow(admis)) {
      for(j in 1:ncol(admis)) {
        ## polygon(j - c(0.3,0.5,0.5,0.3,0.3) ,  c(0.6,0.6,0.4,0.4,0.6) - i , col=c("RED","GREEN")[1+as.numeric(admis[i,j])])
        polygon(i-c(0.3,0.5,0.5,0.3,0.3), j - c(0.6,0.6,0.4,0.4,0.6) ,  col=c("RED","GREEN")[1+as.numeric(admis[i,j])])  
      }
    }
  }
  if(!is.null(n) & !is.null(r)) {
    for(i in 1:nrow(n)) {
      for(j in 1:ncol(n)) {
        ##text(j - 0.5 ,  0.7 - i , paste(r[i,j],"/",n[i,j]) )
        text(i - 0.5 ,  j - 0.7, paste(r[i,j],"/",n[i,j]) )  
        if(i==j) mtext(i,at=j-0.5,side=1,line=1)
        if(i==j) mtext(j,at=i-0.5,side=2,line=1)
      }
    }
  }
  if(!is.null(mat.lik)) {
    for(i in 1:nrow(mat.lik)) {
      for(j in 1:ncol(mat.lik)) {
        ##text(j-0.5 , 0.3-i , paste("p<theta","=",round(1-mat.lik[i,j],2),sep="" ) )
        text(i-0.5 , j-0.3 , paste("p<theta","=",round(1-mat.lik[i,j],2),sep="" ) )   
      }
    }
  }  
  if(!is.null(means)) {
    for(i in 1:nrow(means)) {
      for(j in 1:ncol(means)) {
        #text(j-0.5 , 0.1-i , paste("p(DLT)","=",round(1-means[i,j],2),sep="" ) )
        text(i-0.5 , j-0.1 , paste("p(DLT)","=",round(1-means[i,j],2),sep="" ) )  
      }
    }
  }  
  mtext("Drug A level", side=1, outer=TRUE,line=2)
  mtext("Drug B level", side=2, outer=TRUE,line=2)
  par(opar)
  list(h.med , v.med , h.mod , v.mod)
}

w.median <- function (x, w) {
  if (missing(w)) 
    w <- rep(1, length(x))
  ok <- complete.cases(x, w)
  x <- x[ok]
  w <- w[ok]
  ind <- sort.list(x)
  x <- x[ind]
  w <- w[ind]
  ind1 <- min(which(cumsum(w)/sum(w) >= 0.5))
  ind2 <- if ((w[1]/sum(w)) > 0.5) {
    1
  }
  else {
    max(which(cumsum(w)/sum(w) <= 0.5))
  }
  max(x[ind1], x[ind2])
}