
# turn sampling probabilities and sample counts into a synthetic population
probtodist<-function(classes,nums,prob,n,hajek=TRUE){
	#classes: size value of each class 
	#nums: number of members of that class in the sample
	#prob: Probabilities of selection for a member of that class
	props<-rep(0,length(classes))
	for(k in 1:length(classes)){
		props[k]<-nums[k]/prob[k]
	}
	# Hajek is the one in Krista's thesis and it implicitly multiplies the
	# HT estimate of the proportions by N as an estimate of the population size.
	# Otherwise use the HT estimator which implicitly multiples by the HT
	# estimate of the population size)
	#
	# The Hajek estimator is significantly less efficient than the HT estimator
	# under PPSWOR when the outcome is approximately proportional to degree.
	# When the disease status is unrelated to the degree the Hajek is slightly
	# more efficient (all assuming that the population size is actually N)
	#
	# In this case we are estimating the counts in each class not the disease
	# outcome
	#
	if(!hajek){
		n<-round(sum(props))
	}
	props<-props/sum(props)
	if(any(is.na(props))){
		print("Error in proptodist")
	}	
	#classes: size value of each class: 
	#props: proportion of the population in that class
	list(classes=classes, props=props, n=n)	
}




getestCstacked<-function(samp,n,nit=5, nsampsamp=1000, trace=FALSE, hajek=TRUE,SS.infinity=0.04){
	# samp: sample
	# n: population size
	classes<-sort(unique(samp))
	nums<-table(samp)
	nsamp<-length(samp)
	prob<-classes/sum(samp)
	prob <- prob * sum(nums/prob) / n
	
	#temp$classes: size value of each class 
	#temp$nums: number of members of that class 
	#temp$prob: Probabilities of selection for a member of that class
	temp<-probtodist(classes=classes,nums=nums,prob=prob,n=n,hajek=hajek)
	#temp$classes: size value of each class 
	#temp$props: proportion of the population in that class
	if(n*SS.infinity < nsamp){
		newprobs<-getinclCstacked(classes=temp$classes,props=temp$props,n=n,
				nsamp=nsamp,
				nsampsamp=nsampsamp,
				nums=nums)
		#newprobs$degvec: size value of each class
		#newprobs$pvec: Probabilities of selection for a member of that class
		
		for(i in 1:nit){
			if(!hajek){
				cat(sprintf("The estimated population size is %d.\n",temp$n))
			}
			temp<-probtodist(classes=newprobs$degvec,nums=nums,prob=newprobs$pvec,
					n=temp$n,hajek=hajek)
			newprobs<-getinclCstacked(classes=temp$classes,props=temp$props,
					n=temp$n,
					nsamp=nsamp,
					nsampsamp=nsampsamp,
					nums=nums)
		}
	}else{
		if(nsamp < n){
			#    VH
			pvec=nsamp*(temp$classes)/(temp$n*sum(temp$props*temp$classes))
		}else{
			#    Arithmetic Mean
			pvec=rep(nsamp,length=temp$classes)/temp$n
		}
		newprobs <- list(degvec=temp$classes, pvec=pvec,
				nbyclass=round(temp$n*temp$props/sum(temp$props)))
	}
	temp<-probtodist(classes=newprobs$degvec,nums=nums,prob=newprobs$pvec,
			n=temp$n,hajek=hajek)
	if(!hajek){
		cat(sprintf("The estimated population size is %d.\n",temp$n))
	}
	#classes: size value of each class 
	#props: proportion of the population in that class
	#probs: Probabilities of selection for a member of that class
	list(probs=newprobs$pvec,props=temp$props,classes=temp$classes,n=temp$n)
}


getinclCstacked<-function(classes,props,n,nsamp,nsampsamp,nums,hajek=TRUE){
	props2<-props/sum(props)
	#artificially inflate each value to at least 1 node
	nbyclass<-round(n*props2)
	nbyclass[nbyclass==0]<-1
	offby<-n-sum(nbyclass)
        if(is.na(offby)){
                print("Error in getincl: offby")
        }
	# based on estimate of degrees, estimate true inclusion probs by class
	Cret <- .C("getinclCstacked",
			nbyclass=as.integer(nbyclass),
			size=as.double(classes),
			props2=as.double(props2),
			offby=as.integer(offby),
			N=as.integer(n),
			K=as.integer(length(classes)),
			n=as.integer(nsamp),
			samplesize=as.integer(nsampsamp),
			Nk=as.integer(classes),
			PACKAGE="RDS")
	Ninf <- (Cret$Nk+nums)/(sum(Cret$Nk)+nsamp)
	Ninf <- nsamp*Ninf/Cret$nbyclass
	list(degvec=classes,pvec=Ninf,nbyclass=Cret$nbyclass)
}

