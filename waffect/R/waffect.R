waffect <- function(prob, count, label, method=c("backward","mcmc","reject"), burnin){
	
	if(missing(count)){
		stop('count is missing')
	}
	
	if(missing(prob)){
		if(length(count)==1){
			stop('prob is missing and the total number of individuals is unknown')
		}
		if(length(count)==2){
			warning('prob is missing: for all individuals the probability to be a case will be set to 0.1')
			prob=rep(0.1,sum(count))
		}
		if(length(count)>=3){
			warning('prob is missing: for all individuals the probability to be in class i=1,...,K will be set to 0.1 for i=1 and to 0.9/(K-1) for i=2,...,K where K = total number of classes')
			prob=matrix(rep(c(0.1,rep(0.9/(length(count)-1),length(count)-1)),sum(count)),ncol=sum(count))
		}
	}
		
	if(missing(label)){
		if(length(count)==1){
			warning('label is missing: cases/controls will be coded by 1/0')
			label=c(1,0) 
		}
		if(length(count)==2){
			warning('label is missing: by default label=c(1,0)') 
			label=c(1,0)
		}
		if(length(count)>=3){
			warning('label is missing: by default label=1:K, where K is the lenght of count (i.e. the total number of classes)')
			label=1:length(count)
		}
	}	
		
	if(is.vector(prob)){
		if(length(count)>=3){
			stop('prob is a vector: in this case count must be an integer (the number of cases) or a length 2 vector (the number of cases and the one of controls)')
		}
		if(length(count)==1){
			if(length(prob)<count){
				stop('count is an integer: in this case the length of prob must be >= count')
			}
		}
		if(length(count)==2){
			if(length(prob)!=sum(count)){
				stop('count is a length 2 vector: in this case the length of prob must be equal to the sum of the entries of count (i.e. the total number of individuals)')
			}
		}
		if(sum(prob>1 | prob<0)>0){
			stop('Entries in prob must be probabilities')
		}
		if(length(label)>=3 | length(label)==1){
			stop('prob is a vector: in this case label must be a length 2 vector (codes for cases and controls)') 
		}
	}
	
	if(is.matrix(prob)){
		if(nrow(prob)!=length(count)){
			stop('prob is a matrix: in this case the length of count must be equal to the number of rows of prob (i.e. the number of classes)')
		}
		if(ncol(prob)!=sum(count)){
			stop('prob is matrix: in this case the number of columns of prob must be equal to the sum of the entries of count (i.e. the total number of individuals)')
		}
		if(sum(prob>1 | prob<0)>0){ 
			stop('Entries in prob must be probabilities')
		}
		if(sum(apply(prob,2,sum)!=1)>0){
			stop('Entries in prob must be probabilities (entries in coulumns must add up to one)')
		}
		if(length(label)!=length(count)){
			stop('prob is a matrix: in this case label and count must have the same length') 
		}
	}	
	
	n=sum(count)
	K=length(label) 
		
	#method:
	if(missing(method)){
		#warning('Backward sampling is the method by default')
		method='backward'
	}
	if(method == 'mcmc'){
		if(missing(burnin)){
			warning('burnin for mcmc is missing: by default burnin = 10^5 * total number of individuals')
			burnin = 100000*n
		}
	}
	if(method == 'reject'){
		warning('Rejection algorithm is deprecated: expect very slow running time and possibly no answer at all')
	}
	
	#call R functions:
	if(K==2){
		res <- waffectbin(prob = prob, count = count, label = label, method = method, burnin = burnin)
	}
	if(K>2){
		res=rep(NA,n)
		# main loop
		for (k in 1:(K-1)) {
			# affect k versus the remaining unaffected
			# prepare data for waffectbin call
			p=prob[k,]/apply(prob[k:K,],2,sum)
			
			res[is.na(res)]=waffectbin(prob=p[is.na(res)],count=count[k],label=c(label[k],NA),method=method,burnin=burnin)
		}
	res[is.na(res)]=label[K]
	}       
    
	return(res)    
}    
    
#Message d'erreur/warnings?    
#Verificare risultati non cambiano
#debugging?
#Documentation   
#Compilare
#Inviare a Gregory
    
    
    