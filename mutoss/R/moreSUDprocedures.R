# 
# Author: JonathanRosenblatt
###############################################################################



#---------- Service Functions---------#

reject<- function(sorted, criticals){
	m<- length(sorted)
	stopifnot( length(criticals) == m )
	indicators<-  sorted<criticals # Marking p-values below the critical values
	
	if(!any(indicators)) {
		return(list(cutoff=0,cut.index=0))
	}
	
	cut.index<- max((1:m)[indicators])
	
	cutoff<- sorted[cut.index] #The highest rejected p-value
	
	return( list(cutoff=cutoff,cut.index=cut.index) )
}

bh.adjust<- function(sorted, m, m0, constant=1){
	adjusted<- rep(NA,m)
	temp.min<- sorted[m]
	min.ind<- rep(0,m)    
	for (i in m:1) {
		temp<- min(m0*sorted[i]*constant / i, 1)
		if  ( temp <= temp.min  ) {
			temp.min <- temp
			min.ind[i]<- 1
		}
		adjusted[i]<- temp.min
	}
	return(adjusted)
}




linearStepUp<- function(sorted, q, m, adjust=FALSE, m0=m, pi0, constant=1){
	if(missing(m0) & !missing(pi0)) {
		m0=pi0*m
	}
	else{
		
		criticals<- (1:m)*q/(m0*constant)
		cutoff<- reject(sorted,criticals)
		rejected<-  sorted<=cutoff$cutoff
		
		adjusted=rep(NA,m)
		
		if(adjust) {
			adjusted<-bh.adjust(sorted,m=m,m0=m0,constant=constant)
		}
		
		multiple.pvals<- data.frame(
				original.pvals=sorted,
				criticals=criticals,
				rejected=rejected,
				adjusted.pvals=adjusted)
		output<- list(Cutoff=cutoff,Pvals=multiple.pvals)
		return(output)
	}
	
}

#---------------Two Stage------------------#
solve.q<- function(sorted, m, j, r){
	##TODO: [JR] correct problem: when all pvalues are rejected a might be negative
	a<- sorted*(m-r)/(1:m)
	#stopifnot(a>0 , j>=1 , j<=m , r>=0 , r<=m) 
	adjusted<- ifelse(a>0.5, 1 , a/(1-a) )
	temp.min<- adjusted[m]
	for(i in m:j){
		if(adjusted[i]<=temp.min) temp.min<- adjusted[i]
		else adjusted[i]<- temp.min
	}
	return(adjusted)
}# Close solve.q function

two.stage.adjust<- function(sorted, r=0, patience=4, m){
	adjusted<- rep(0,m)  
	
	# Adjusting sorted p-values
	adjusted.q<- solve.q(sorted=sorted,m=m,j=1,r=0)
	checking<- adjusted.q
	
	#Has the procedure rejected everything at the first stage?
	if(sum(linearStepUp(sorted,adjusted.q[1]/(1+adjusted.q[1]),m=m)$Pvals[['rejected']])==m){
		adjusted.q<- rep(adjusted.q[1],m)
		return(adjusted.q)
	}
	
	else{
		for (j in 1:m) { 
			delta.r<- 1
			delta.q<- 1
			new.q<- adjusted.q[j]
			r.new<- sum(linearStepUp(sorted,new.q/(1+new.q),m=m)$Pvals[['rejected']])
			counter<- 0
			max.q<- 0
			
			while(delta.r>0 & delta.q>0){
				old.q<- new.q
				r.old<- r.new
				new.q<- solve.q(sorted=sorted,m=m,j=j,r=r.old)[j]
				r.new<- sum(linearStepUp(sorted,new.q/(1+new.q),m=m)$Pvals[['rejected']])
				delta.r<- abs(r.new-r.old)  
				delta.q<- abs(new.q-old.q)
				counter<- counter+1
				if(counter>patience & max.q!=new.q) max.q<- max(max.q,new.q)
				else if(counter>patience & max.q==new.q ) break
			} #Close interations inside q[j]
			
			adjusted.q[j]<- min(new.q,1)
			adjusted.q[min(j+1,m)]<- adjusted.q[j]
			stopifnot(any(adjusted.q[j]<=checking[j]))
		}#Close looping over j.    
		
		temp.min<- adjusted.q[m]
		for(i in m:1){
			if(adjusted.q[i]<=temp.min) temp.min<- adjusted.q[i]
			else adjusted.q[i]<- temp.min
		}
		return(adjusted.q)      
	}#Close 'else' clause
	
}# Close two.stage.adjust


two.stage<- function(pValues, alpha){
	ranks<- rank(pValues)
	sorted<-sort(pValues)
	m<- length(sorted)
	
	#Stage I- estimating m0
	q1<- alpha/(1+alpha) 
	stage.one<- linearStepUp(sorted, q1, adjust=TRUE, m=m)
	r<- sum(stage.one$Pvals[['rejected']]) #count rejection
	if (r==0) {  #if nothing is rejected, return the results of the linear step up
		stage.one$Pvals[['adjusted.pvals']]<- 1
		return(stage.one)  
	}
	else if (r==m) {
		stage.one$Pvals[['adjusted.pvals']]<- stage.one$Pvals[['adjusted.pvals']][1]
		return(stage.one)  
	}
	
	#Stage II- updating alpha using m0
	else {
		m0<- m-r
		output<- linearStepUp(sorted=sorted,q=q1,m0=m0,m=m)
		output$Pvals[['adjusted.pvals']]<- two.stage.adjust(sorted, alpha, m=m)
		output<-output$Pvals[ranks,]
		output.2<- list(
				criticalValues=output$criticals,
				rejected=output$rejected,
				adjPValues=output$adjusted.pvals,
				errorControl=new(Class='ErrorControl',type="FDR",alpha=alpha),
				pi0= m0/m
		)
		return(output.2)
	}
}

#pvals<- runif(100,0,0.1)
#two.stage(pvals,0.1)

mutoss.two.stage<- function() { return(new(Class="MutossMethod",
					label="B.K.Y. (2006) Two-Stage Step-Up",
					errorControl="FDR",
					callFunction="two.stage",
					output=c("adjPValues", "criticalValues", "rejected", "pi0", "errorControl"),
					assumptions=c("Independent test statistics"),
					info="<h2>Benjamini-Krieger-Yekutieli (2006) Two-Stage Step-Up Procedure</h2>\n\n							
							<p>A p-value procedure which controls the FDR at level <i>&alpha;</i> for independent test statistics, in which case it is more powerful then non adaptive procedures such as the Linear Step-Up (BH). On the other hand, when this is not the case, no error control is guaranteed.
							The linear step-up procedure is used in he first stage to estimate the number of true null hypotheses (mo) which is plugged in a linear step-up
							procedure at the second stage. 
							<h3>Reference:</h3>
							<ul>
							<li>Benjamini, Y., Krieger, A. and Yekutieli, D. \"<i> Adaptive linear step-up procedures that control the false
							discovery rate. </i>\" Biometrika, 93(3):491-507, 2006. </li>\n
							</ul>",
					parameters=list(pValues=list(type="numeric"), alpha=list(type="numeric"))
			)) }

#---------------------Multistage Step-Down-------------------#

multiple.down.adjust<- function(sorted, m){
	adjusted<- rep(NA,m)
	temp.max<- sorted[1]
	max.ind<- rep(0,m)
	
	for (i in 1:m) {
		temp<-  min(sorted[i]*(m+1-i)/(i*(1-sorted[i])),1)
		if  ( temp >= temp.max  ) {
			temp.max <- temp
			max.ind[i] <- 1
		}
		adjusted[i]<- temp.max
	}
	
	return(adjusted)
}





multiple.down=function(pValues, alpha){
	sorted<- sort(pValues)
	ranks<- rank(pValues)
	m<- length(pValues)
	
	if(alpha>0.5) warning('FDR is not controlled when q>0.5')
	
	criticals<- sapply(1:m,function(i) alpha*i/(m-i*(1-alpha)+1))
	indicators<-  sorted<criticals # Marking p-values below the critical values
	
	if(!indicators[1]) cutoff<-list(cutoff=0,cut.index=0)
	else if(all(indicators)) cutoff<- list(cutoff=sorted[m],cut.index=m)
	else{ 
		cut.index<- min((1:m)[!indicators])-1
		cutoff<- list(cutoff=sorted[cut.index],cut.index=cut.index) 
	}
	
	rejected<-  sorted<=cutoff$cutoff
	
	adjusted<-multiple.down.adjust(sorted,m)  
	
	output<- data.frame(
			criticals=criticals,
			rejected=rejected,
			adjusted.pvals=adjusted)
	output<- output[ranks,]
	output.2<-list(
			criticalValues=output$criticals,
			rejected=output$rejected,
			adjPValues=output$adjusted.pvals,
			errorControl=new(Class='ErrorControl',type="FDR",alpha=alpha)
	)
	return(output.2)	
}

mutoss.multiple.down <- function() { return(new(Class="MutossMethod",
					label="B.K.Y. (2006) Multi-Stage Step-down",
					errorControl="FDR",
					callFunction="multiple.down",
					output=c("adjPValues", "criticalValues", "rejected", "errorControl"),
					assumptions=c("Independent test statistics"),
					info="<h2>Benjamini-Krieger-Yekutieli (2006) multi-stage step-down procedure</h2>\n\n\
							<p>A non-linear step-down p-value procedure which control the FDR for independent test statistics and enjoys more power then other non-adaptive procedure such as the linear step-up (BH).
							For the case of non-independent test statistics, non-adaptive procedures such as the linear step-up (BH) or the all-purpose conservative Benjamini-Yekutieli (2001) are recommended.</p>\n
							<h3>Reference:</h3>
							<ul>
							<li>Benjamini, Y., Krieger, A. and Yekutieli, D. \"<i> Adaptive linear step-up procedures that control the false
							discovery rate. </i>\" Biometrika, 93(3):491-507, 2006. </li>\n
							</ul>",
					parameters=list(pValues=list(type="numeric"), alpha=list(type="numeric"))
			)) }

