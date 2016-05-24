	
	
	################################################################################################################################################################################
	############################################################Hybrid search for solution using single protein#####################################################################
	############################################################The following program is to privide optimal design solutions for a multi-stage proteomic study######################
	############################################################Irene SL Zeng supervised by Thomas Lumley, March 2011###############################################################
	################################################################################################################################################################################
	
	#Function to generate the number of detected true positive at the final stage 
	#A protein file with beta(fold changes or any kind of efficacy measures from a paired sample design) and sigma (variance of the efficacy measures))  (ranges of both parameters are estimated from the actual pilot data, efficacy can be close, the variance is 2 times of the pilot 
	#make.some.data is a simulation of selecting proteins from its population and output the group f test p value and individual t test p value 
	#do.one.experiment is a simulation of three stage design in which groups of candidates are selected based on their group F test results as well as the individual t test result. 
	
	ots.env <- new.env()
	#protein<-read.csv("c:/selected_proteins_nooverlap.csv",sep=",",header=T)
	#library(MASS)

	#FUNTION 1: Ttest 
	Ttest<-function(sample,n,m)
	
		#Step 1: 
		# A random selection of sample means from multivariate normal distribution with beta and sigma as the mean and standard deviation of mean difference.
	
		{	
			
		#Retrieving variables from the wrap-up function

			budget<-get("BUDGET",envir=ots.env )
			sigma <- get("SIGMA", envir=ots.env )
			beta.artifact <- get("BETA.ARTIFACT", envir=ots.env )
			proteinID<-get("PROTEINID",envir=ots.env )
			sample.matrix<-get("SAMPLE.MATRIX",envir=ots.env )
			n1<-get("N1",envir=ots.env )
			s<-get("S",envir=ots.env )

			proteingr2=sample;proteinsample=t(sample)[,1:m]  					#sample is the protein sample with m proteins of row and n patients of column and a group variable 
			
		#Step 2: 
		#calculates the individual t statistics and provides the p values for g individual protein
		
			beta.sample=rowMeans(proteingr2[,1:n]); sigma.sample=apply(proteinsample,2,sd);tstat=beta.sample/sigma.sample*sqrt(n)
			proteingr3<-cbind(proteingr2,beta.sample,sigma.sample)				#attach mean difference beta and standard deviation sigma to the new simulated protein file 
			prob.t=pt(abs(tstat),df=n-1,lower.tail=FALSE)*2
			
		#Provides the p values for each protein
			c(prob.t)
		}

	
	#FUNTION 2: calculate.n3: Function to find n3 
	calculate.n3<-function(n2,p2,p3)
			{
	     		
			#Retrieving variables from the wrap-up function
			budget<-get("BUDGET",envir=ots.env )
			s<-get("S",envir=ots.env)
			cost2<-get("COST2.FUNCTION",envir=ots.env)
			cost3<-get("COST3.FUNCTION",envir=ots.env)
			recruit<-get("RECRUIT",envir=ots.env )

			return(((budget-s)-match.fun(cost2)(n2,p2)-recruit*n2)/(match.fun(cost3)(p3)+recruit))
			}



	#FUNTION 3: do.one.experiment
	#Screen each stage to select those with t test P value < stage I threshold 
	
	
	do.one.experiment.t<-function(initial)

		{
	
			# Retrieving variables from the wrap-up function
			budget<-get("BUDGET",envir=ots.env )
			sigma <- get("SIGMA", envir=ots.env )
			beta.artifact <- get("BETA.ARTIFACT", envir=ots.env )
			proteinID<-get("PROTEINID",envir=ots.env )
			sample.matrix<-get("SAMPLE.MATRIX",envir=ots.env )
			n1<-get("N1",envir=ots.env )
			recruit<-get("RECRUIT",envir=ots.env )

			n2=round(initial[3]);p=dim(sample.matrix)[1];

			#Sampling from data frame 
			sample.stage1<-sample(c(1:10000),n1)
			sample.data=sample.matrix[,sample.stage1]
			t.pvalue=Ttest(sample=sample.data,n=n1,m=p)
													
			in.stage2<-(t.pvalue<initial[1]) 						#select proteins and groups
			p2<-sum(in.stage2)
			beta2<-beta.artifact[in.stage2]
			sigma2<-sigma[in.stage2]								
		#number of proteins selected from stage I 
			proteinid2<-proteinID[in.stage2]
			
		#Stage 2 is a technical verification , so the selection is more favorite towards selecting individual protein
			#set.seed(20)
			sample.stageII=sample(c(1:10000),n2)
			sample.data2=sample.matrix[in.stage2,sample.stageII]
			
			t.pvalue2=Ttest(sample=sample.data2,n=n2,m=p2)
			in.stage3<-(t.pvalue2<initial[2]) 						#select proteins and groups; of note the selection criteria is different from stage I
			beta3<-beta2[in.stage3]
			sigma3<-sigma2[in.stage3]
			p3<-sum(in.stage3)								#number of proteins selected from stage II 		
			proteinid3<-proteinid2[in.stage3]
			
		#Stage 3 using the cost constraint with slake term S, to calculate n3, usig n3 to calculate the E(num of true positive) 
				
			n3=calculate.n3(n2=n2,p2=p2,p3=p3)
						
		#Calculate the expected number of detected positive  
			if (n3>100)
			{
			final.power<-power.t.test(n3,delta=beta3,sd=sigma3,sig.level=0.05,type="paired",alternative="two.sided")$power
			final.stage<-(final.power>0.85)
			detect.positive=sum(final.power)
			p4=sum(final.stage)
			proteinid4<-proteinid3[final.stage]
				
			return(c(detect.positive=p4,n3=n3,p2=p2,p3=p3,p4=p4,proteinid2=proteinid2,proteinid3=proteinid3,proteinid4=proteinid4))
			}
			else return(c(0,0,p2,p3,1,rep(0,p2),rep(0,p3),rep(0,1)))
		}

	#FUNCTION 4 power : Simulated power measured by the expected number of 1000 experiments with the same stage I, II cut off p values and sample size  
	power.t<-function(initial,optimize=TRUE)
		{ 	

			#retrive variable from the wrap-up function
			budget<-get("BUDGET",envir=ots.env )
			sigma <- get("SIGMA", envir=ots.env )
			beta.artifact <- get("BETA.ARTIFACT", envir=ots.env )
			proteinID<-get("PROTEINID",envir=ots.env )
			sample.matrix<-get("SAMPLE.MATRIX",envir=ots.env )
			n1<-get("N1",envir=ots.env )
			recruit<-get("RECRUIT",envir=ots.env )
			no.protein=length(beta.artifact);num.col=6+no.protein
			expect.positive=c(rep(0,1000))
			n3.vect=c(rep(0,1000))
			SELECT1<-matrix(nrow=1000,ncol=num.col,rep(0,1000*num.col))
			SELECT2<-matrix(nrow=1000,ncol=num.col,rep(0,1000*num.col))
			SELECT3<-matrix(nrow=1000,ncol=num.col,rep(0,1000*num.col))

			for (i in 1:1000) 
			{
			set.seed(i*100)
			esolution=do.one.experiment.t(initial)
			
			expect.positive[i]=esolution[1]
			n3.vect[i]=esolution[2]
			p2=esolution[3]
			p3=esolution[4]
			p4=esolution[5]
			SELECT1[i,1:p2]=esolution[6:(6+p2-1)]
			SELECT2[i,1:p3]=esolution[(6+p2):(6+p2+p3-1)]
			SELECT3[i,1:p4]=esolution[(6+p2+p3):(6+p2+p3+p4-1)]

			}
			
			mean.expect.pos=mean(expect.positive[expect.positive>=0])
			mean.n3<-mean(n3.vect[n3.vect>=0])

			if(optimize) return(-mean.expect.pos) else {
			print(expect.positive); print(n3.vect);print(summary(expect.positive));print(summary(n3.vect));print(table(SELECT1));print(table(SELECT2));print(table(SELECT3))

			par(mfcol=c(3,1))
			hist(expect.positive,main="expected positive");hist(n3.vect,main="n3");plot(expect.positive,n3.vect,xlab="expected positive",ylab="n3")
			return(mean.n3)}
			
		}
	
	genseq.single<-function(initial)
	{
		R=initial[4]
		# The following program is to set the local boundary of the geometry searching area
		llim.a1.t<-initial[1]-R*(initial[6]-initial[5])
		ulim.a1.t<-initial[1]+R*(initial[6]-initial[5])
		if (llim.a1.t<initial[5]) llim.a1.t=initial[5]
		if (ulim.a1.t>initial[6]) ulim.a1.t=initial[6]

		llim.a2.t<-initial[2]-R*(initial[9]-initial[8])
		ulim.a2.t<-initial[2]+R*(initial[9]-initial[8])
		if (llim.a2.t<initial[8]) llim.a2.t=initial[8]
		if (ulim.a2.t>initial[9]) ulim.a2.t=initial[9]
	
		llim.n2<-round(initial[3]-R*(initial[12]-initial[11]))
		ulim.n2<-round(initial[3]+R*(initial[12]-initial[11]))
		if (llim.n2<initial[11]) llim.n2=initial[11]
		if (ulim.n2>initial[12]) ulim.n2=initial[12]

		a1.t<-seq (llim.a1.t,ulim.a1.t,by=initial[7])
		a2.t<-seq (llim.a2.t,ulim.a2.t,by=initial[10])
		n2<-seq (llim.n2,ulim.n2,by=initial[13])
		
		alpha.seed<-expand.grid(a1.t,a2.t,n2)
		names(alpha.seed)<-c("a1.t","a2.t","n2")
		
		#insert a term here to control for 3 stage false positive of single protein
		select<-(alpha.seed[,1]*alpha.seed[,2]*0.05<0.01)
		alpha.seed=alpha.seed[select,]
		rown=dim(alpha.seed)[1]
		#print(rown)
		#Random local search with uniform probability: select one sample from the combinations,can consider a different distribution for selection prob.  
		if (rown ==0) return(initial)
		else {changepoints<-sample(c(1:rown),size=1,replace=FALSE)
			alpha.n<-alpha.seed[changepoints,]
			initial<-c(alpha.n[1,1],alpha.n[1,2],alpha.n[1,3],R,initial[5],initial[6],initial[7],initial[8],initial[9],initial[10],initial[11],initial[12],initial[13])
			#print(initial)
			return(initial)}
	}
	
	#select neighbourhood in the global cube and between a time interval(The opitimization surface is like each joints of a big cube travel between different times !!) 
	#local search within the defined neighbour hood cube in a defined time interval
	#Beta distribution (alpha=4,beta=6) Random start 
	#Wrap up function with global parameters budget, protein parameters

	
	optim.two.stage.single<-function(budget,protein,n1,artifact,iter.number,assaycost2.function,assaycost3.function,recruit=100,s=1000,a1.t.min=0.01,a1.t.max=0.25,a1.step=0.025,a2.t.min=0.01,a2.t.max=0.05,a2.step=0.025,n2.min=100,n2.max=1000,n2.step=100) 
	{
		#Insert programs to generate 1000 datasets for power function. And generate the dataset from first stage study with artifact 	
		#Fixed a random seed
		set.seed(100)
		sigma_m=diag(protein$sigma^2)#when we assume protein are independant, the covariance matrix is an diagnal matrix 
		BETA.ARTIFACT=protein$beta*artifact
		SIGMA=protein$sigma
		PROTEINID=protein$proteinid
		sample.frame=mvrnorm(n=10000,BETA.ARTIFACT,sigma_m)
		sample.t=t(sample.frame); 
		SAMPLE.MATRIX=cbind(sample.t,protein$group)
		N1=n1
		BUDGET<-budget;RECRUIT<-recruit;S<-s
		COST2.FUNCTION<-match.fun(assaycost2.function)
		COST3.FUNCTION<-match.fun(assaycost3.function)

	
		#Assigning variables for the wrap-up function environment
		assign("BUDGET",BUDGET,envir =ots.env )
		assign("SIGMA", SIGMA,envir = ots.env )
		assign("BETA.ARTIFACT",BETA.ARTIFACT,envir =ots.env )
		assign("PROTEINID",PROTEINID, envir =ots.env )
		assign("SAMPLE.MATRIX",SAMPLE.MATRIX,envir =ots.env )
		assign("N1",N1,envir =ots.env )
		assign("S",S,envir =ots.env )
		assign("RECRUIT",RECRUIT,envir =ots.env)
		assign("COST2.FUNCTION",COST2.FUNCTION,envir =ots.env)
		assign("COST3.FUNCTION",COST3.FUNCTION,envir =ots.env)

		#Generate the global search geometry areas, obtain the total number of combinations rown
		a1.t<-seq (a1.t.min,a1.t.max,by=a1.step)
		a2.t<-seq (a2.t.min,a2.t.max,by=a2.step)
			
		n2<-seq(n2.min,n2.max,by=n2.step)
		alpha.seed<-expand.grid(a1.t,a2.t,n2)
		names(alpha.seed)<-c("a1.t","a2.t","n2")
		rown=dim(alpha.seed)[1]
	
		solution.previous=0
		start.previous<-round(runif(1,min=0,max=1)*rown)
       	iter=1
		p.par=c(rep(0,13))
		solution.matrix=matrix(nrow=iter.number,ncol=14)

		while(iter <iter.number)
		{
			#Assign the next address based on the current address and obtain the radius of next local search geometry area. The distance between current and next address is beta distributed
			set.seed(100+iter)
			radius<-runif(1,min=0,max=1)
			if (radius > 0.5) start<-start.previous+round((rown-start.previous)*pbeta(runif(1),4,20))+1
     			if (radius <= 0.5) start<-start.previous-round(start.previous*pbeta(runif(1),4,20))+1
			if (start > rown) start<-rown 
	
			alpha.n<-alpha.seed[start,]
			#print(start)
			#print(alpha.n)
	
			initial<-c(alpha.n[1,1],alpha.n[1,2],alpha.n[1,3],radius,a1.t.min,a1.t.max,a1.step,a2.t.min,a2.t.max,a2.step,n2.min,n2.max,n2.step)
	
			solution<-optim(initial,power.t,genseq.single,method="SANN",control=list(maxit=1000,temp=10,trace=TRUE,REPORT=1000))

			if (-solution.previous < -solution$value){solution.previous=solution$value
										p.par=solution$par}
      		start.previous<-start
			solution.matrix[iter,1]=solution.previous
			solution.matrix[iter,2:14]=p.par
			iter=iter+1
			print(solution.previous)
			print(p.par)
	
		}
		return(solution.matrix)
		}

#assaycost2=function(n,p){280*p+1015*n}
#assaycost3=function(p){200*p}
#optim.two.stage.single(budget=1200000,artifact=rep(1,52),protein=protein,n1=30,iter.number=10,assaycost2.function=assaycost2,assaycost3.function=assaycost3,n2.min=30,n2.max=100,n2.step=10) 
	


