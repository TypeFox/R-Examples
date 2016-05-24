	
	################################################################################################################################################################################
	############################################################Hybrid search with group information of different proteins group pattern and budget###################################
	############################################################The following program is to privide optimal design solutions for a multi-stage proteomic study######################
	############################################################Irene SL Zeng supervised by Thomas Lumley, March 2011###############################################################
	################################################################################################################################################################################
	
	#protein<-read.csv("c:/selected_proteins_nooverlap.csv",sep=",",header=T)
	#library(MASS)
	#Function to generate the number of detected true positive at the final stage 
	#Function 1: Ftest.Ttest() is a function to calculate the  group f test p value and individual t test p value . It is a sub-function in do.one.experiment
	#Function 2: do.one.experiment() is function to simulate the three stage design in which groups of candidates are selected based on their group F test results as well as the individual t test result. It is a sub-
	#function in power(). It returns the number of true positives discovered through the three stage design. 
	#Function3:  calculate.n3() is a Function to find stage III sample size n3 based on the (maximum cost - slake term). It is a sub-function in do.one.experiment
	#Function4:  genseq() is a function to generate the sub-space of solutions of stage I,II, t and f test p values and stage II sample size n2. It is used for optim() and have a constraint for controlling  false positive. 
	#Function 5: power() is a function to simulate the 3 stage study process and calculate the expected  number of true positive 
	
	#FUNTION 1: Ftest.Ttest 
	#Function to generate the number of detected true positive at the final stage 
	#A protein file with beta(fold changes or any kind of efficacy measures from a paired sample design) and sigma (variance of the efficacy measures))(ranges of both parameters are estimated from the actual pilot data, efficacy can be close, the variance is 2 times of the pilot 
	
    	
	#FUNTION 1: Ftest.Ttest 
	Ftest.Ttest<-function(sample,n,m)
	
		#Step 1: 
		# A random selection of sample means from multivariate normal distribution with beta and sigma as the mean and standard deviation of mean difference.
	                                                       
		{
			#Retrieving variables from the wrap-up function

			budget<-get("BUDGET",envir=ots.env)
			sigma <- get("SIGMA", envir=ots.env)
			beta.artifact <- get("BETA.ARTIFACT", envir=ots.env)
			proteinID<-get("PROTEINID",envir=ots.env)
			sample.matrix<-get("SAMPLE.MATRIX",envir=ots.env)
			n1<-get("N1",envir=ots.env)
			s<-get("S",envir=ots.env)
			
			proteingr2=sample;proteinsample=t(sample[,1:n]);group=sample[,n+1] 		#proteinsample is the transposed sample with columns of proteins one row per patient 					#sample is the protein sample with m proteins of row and n patients of column and a group variable 
			
		#Step 2: 
		#calculates the individual t statistics and provides the p values for g individual protein
		
			beta.sample=rowMeans(proteingr2[,1:n]); sigma.sample=apply(proteinsample,2,sd);tstat=beta.sample/sigma.sample*sqrt(n)
			proteingr3<-cbind(proteingr2,beta.sample,sigma.sample)				#attach mean difference beta and standard deviation sigma to the new simulated protein file 
			prob.t=pt(abs(tstat),df=n-1,lower.tail=FALSE)*2
			
		#calculates the group F statistics and provides the p values for each group
			#matrix version of calculating F statistics, data need to sort by group 
			p=table(group);no.group=length(p);Fstat2<-c(rep(0,no.group));index<-c(rep(0,(no.group+1)))
			jj=1;group.name=as.numeric(names(p))
			
			while (jj <= no.group) 
			{
			group.sample=proteingr2[group==group.name[jj],1:n] 					#select the sample with n numbers of observations, do not exclude the group variable
			s0=var(t(group.sample)); 
			Fstat2[jj]=tryCatch(	{inv.s=solve(s0)
							tsqr=n*t(beta.sample[group==group.name[jj]])%*%inv.s%*%beta.sample[group==group.name[jj]]      #Hotelling t statistics	
							tsqr*(n-p[jj])/(p[jj]*(n-1))},error=function(e) 0)
			jj=jj+1
			}		#Provides the p values for each group using the derived F statisics above
					
			prob.f2=pf(Fstat2,df1=p,df2=n-p,lower.tail=FALSE)
			c(prob.f2,prob.t)
		}


	#FUNTION 3: calculate.n3: Function to find n3 
	calculate.n3<-function(n2,p2,p3)
			{
	      	#Retrieving variables from the wrap-up function
			
			budget<-get("BUDGET",envir=ots.env)
			s<-get("S",envir=ots.env)
			cost2<-get("COST2.FUNCTION",envir=ots.env)
			cost3<-get("COST3.FUNCTION",envir=ots.env)
			recruit<-get("RECRUIT",envir=ots.env )

			return(((budget-s)-match.fun(cost2)(n2,p2)-recruit*n2)/(match.fun(cost3)(p3)+recruit))
			}


	#FUNTION 2: do.one.experiment

	#Screen each stage to select those with Group F p value < 0.05 and individual T p value< alpha.t(initial[1]) Or Group F p value < alpha.f (initial[2])
	
	do.one.experiment<-function(initial,optimize=TRUE)
		{

			# Retrieving variables from the wrap-up function
			
			budget<-get("BUDGET",envir=ots.env)
			sigma <- get("SIGMA", envir=ots.env)
			beta.artifact <- get("BETA.ARTIFACT", envir=ots.env)
			proteinID<-get("PROTEINID",envir=ots.env)
			sample.matrix<-get("SAMPLE.MATRIX",envir=ots.env)
			n1<-get("N1",envir=ots.env)
			
			group=sample.matrix[,10001] 								   	#sample matrix is the simulated popultion protein data generated with columns of proteins and a group variable 
			N2=round(initial[5]);p=dim(sample.matrix)[1];no.group1=length(table(group))	#p is the number of protein, no.group1 records num of proteins in each group 
					
			#Sampling from data frame 
			sample.stage1<-sample(c(1:10000),n1)
			sample.data=cbind(sample.matrix[,sample.stage1],group)
			pvalue=Ftest.Ttest(sample=sample.data,n=n1,m=p)
			f.pvalue=pvalue[1:no.group1]
			t.pvalue=pvalue[(no.group1+1):(no.group1+p)]
										
			index1<- match( f.pvalue[f.pvalue<0.05],f.pvalue);index2<-match(f.pvalue[f.pvalue<initial[2]],f.pvalue)						#match indexs of groups that satisfied f statiscs < 0.05 or f statistics < alpha.f1
			instage2a<-(group %in% index1);instage2b<-t.pvalue %in% t.pvalue[t.pvalue<initial[1]]; instage2c<-(group %in% index2) 	#generate indexs of groups that satisfied f statiscs < 0.05 or f statistics < alpha.f1

			in.stage2<-((instage2b & instage2a)|instage2c) 						#select proteins and groups
			no.group2=length(table(group[in.stage2]))							#number of groups	
			P2<-sum(in.stage2)
			beta2<-beta.artifact[in.stage2]
			sigma2<-sigma[in.stage2]	
			proteinid2<-proteinID[in.stage2]
			if (P2 ==0) {print("no protein selected");return(c(0,0,1,1,1,rep(0,1),rep(0,1),rep(0,1)))}
			else {
			#number of proteins selected from stage I 
	
			#Stage 2 is a technical verification , so the selection is more favorite towards selecting individual protein
			#set.seed(20)
			sample.stageII=sample(c(1:10000),N2)
			group2=group[in.stage2]
			sample.data2=cbind(sample.matrix[in.stage2,sample.stageII],group2)
			
			pvalue2=Ftest.Ttest(sample=sample.data2,n=N2,m=P2)
			fstat2<-pvalue2[1:no.group2]; tstat2<-pvalue2[(no.group2+1):(no.group2+P2)]
			
			index3<-match(fstat2[fstat2<initial[4]],fstat2)
			instage3a<-(group2 %in% index3);instage3b<-tstat2 %in% tstat2[tstat2<initial[3]]; instage3c<-tstat2 %in% tstat2[tstat2<0.05]
			
			in.stage3<-((instage3a & instage3c)|instage3b) 						#select proteins and groups; of note the selection criteria is different from stage I
			
			no.group3=length(table(group2[in.stage3]))
			P3<-sum(in.stage3)										#number of proteins selected from stage II 		
			if (P3 ==0) {print("no protein selected");return(c(0,0,1,1,1,rep(0,1),rep(0,1),rep(0,1)))}	
			else{	beta3<-beta2[in.stage3]
			sigma3<-sigma2[in.stage3]
			proteinid3<-proteinid2[in.stage3]

		#Stage 3 using the cost constraint with slake term S, to calculate n3, usig n3 to calculate the E(num of true positive) 
		#set.seed(30)
		#stage3<-F test. t test(beta=beta2[in.stage3],sigma=sigma2[in.stage3],group=group2[in.stage3],n=n3,m=p3)
		#fstat3<-stage3[1:no.group3]; tstat3<-stage3[(no.group3+1):(no.group3+p3)]
		
			group3=group2[in.stage3]
			N3=calculate.n3(n2=N2,p2=P2,p3=P3)
					
	#Calculate the expected number of detected positive  
			if (N3>100)
			{
			final.power<-power.t.test(N3,delta=beta3,sd=sigma3,sig.level=0.05,type="paired",alternative="two.sided")$power
			final.stage<-(final.power>0.85)
			detect.positive=sum(final.power)
			P4=sum(final.stage)
			proteinid4<-proteinid3[final.stage]
					
			return(c(detect.positive,n3=N3,p2=P2,p3=P3,p4=P4,proteinid2=proteinid2,proteinid3=proteinid3,proteinid4=proteinid4))
			}
			else return(c(0,0,1,1,1,rep(0,1),rep(0,1),rep(0,1)))
		}
	}	
	}

	#FUNCTION 5 power : Simulated power measured by the expected number of 1000 experiments with the same stage I, II cut off p values and sample size
	#initial[1]: stage I t test p value cut-off; initial[2]: stage I f test p value cut-off; initial[3]: stage II t test p value cut-off; initial[4]: stage II f test p value cut-off
	power<-function(initial,optimize=TRUE)
		{ 	
			#retrive variable from the wrap-up function
			budget<-get("BUDGET",envir=ots.env )
			sigma <- get("SIGMA", envir=ots.env )
			beta.artifact <- get("BETA.ARTIFACT", envir=ots.env )
			proteinID<-get("PROTEINID",envir=ots.env )
			sample.matrix<-get("SAMPLE.MATRIX",envir=ots.env )
			n1<-get("N1",envir=ots.env )
			cost2<-get("COST2.FUNCTION",envir=ots.env)
			cost3<-get("COST3.FUNCTION",envir=ots.env)
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
			esolution=do.one.experiment(initial)
			#print(esolution)
			expect.positive[i]=esolution[1]
			n3.vect[i]=esolution[2]
			p2=esolution[3]
			p3=esolution[4]
			p4=esolution[5]
			#print(p2);print(p3);print(p4)
			if(p2 !=0 & p3 !=0 & p4!=0)
				{SELECT1[i,1:p2]=esolution[6:(6+p2-1)]
				SELECT2[i,1:p3]=esolution[(6+p2):(6+p2+p3-1)]
				SELECT3[i,1:p4]=esolution[(6+p2+p3):(6+p2+p3+p4-1)]}
			else print ("NO protein selected at final stage")

			}
			
			mean.expect.pos=mean(expect.positive[expect.positive>=0])
			mean.n3<-mean(n3.vect[n3.vect>=0])

			if(optimize) return(-mean.expect.pos) 
			else {
			print("no. of expected positive:");print(expect.positive)
			print("n3:");print(n3.vect)
			print("summary of expected positive:");print(summary(expect.positive))
			print("n3 vector"); print(summary(n3.vect));
			print(table(SELECT1));print(table(SELECT2));print(table(SELECT3))
			par(mfcol=c(3,1))
			hist(expect.positive,main="expected positive");hist(n3.vect,main="n3");plot(expect.positive,n3.vect,xlab="expected positive",ylab="n3")
			n2=round(initial[5]);
			stage2.cost[i]=match.fun(cost2)(n2,p2)+recruit*n2;stage3.cost[i]=match.fun(cost3)(p3)*n3.vect[i]+recruit*n3.vect[i]
			print(paste("cost at stage II:",stage2.cost),quote=FALSE)
			print(paste("cost at stage III:",stage3.cost),quote=FALSE)
			return(mean.n3)
				}
		}

	
	#Function4: genseq() is a function to generate the sub-space of solutions of stage I,II, t and f test p values and stage II sample size n2.
	genseq.group<-function(initial)
	{
		R=initial[6]
		# The following program is to set the local boundary of the geometry searching area
		llim.a1.t<-initial[1]-R*(initial[8]-initial[7])
		ulim.a1.t<-initial[1]+R*(initial[8]-initial[7])
		if (llim.a1.t<initial[7]) llim.a1.t=initial[7]
		if (ulim.a1.t>initial[8]) ulim.a1.t=initial[8]

		llim.a1.f<-initial[2]-R*(initial[10]-initial[9])
		ulim.a1.f<-initial[2]+R*(initial[10]-initial[9])
		if (llim.a1.f<initial[9]) llim.a1.f=initial[9]
		if (ulim.a1.f>initial[10]) ulim.a1.f=initial[10]
	
		llim.a2.t<-initial[3]-R*(initial[13]-initial[12])
		ulim.a2.t<-initial[3]+R*(initial[13]-initial[12])
		if (llim.a2.t<initial[12]) llim.a2.t=initial[12]
		if (ulim.a2.t>initial[13]) ulim.a2.t=initial[13]
	
		llim.a2.f<-initial[4]-R*(initial[15]-initial[14])
		ulim.a2.f<-initial[4]+R*(initial[15]-initial[14])
		if (llim.a2.f<initial[14]) llim.a2.f=initial[14]
		if (ulim.a2.f>initial[15]) ulim.a2.f=initial[15]

		llim.n2<-round(initial[5]-R*(initial[18]-initial[17]))
		ulim.n2<-round(initial[5]+R*(initial[18]-initial[17]))
		if (llim.n2<initial[17]) llim.n2=initial[17]
		if (ulim.n2>initial[18]) ulim.n2=initial[18]

		a1.t<-seq (llim.a1.t,ulim.a1.t,by=initial[11])
		a1.f<-seq (llim.a1.f,ulim.a1.f,by=initial[11])
	
		a2.t<-seq (llim.a2.t,ulim.a2.t,by=initial[16])
		a2.f<-seq (llim.a2.f,ulim.a2.f,by=initial[16])

		n2<-seq (llim.n2,ulim.n2,by=initial[19])
		alpha.seed<-expand.grid(a1.t,a1.f,a2.t,a2.f,n2)
		names(alpha.seed)<-c("a1.t","a1.f","a2.t","a2.f","n2")
		
		#insert a term here to control for 3 stage false positive of single protein
		select<-(alpha.seed[,1]*alpha.seed[,3]*0.05<0.01)
		alpha.seed=alpha.seed[select,]
		rown=dim(alpha.seed)[1]

		#print(rown)
		#Random local search with uniform probability: select one sample from the combinations,can consider a different distribution for selection prob.  
		if (rown ==0) return(initial)
		else 
		{
		changepoints<-sample(c(1:rown),size=1,replace=FALSE)
		alpha.n<-alpha.seed[changepoints,]
		initial<-c(alpha.n[1,1],alpha.n[1,2],alpha.n[1,3],alpha.n[1,4],alpha.n[1,5],R,initial[7],initial[8],initial[9],initial[10],initial[11],initial[12],initial[13],initial[14],initial[15],initial[16],initial[17],
		initial[18],initial[19])
		#print(initial)
		return(initial)}
	}
	
	#select neighbourhood in the global cube and between a time interval(The opitimization surface is like each joints of a big cube travel between different times !!) 
	#local search within the defined neighbourhood cube in a defined time interval
	#Beta distribution (alpha=4,beta=6) Random start 
	#Wrap up function with global parameters BUDGET, protein parameters
	ots.env <- new.env()
	optim.two.stage.group<-function(budget,protein,n1,artifact,iter.number,assaycost2.function,assaycost3.function,recruit=100,s=1000,a1.t.min=0.01,a1.t.max=0.25,a1.f.min=0.01,a1.f.max=0.25,a1.step=0.025,a2.t.min=0.01,a2.t.max=0.05,a2.f.min=0.05,a2.f.max=0.05,a2.step=0.025,n2.min=100,n2.max=1000,n2.step=100)
	{
		#Insert programs to generate 1000 datasets for power function, and generate the dataset from first stage study with artifact 	
		#Fixed a random seed
		set.seed(100)
		sigma_m=diag(protein$sigma^2)				#when we assume protein are independant, the covariance matrix is an diagnal matrix 
		BETA.ARTIFACT<-protein$beta*artifact
		SIGMA<-protein$sigma
		PROTEINID<-protein$proteinid
		sample.frame=mvrnorm(n=10000,BETA.ARTIFACT,sigma_m)
		sample.t=t(sample.frame)
		SAMPLE.MATRIX<-cbind(sample.t,protein$group)
		N1=n1
		BUDGET<-budget;RECRUIT<-recruit;S<-s
		COST2.FUNCTION<-match.fun(assaycost2.function)
		COST3.FUNCTION<-match.fun(assaycost3.function)

		#Assigning variables for the wrap-up function environment
		assign("BUDGET",BUDGET,envir =ots.env)
		assign("SIGMA", SIGMA,envir = ots.env)
		assign("BETA.ARTIFACT",BETA.ARTIFACT,envir =ots.env)
		assign("PROTEINID",PROTEINID, envir =ots.env)
		assign("SAMPLE.MATRIX",SAMPLE.MATRIX,envir =ots.env)
		assign("N1",N1,envir =ots.env)
		assign("S",S,envir =ots.env)
		assign("RECRUIT",RECRUIT,envir =ots.env)
		assign("COST2.FUNCTION",COST2.FUNCTION,envir =ots.env)
		assign("COST3.FUNCTION",COST3.FUNCTION,envir =ots.env)

		#Generate the global search geometry areas, obtain the total number of combinations rown
		a1.t<-seq (a1.t.min,a1.t.max,by=a1.step)
		a1.f<-seq (a1.f.min,a1.f.max,by=a1.step)
		a2.t<-seq (a2.t.min,a2.t.max,by=a2.step)
		a2.f<-seq (a2.f.min,a2.f.max,by=a2.step)
	
		n2<-seq(n2.min,n2.max,by=n2.step)
		alpha.seed<-expand.grid(a1.t,a1.f,a2.t,a2.f,n2)
		names(alpha.seed)<-c("a1.t","a1.f","a2.t","a2.f","n2")
		rown=dim(alpha.seed)[1]
	
		solution.previous=0
		start.previous<-round(runif(1,min=0,max=1)*rown)
       	iter=1
		p.par=c(rep(0,19))
		solution.matrix=matrix(nrow=iter.number,ncol=20)

		while(iter <iter.number)
		{
			#Assign the next address based on the current address and obtain the radius of next local search area. The distance between current and next address is beta distributed
			set.seed(100+iter)
			radius<-runif(1,min=0,max=1)
			if (radius > 0.5) start<-start.previous+round((rown-start.previous)*pbeta(runif(1),4,20))+1
     			if (radius <= 0.5) start<-start.previous-round(start.previous*pbeta(runif(1),4,20))+1
			if (start > rown) start<-rown 
	
			alpha.n<-alpha.seed[start,]
			print(start)
			print(alpha.n)
			initial<-c(alpha.n[1,1],alpha.n[1,2],alpha.n[1,3],alpha.n[1,4],alpha.n[1,5],radius,a1.t.min,a1.t.max,a1.f.min,a1.f.max,a1.step,a2.t.min,a2.t.max,a2.f.min,a2.f.max,a2.step,n2.min,n2.max,n2.step)
	
			#print(power(initial))
			solution<-optim(initial,power,genseq.group,method="SANN",control=list(maxit=1000,temp=10,trace=TRUE,REPORT=1000))
			if (-solution.previous < -solution$value){solution.previous=solution$value
								p.par=solution$par}
      		start.previous<-start
			solution.matrix[iter,1]=-solution.previous
			solution.matrix[iter,2:20]=p.par
			iter=iter+1
			print(solution.previous)
			print(p.par)
			
		}
		return(solution.matrix)
		}
#assaycost2=function(n,p){280*p+1015*n}
#assaycost3=function(p){200*p}
#optim.two.stage.group(budget=6000000,protein=protein,n1=60,artifact=rep(1,52),iter.number=2,assaycost2.function=assaycost2,assaycost3.function=assaycost3)

