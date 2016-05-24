
# A power function that can be used to obtain the average cost of stage II and III from the simulation function

		power.single.cost<-function(initial,protein,n1,artifact,budget,s,assaycost2.function,assaycost3.function,recruit,optimize=FALSE)
		{ 	
			set.seed(100)
			#library(MASS)
			sigma_m=diag(protein$sigma^2)				
			BETA.ARTIFACT<-protein$beta*artifact
			SIGMA<-protein$sigma
			PROTEINID<-protein$proteinid
			sample.frame=mvrnorm(n=10000,BETA.ARTIFACT,sigma_m)
			sample.t=t(sample.frame)
			SAMPLE.MATRIX<-cbind(sample.t,protein$group)
			N1=n1
			BUDGET<-budget;RECRUIT<-recruit;S=s;
			COST2.FUNCTION<-match.fun(assaycost2.function)
			COST3.FUNCTION<-match.fun(assaycost3.function)

	#FUNTION 1: Ttest 
		Ttest<-function(sample,n,m,beta.artifact=BETA.ARTIFACT,sigma=SIGMA,sample.matrix=SAMPLE.MATRIX,n1=N1,recruit=RECRUIT,cost2=COST2.FUNCTION,cost3=COST3.FUNCTION,proteinid=PROTEINID)

	
		#Step 1: 
		# A random selection of sample means from multivariate normal distribution with beta and sigma as the mean and standard deviation of mean difference.
	
		{	
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

		calculate.n3<-function(n2,p2,p3,budget=BUDGET,s=S,cost2=COST2.FUNCTION,cost3=COST3.FUNCTION,recruit=RECRUIT)
			{
	      	
			return(((budget-s)-match.fun(cost2)(n2,p2)-recruit*n2)/(match.fun(cost3)(p3)+recruit))
			}

#FUNTION 3: do.one.experiment
	#Screen each stage to select those with t test P value < stage I threshold 
		
	do.one.experiment.t<-function(initial,budget=BUDGET,s=S,sigma=SIGMA,beta.artifact=BETA.ARTIFACT,proteinID=PROTEINID,sample.matrix=SAMPLE.MATRIX,n1=N1,cost2=COST2.FUNCTION,cost3=COST3.FUNCTION,proteinid=PROTEINID,recruit=RECRUIT)

		{
	
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
				
			n3=calculate.n3(n2=n2,p2=p2,p3=p3,budget=BUDGET,s=S,cost2=COST2.FUNCTION,cost3=COST3.FUNCTION,recruit=RECRUIT)
						
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



			expect.positive=c(rep(0,1000))
			n3.vect=c(rep(0,1000))
			n2=round(initial[3])
			no.protein=length(BETA.ARTIFACT);num.col=6+no.protein

			stage2.cost=c(rep(0,1000));stage3.cost=c(rep(0,1000))
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
			if(p2 !=0 & p3 !=0 & p4!=0)
				{SELECT1[i,1:p2]=esolution[6:(6+p2-1)]
				SELECT2[i,1:p3]=esolution[(6+p2):(6+p2+p3-1)]
				SELECT3[i,1:p4]=esolution[(6+p2+p3):(6+p2+p3+p4-1)]}
			else print ("NO protein selected at final stage")
			stage2.cost[i]=match.fun(COST2.FUNCTION)(n2,p2)+recruit*n2;stage3.cost[i]=match.fun(COST3.FUNCTION)(p3)*n3.vect[i]+recruit*n3.vect[i]

			}
			
			mean.expect.pos=mean(expect.positive[expect.positive>=0])
			mean.n3<-mean(n3.vect[n3.vect>=0])
			mean.stage2.cost=mean(stage2.cost);mean.stage3.cost=mean(stage3.cost)
			if(optimize) return(-mean.expect.pos) 
			else {
			print("no. of expected positive:");print(expect.positive)
			print("n3:");print(n3.vect)
			print("summary of expected positive:");print(summary(expect.positive))
			print("n3 vector"); print(summary(n3.vect));
			print(table(SELECT1));print(table(SELECT2));print(table(SELECT3))
			par(mfcol=c(3,1))
			hist(expect.positive,main="expected positive",xlim=c(0,(mean.expect.pos+1)));hist(n3.vect,main="n3",xlim=c(0,(mean.n3+100)));plot(expect.positive,n3.vect,xlab="expected positive",ylab="n3")
			print(paste("cost at stage II:",mean.stage2.cost),quote=FALSE)
			print(paste("cost at stage III:",mean.stage3.cost),quote=FALSE)
			return(c(mean.n3,mean.stage2.cost,mean.stage3.cost))
				}
		}

#assaycost2=function(n,p){280*p+1015*n}
#assaycost3=function(p){200*p}
#initial=c(0.01,0.01,100)
#power.single.cost(initial,protein=protein,artifact=rep(1,52),n1=30,budget=1200000,s=1000,assaycost2.function=assaycost2,assaycost3.function=assaycost3,recruit=100,optimize=FALSE)
	

