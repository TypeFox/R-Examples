# A power function that can be used to obtain the average cost of stage II and III from the simulation function


power.group.cost<-function(initial,protein,n1,artifact,budget,s,assaycost2.function,assaycost3.function,recruit,optimize=FALSE)
{ 	
			set.seed(100)
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

	
	#FUNTION 1: Ftest.Ttest 
	Ftest.Ttest<-function(sample,n,m,beta.artifact=BETA.ARTIFACT,sigma=SIGMA,sample.matrix=SAMPLE.MATRIX,n1=N1,recruit=RECRUIT,cost2=COST2.FUNCTION,cost3=COST3.FUNCTION,proteinid=PROTEINID)
	
		#Step 1: 
		# A random selection of sample means from multivariate normal distribution with beta and sigma as the mean and standard deviation of mean difference
	                                                       
		{
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

	#FUNTION 2: do.one.experiment

	#Screen each stage to select those with Group F p value < 0.05 and individual T p value< alpha.t(initial[1]) Or Group F p value < alpha.f (initial[2])
	
	do.one.experiment<-function(initial,optimize=TRUE,budget=BUDGET,s=S,beta.artifact=BETA.ARTIFACT,sigma=SIGMA,proteinid=PROTEINID,sample.matrix=SAMPLE.MATRIX,n1=N1,recruit=RECRUIT,cost2=COST2.FUNCTION,cost3=COST3.FUNCTION)

		{					
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
			proteinid2<-proteinid[in.stage2]
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
			N3=calculate.n3(n2=N2,p2=P2,p3=P3,budget=BUDGET,s=S,cost2=COST2.FUNCTION,cost3=COST3.FUNCTION,recruit=RECRUIT)
					
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
	
#FUNTION 3: calculate.n3: Function to find n3 

		calculate.n3<-function(n2,p2,p3,budget=BUDGET,s=S,cost2=COST2.FUNCTION,cost3=COST3.FUNCTION,recruit=RECRUIT)
			{
	      	
			return(((budget-s)-match.fun(cost2)(n2,p2)-recruit*n2)/(match.fun(cost3)(p3)+recruit))
			}


			expect.positive=c(rep(0,1000))
			n3.vect=c(rep(0,1000))
			n2=round(initial[5])
			no.protein=length(BETA.ARTIFACT);num.col=6+no.protein
			stage2.cost=c(rep(0,1000));stage3.cost=c(rep(0,1000))
			SELECT1<-matrix(nrow=1000,ncol=num.col,rep(0,1000*num.col))
			SELECT2<-matrix(nrow=1000,ncol=num.col,rep(0,1000*num.col))
			SELECT3<-matrix(nrow=1000,ncol=num.col,rep(0,1000*num.col))
			
			for (i in 1:1000) 
			{
			set.seed(i*100)
			esolution=do.one.experiment(initial)
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
#initial=c(0.11,0.18,0.01,0.05,92)
#initial=c(0.09,0.13,0.01,0.05,70)
#initial=c(0.01,0.18,0.01,0.05,100)
#power.group.cost(initial,protein=protein,artifact=rep(1,5),n1=30,budget=1200000,s=1000,assaycost2.function=assaycost2,assaycost3.function=assaycost3,recruit=100,optimize=FALSE)
	


