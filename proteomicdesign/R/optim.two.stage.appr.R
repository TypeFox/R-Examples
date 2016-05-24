	
	
	################################################################################################################################################################################
	############################################################Hybrid search with group information with different proteins group pattern and cost#################################
	############################################################Use approximation for the analytical function#######################################################################
	############################################################Irene SL Zeng supervised by Thomas Lumley, Dec, 2011################################################################
	################################################################################################################################################################################
	
	#Function 2: calculate.cost() is a function to estimate the cost
	#function 3: power.appr(). A function returns the number of true positives discovered through the three-stage design using approximation of the analytical funtion 
	#Function 4: genseq() is a function to generate the sub-space of solutions of stage I,II, T and F test p values and stage II sample size n2. It is used for optim() and have a constraint for controlling false positive. 
	#A protein file with beta(fold changes or any kind of efficacy measures from a paired sample design) and sigma(variance of the efficacy measures) 
	#protein<-read.csv("c:/selected_proteins_top5.csv",sep=",",header=T)
	#library(MASS)
	#FUNTION 2: calculate.cost: Function to estimate the cost
	calculate.cost<-function(n2,n3,p1,p2)
			{
			cost2<-get("COST2.FUNCTION",envir=ots.env)
			cost3<-get("COST3.FUNCTION",envir=ots.env)
			recruit<-get("RECRUIT",envir=ots.env )

			return((n2+n3)*recruit+match.fun(cost2)(n2,p1)+match.fun(cost3)(p2)*n3)
			}

	#FUNCTION 3 power: analytical power function   
	power.appr<-function(initial)
		{ 	
			
			#Retrieve variable from wrap-up function
			#utils::globalVariables(BUGET,SIGMA,BETA.ARTIFACT,PROTEINID,SAMPLE.MATRIX,N1,RECRUIT,PROTEINGR2,GROUP,NO.GROUP,HOTEL.T2,P)
			budget<-get("BUDGET",envir=ots.env)
			sigma<-get("SIGMA",envir=ots.env)
			beta.artifact<-get("BETA.ARTIFACT",envir=ots.env)
			proteinID<-get("PROTEINID",envir=ots.env)
			sample.matrix<-get("SAMPLE.MATRIX",envir=ots.env)
			n1<-get("N1",envir=ots.env)
			recruit<-get("RECRUIT",envir=ots.env)
			proteingr2<-get("PROTEINGR2",envir=ots.env)		
			group<-get("GROUP",envir=ots.env)
			no.group<-get("NO.GROUP",envir=ots.env)
			hotel.t2<-get("HOTEL.T2",envir=ots.env)
			p<-get("P",envir=ots.env)		
			#parameters for optimization
			n2=initial[5];n3=initial[6];a1.t=initial[1];a1.f=initial[2];a2.t=initial[3];a2.f=initial[4]
			
			#Stage I: Step 2 
			#Calculates the individual t statistics(ncp) and type II error of stage I
			#Attach mean difference beta and standard deviation sigma to the new simulated protein file 
				 							     		
			p1=length(beta.artifact) 
			tstat=beta.artifact/sigma*sqrt(n1)
			d.t1=qt(a1.t,df=n1-1,lower.tail=FALSE,ncp=0)
			prob.t=pt(d.t1,df=n1-1,lower.tail=TRUE,ncp=abs(tstat))      						#Type II error 
			
			#Calculates the group F statistics and type II error for each group
			#Matrix version of calculating F statistics, need to sort data by group 
			
			Fstat2<-c(rep(0,no.group))
			jj=1
			while (jj <= no.group) 
			{
							Fstat2[jj]=tryCatch(	{
									tsqr=n1*hotel.t2[jj]  							#Hotelling t statistics	
										 },error=function(e) 0)
							jj=jj+1
			}																#Provides the p value for each group using the derived F statisics above, the F statistics is also the ncp
			
			Fstat22<-c(rep(0,p1))
			jj=1
			
			while (jj <= p1) #each protein will have a F statistics T tau associated with it
			{
							group.sample=proteingr2[group==group[jj],1:n1] 					#Select the sample with n numbers of observation, do not include the group variable
							s0=var(t(group.sample)); 
							beta.artifact2=abs(beta.artifact)
							beta.artifact2[jj]=0
							Fstat22[jj]<-tryCatch(	{inv.s=solve(s0)
							tsqr=n1*t(beta.artifact2[group==group[jj]])%*%inv.s%*%beta.artifact2[group==group[jj]]  #Hotelling t statistics	
							},error=function(e) 0)
							jj=jj+1
			}		
							
									
							d.f1=qf(a1.f,df1=p,df2=n1-p,lower.tail=FALSE,ncp=0) 
							#print(d.f1)  
							d.f.05=qf(0.05,df1=p,df2=n1-p,lower.tail=FALSE,ncp=0)
							#print(d.f.05)
							prob.f1.prob.t=rep(0,p1) 												#F conditional on T 
							low.bound=min((d.t1-abs(tstat)))-10
			for (i in 1:p1)
			{
							t1<-seq(low.bound,(d.t1-abs(tstat[i])),by=0.1)
							#print(t1)
							d1.f.appr=d.f1[group[i]]-(n1-p[group[i]])/(p[group[i]]*(n1-1))*(t1)^2
							#print(d1.f.appr)
							ft=pf(d1.f.appr,df1=p[group[i]],df2=n1-p[group[i]],lower.tail=TRUE,ncp=Fstat22[i])*dt(t1,df=n1-1,ncp=abs(tstat[i]))
							#print(ft)
							spline.f<-splinefun(t1,ft) 
							prob.f1.prob.t[i]=integrate(spline.f,min(t1),max(t1))$value
							#Approximation of the intergral of marginal distribution of the group statistics under condition of individual stats threshold
			}	

							prob.f1=pf(d.f1,df1=p,df2=n1-p,lower.tail=TRUE,ncp=Fstat2)   						#Type II error of F test at d.f1 significant level
							prob.f.05=pf(d.f.05,df1=p,df2=n1-p,lower.tail=TRUE,ncp=Fstat2)   					#Type II error of F test at 0.05 significant level
							prob.f.com=pmin(prob.f1,prob.f.05)
						
			#Step 3; 
			#Overall stage I power (1-type II error)
			type.II.error1=c(rep(0,p1));
			for (i in 1:p1)
 			{type.II.error1[i]=prob.t[i]*prob.f1.prob.t[i]+prob.f.com[group[i]]-prob.t[i]*prob.f.com[group[i]]}
			
			power.1=1-type.II.error1
			power.1[power.1<0]=0;power.1[power.1>1]=1
	
			
			#Step 4: Stage II 
			tstat2=beta.artifact/sigma*sqrt(n2)
			d.t2=qt(a2.t,df=n2-1,lower.tail=FALSE,ncp=0)
			d.t.05=qt(0.05,df=n2-1,lower.tail=FALSE,ncp=0)

			prob.t2=pt(d.t2,df=n2-1,lower.tail=TRUE,ncp=abs(tstat2))     
			prob.t.05=pt(d.t.05,df=n2-1,lower.tail=TRUE,ncp=abs(tstat2)) 
			prob.t.com=pmin(prob.t2,prob.t.05)
				
			#group statistics
			Fstat3<-c(rep(0,no.group)); jj=1
			while (jj <= no.group) 
			{
							Fstat3[jj]=tryCatch(	{
							tsqr=n2*hotel.t2[jj]  									#Hotelling t statistics	
							},error=function(e) 0)
			jj=jj+1
			}																#Provides the f statistics and the F statistics is also the ncp
			
			#group statistics given individual statisitics threshold
			Fstat32<-c(rep(0,p1))
			jj=1
			
			while (jj <= p1) 
			{
			group.sample=proteingr2[group==group[jj],1:n1] 									#Select the sample with n numbers of observations, do not include the group variable
			s0=var(t(group.sample)); 
			beta.artifact2=beta.artifact
			beta.artifact2[jj]=0
			Fstat32[jj]<-tryCatch(	{inv.s=solve(s0)
							tsqr=n2*t(beta.artifact2[group==group[jj]])%*%inv.s%*%beta.artifact2[group==group[jj]]  #Hotelling t statistics	
							},error=function(e) 0)
			jj=jj+1
			}		

			#F conditional on T 
			#Conditional group probability based on individual test statistics > individual threshold
			d.f2=qf(a2.f,df1=p,df2=n2-p,lower.tail=FALSE,ncp=0)
			prob.f2=pf(d.f2,df1=p,df2=n2-p,lower.tail=TRUE,ncp=Fstat3)      	
			prob.f2.prob.t=rep(0,p1) 							
			low.bound.2=min((d.t2-abs(tstat)))-10

			for (i in 1:p1)
			{
			t1<-seq(low.bound.2,(d.t2-abs(tstat[i])),by=0.1)
			d.f2.appr=d.f2[group[i]]-(n2-p[group[i]])/(p[group[i]]*(n2-1))*(t1)^2
			ft2=pf(d.f2.appr,df1=p[group[i]],df2=n2-p[group[i]],lower.tail=TRUE,ncp=Fstat32[i])*dt(t1,df=n2-1,ncp=abs(tstat[i]))
			spline.f2<-splinefun(t1,ft2) 
			prob.f2.prob.t[i]=integrate(spline.f2,min(t1),max(t1))$value
			}	
						
			#Overall stage II power (1-type II error)
			type.II.error2=c(rep(0,p1));
			for (i in 1:p1)
 			{type.II.error2[i]=prob.t2[i]*prob.f2.prob.t[i]+prob.t.com[i]-prob.t.com[i]*prob.f2[group[i]]}
			#print("type.II.error stage II ");#print(type.II.error2)

			power.2=1-type.II.error2
			power.2[power.2<0]=0;power.2[power.2>1]=1

			#Step 5
			
			tstat3=beta.artifact/sigma*sqrt(n3)
			d.t3.05=qt(0.05,df=n3-1,lower.tail=FALSE,ncp=0)
			prob.t3=pt(d.t3.05,df=n3-1,lower.tail=TRUE,ncp=abs(tstat3))  
			power.3=1-prob.t3
			power.3[power.3<0]=0
			
			#Step 6
			power.overall=sum(power.1*power.2*power.3)
			cost=calculate.cost(n2,n3,sum(power.1),sum(power.1*power.2))
			#print(cost);#print(power.overall)
			if (cost<budget) return(-power.overall)
			else return(0)
		}

	
	genseq.appr<-function(initial)
	{
		R=initial[7]
		# The following program is to set the local boundary of the geometry searching area
		llim.a1.t<-initial[1]-R*(initial[9]-initial[8])
		ulim.a1.t<-initial[1]+R*(initial[9]-initial[8])
		if (llim.a1.t<initial[8]) llim.a1.t=initial[8]
		if (ulim.a1.t>initial[9]) ulim.a1.t=initial[9]

		llim.a1.f<-initial[2]-R*(initial[11]-initial[10])
		ulim.a1.f<-initial[2]+R*(initial[11]-initial[10])
		if (llim.a1.f<initial[10]) llim.a1.f=initial[10]
		if (ulim.a1.f>initial[11]) ulim.a1.f=initial[11]
	
		llim.a2.t<-initial[3]-R*(initial[14]-initial[13])
		ulim.a2.t<-initial[3]+R*(initial[14]-initial[13])
		if (llim.a2.t<initial[13]) llim.a2.t=initial[13]
		if (ulim.a2.t>initial[14]) ulim.a2.t=initial[14]
	
		llim.a2.f<-initial[4]-R*(initial[16]-initial[15])
		ulim.a2.f<-initial[4]+R*(initial[16]-initial[15])
		if (llim.a2.f<initial[15]) llim.a2.f=initial[15]
		if (ulim.a2.f>initial[16]) ulim.a2.f=initial[16]

		llim.n2<-round(initial[5]-R*(initial[19]-initial[18]))
		ulim.n2<-round(initial[5]+R*(initial[19]-initial[18]))
		if (llim.n2<initial[18]) llim.n2=initial[18]
		if (ulim.n2>initial[19]) ulim.n2=initial[19]


		llim.n3<-round(initial[6]-R*(initial[22]-initial[21]))
		ulim.n3<-round(initial[6]+R*(initial[22]-initial[21]))
		if (llim.n3<initial[21]) llim.n3=initial[21]
		if (ulim.n3>initial[22]) ulim.n3=initial[22]

		a1.t<-seq (llim.a1.t,ulim.a1.t,by=initial[12])
		a1.f<-seq (llim.a1.f,ulim.a1.f,by=initial[12])
	
		a2.t<-seq (llim.a2.t,ulim.a2.t,by=initial[17])
		a2.f<-seq (llim.a2.f,ulim.a2.f,by=initial[17])

		n2<-seq (llim.n2,ulim.n2,by=initial[20])
		n3<-seq (llim.n3,ulim.n3,by=initial[23])
		alpha.seed<-expand.grid(a1.t,a1.f,a2.t,a2.f,n2,n3)
		names(alpha.seed)<-c("a1.t","a1.f","a2.t","a2.f","n2","n3")
		
		#insert a term here to control for 3 stage false positive of single protein
		select<-(alpha.seed[,1]*alpha.seed[,3]*0.05<0.01)
		alpha.seed=alpha.seed[select,]
		rown=dim(alpha.seed)[1]

		
		#Random local search with uniform probability: select one sample from the combinations,can consider a different distribution for selection prob.  
		if (rown ==0) return(initial)
		else 
		{
		changepoints<-sample(c(1:rown),size=1,replace=FALSE)
		alpha.n<-alpha.seed[changepoints,]
		initial<-c(alpha.n[1,1],alpha.n[1,2],alpha.n[1,3],alpha.n[1,4],alpha.n[1,5],alpha.n[1,6],R,initial[8],initial[9],initial[10],initial[11],initial[12],initial[13],initial[14],initial[15],initial[16],initial[17],initial[18],
		initial[19],initial[20],initial[21],initial[22],initial[23])
		
		return(initial)}
	}
	
	#select neighbourhood in the global cube and between a time interval(The opitimization surface is like each joints of a big cube travel between different times !!) 
	#local search within the defined neighbour hood cube in a defined time interval
	#Beta distribution (alpha=4,beta=6) Random start 
	#Wrap up function with global parameters BUDGET, protein parameters
	ots.env <- new.env()
	optim.two.stage.appr<-function(budget,protein,n1,artifact,iter.number,assaycost2.function,assaycost3.function,recruit=100,a1.t.min=0.01,a1.t.max=0.25,a1.f.min=0.01,a1.f.max=0.25,a1.step=0.025,a2.t.min=0.01,a2.t.max=0.05,a2.f.min=0.05,a2.f.max=0.05,a2.step=0.025,n2.min=100,n2.max=1000,n2.step=100,n3.min=1000,n3.max=5000,n3.step=100) 
	{
		#Insert programs to generate 1000 datasets for power function, and generate the dataset from first stage study with artifact 	
		#Fixed a random seed
		set.seed(100)
		sigma_m=diag(protein$sigma^2)		    										#When we assume protein are independant, the covariance matrix is an diagnal matrix 
		BETA.ARTIFACT=protein$beta*artifact   										#These are all set as global parameters
		SIGMA=protein$sigma
		PROTEINID=protein$proteinid
		GROUP=protein$group
		N1=n1
		BUDGET<-budget;RECRUIT<-recruit
		COST2.FUNCTION<-match.fun(assaycost2.function)
		COST3.FUNCTION<-match.fun(assaycost3.function)


		sample.frame=mvrnorm(n=10000,BETA.ARTIFACT,sigma_m)
		sample.t=t(sample.frame);SAMPLE.MATRIX<-cbind(sample.t,protein$group)
		sample.stage1<-sample(c(1:10000),N1)
		sample.data=cbind(SAMPLE.MATRIX[,sample.stage1],protein$group)
		PROTEINGR2<-sample.data
		
		P<-table(protein$group); NO.GROUP<-length(P); HOTEL.T2<-c(rep(0,NO.GROUP)); index<-c(rep(0,(NO.GROUP+1)))
		jj=1;group.name=as.numeric(names(P))
			
			while (jj <= NO.GROUP) 
			{
			group.sample=PROTEINGR2[GROUP==group.name[jj],1:N1] 							#Select the sample with n numbers of observation, do not include the group variable
			s0=var(t(group.sample)); 
			HOTEL.T2[jj]<-tryCatch(	{inv.s=solve(s0)
							tsqr=t(BETA.ARTIFACT[GROUP==group.name[jj]])%*%inv.s%*%BETA.ARTIFACT[GROUP==group.name[jj]]  
							},error=function(e) 0)								#Hotelling t statistics	
			jj=jj+1
			}															#Provides the partical ncp
			
		#assign variables in the wrap-up function environment
		assign("BUDGET",BUDGET,envir=ots.env)
		assign("SIGMA", SIGMA,envir=ots.env)
		assign("BETA.ARTIFACT",BETA.ARTIFACT,envir=ots.env)
		assign("PROTEINID",PROTEINID,envir=ots.env)
		assign("SAMPLE.MATRIX",SAMPLE.MATRIX,envir=ots.env)
		assign("N1",N1,envir=ots.env)
		assign("RECRUIT",RECRUIT,envir =ots.env)
		assign("COST2.FUNCTION",COST2.FUNCTION,envir =ots.env)
		assign("COST3.FUNCTION",COST3.FUNCTION,envir =ots.env)
		assign("PROTEINGR2",PROTEINGR2,envir=ots.env)		
		assign("GROUP",GROUP,envir=ots.env)
		assign("NO.GROUP",NO.GROUP,envir=ots.env)
		assign("HOTEL.T2",HOTEL.T2,envir=ots.env)
		assign("P",P,envir=ots.env)

		#Generate the global search geometry areas, obtain the total number of combinations rown
		a1.t<-seq (a1.t.min,a1.t.max,by=a1.step)
		a1.f<-seq (a1.f.min,a1.f.max,by=a1.step)
		a2.t<-seq (a2.t.min,a2.t.max,by=a2.step)
		a2.f<-seq (a2.f.min,a2.f.max,by=a2.step)
	
		n2<-seq(n2.min,n2.max,by=n2.step)
		n3<-seq(n3.min,n3.max,by=n3.step)
		alpha.seed<-expand.grid(a1.t,a1.f,a2.t,a2.f,n2,n3)
		names(alpha.seed)<-c("a1.t","a1.f","a2.t","a2.f","n2","n3")
		rown=dim(alpha.seed)[1]
		
		#Initialization
		solution.previous=0
		start.previous<-round(runif(1,min=0,max=1)*rown)
       	iter=1
		p.par=c(rep(0,23))
		solution.matrix=matrix(nrow=iter.number,ncol=24)

		while(iter < iter.number)
		{
			#Assign the next address based on the current address and obtain the radius of next local search geometry area. The distance between current and next address is beta distributed
			set.seed(100+iter*50)
			radius<-runif(1,min=0,max=1)
			if (radius > 0.5) start<-start.previous+round((rown-start.previous)*pbeta(runif(1),4,20))+1
     			if (radius <= 0.5) start<-start.previous-round(start.previous*pbeta(runif(1),4,20))+1
			if (start > rown) start<-rown 
	
			alpha.n<-alpha.seed[start,]
			#print(start)
			#print(alpha.n)
	
			initial<-c(alpha.n[1,1],alpha.n[1,2],alpha.n[1,3],alpha.n[1,4],alpha.n[1,5],alpha.n[1,6],radius,a1.t.min,a1.t.max,a1.f.min,a1.f.max,a1.step,a2.t.min,a2.t.max,a2.f.min,a2.f.max,a2.step,n2.min,n2.max,n2.step,n3.min,n3.max,n3.step)
			solution<-optim(initial,power.appr,genseq.appr,method="SANN",control=list(maxit=1000,temp=20,trace=TRUE,REPORT=1000))

			if (-solution.previous < -solution$value){solution.previous=solution$value
										p.par=solution$par}
      		start.previous<-start
			solution.matrix[iter,1]=solution.previous
			solution.matrix[iter,2:24]=p.par
			iter=iter+1
			print(solution.previous)
			print(p.par)
	
		}
	return(solution.matrix)
		}

#assaycost2=function(n,p){280*p+1015*n}
#assaycost3=function(p){200*p}
#optim.two.stage.appr(budget=500000,protein=protein,n1=30,artifact=rep(1,5),iter.number=3,assaycost2.function=assaycost2,assaycost3.function=assaycost3,recruit=100,n2.min=10,n2.max=100,n2.step=10,n3.min=100,n3.max=1000,n3.step=100) 



