################################################################ 
### Function for calculating sample sizes for a 2-stage sample
### Author: Susie Jentoft & Johan Heldal
### LAst Edited: 14.11.2012

################################################################ 


#### Variables
#File 1:

#COUNTRY - input dataset   
#SURVEY - input dataset   
#STRATUM_ID - input dataset  
#STRATUM_NAME - input dataset 
#STRATUM_SIZE  (done)   
#ST1_ANT_SSU    (done) was z$n /dem$n /big$n
#ST1_NO_PSU     (done) z$M.PSU
#ST1_SEL_PSU    (done) z$m
#ST1_CV  (done) #Was ST1_SEL_PSU calculated using cost/Var? - 1=yes 2=no 
#ST1_CPSU    - input and calc (done)
#ST1_CSSU    - input and calc (done)
#ST1_VWITHIN  (done)  z$Vwithin
#ST1_VAMONG    (done) z$Vamong
#ST1_COST   (done) z$cost 
#DOMAINS (done)

####
#File 2: 

#STRATUM_ID - input
#PSU_SN - input
#PSU_SIZE - input
#ST1_PROB (done)
#ST1_SEL_PSU 
#DOMAIN_ID (done)
#DOMAIN_SIZE_PSU   (in PSU) (done) domain.pop
#ST2_PROB (done)
#ST2_ANT_SSU (done)
#ST2_SEL_SSU - needs to be a part of stage2Sample.fn (done)

###Add in
#DOMAIN_SIZE_STR (done) stratademo.pop 
#ST1_POPT (done) p.opt
#
##########################################

SampleSizes<-function(df1, df2, PSU1, strata1, size1, demo1, demoyesno, domainalloc, stratayesno, strata2,
		mkvalue, mkvar, ssyesno, ssvar, sstotal, ppstype, nminimum, 
		Parcost, PSUcost, PSUV, PSUmean)
	{
	x<-df1
	PSUname<-names(x[PSU1])

#if df2 missing - create dataset
	if(stratayesno==1) {
		strata.levels<-levels(as.factor(x[,strata1]))
		strata.levels<-strata.levels[order(strata.levels)]
		strata.id<-seq(1:length(strata.levels))
		z<-data.frame(strata.levels, strata.id)
		strata2<-1 
		if (ssyesno=="1" & mkvalue==1) {
			z$ST1_CPSU<-tapply(x[,PSUcost], list(x[,strata1]), mean)[]
			z$ST1_CSSU<-tapply(x[,Parcost], list(x[,strata1]), mean)[]
			PSUcost<-match("ST1_CPSU", names(z))
			Parcost<-match("ST1_CSSU", names(z))
			}
		} else {
		z<-df2
		z<-z[order(z[,strata2]),]
		z$strata.id<-row.names(z)
		}
	
#Standardise
	x<-x[order(x[,strata1]),]
	x$strata.id<-as.numeric(as.factor(x[,strata1]))		### new bit ### put in as.factor as as.numeric not working by itself

#add in calculated size
	z$STRATUM_SIZE<-as.numeric(xtabs(x[,size1]~x$strata.id))
	size2<-match("STRATUM_SIZE", names(z))
			
#add in proportional weight
	z$ww<-z[,size2]/sum(z[,size2])

#add in weighted mean of variable - at the moment only using this with optimal sampling and no demographics
if (mkvalue==1 & demoyesno==1){
	mean.Var<-NULL
	for (i in 1:max(x$strata.id)){
		temp<-subset(x, strata.id==i)
		a<-weighted.mean(temp[,PSUmean], w=temp[,size1])
		mean.Var<-c(mean.Var, a)
	}
	z$mean.Var<-mean.Var
	assign("T.Var", weighted.mean(x[,PSUmean], w=x[,size1]), inherits=TRUE )
	}

#Add in Vwithin and Vamong
if (mkvalue==1){
	temp<-data.frame(xtabs(x[,size1]~x[,strata1]))
	names(temp)<-c("strata","size")
	Nk<-temp$size[x$strata.id]
	mu.calc<-x[,size1]*x[,PSUmean]/Nk
	mu<-xtabs(mu.calc~x$strata.id)

	##Calculating Vwithin
	Vwith.calc<-x[,size1]*x[,PSUV]/Nk
	z$ST1_VWITHIN<-xtabs(Vwith.calc~x$strata.id)
	Vwithin<-match("ST1_VWITHIN",names(z))

	##Calculating Vamong
	Vamong.calc<-x[,size1]*((x[,PSUmean]-mu[x$strata.id])^2)/Nk
	z$ST1_VAMONG<-xtabs(Vamong.calc~x$strata.id)
	Vamong<-match("ST1_VAMONG",names(z))
	}

#n strata calculations
if (ssyesno=="0" & stratayesno==0){z$ST1_ANT_SSU<-z[,ssvar]}
if (ssyesno=="0" & stratayesno==1){
		temp<-cbind(x$strata.id, x[,ssvar])
		z$ST1_ANT_SSU<-unique(temp)[,2]}

if (ssyesno==1){
	if (ppstype==1|ppstype==2) z$ST1_ANT_SSU<-z$ww*sstotal
	if (ppstype==2) z$ST1_ANT_SSU<-ifelse (z$ST1_ANT_SSU < nminimum, nminimum, z$ST1_ANT_SSU)
	if (ppstype==3) {nfixed<-sstotal/nrow(z)
		z$ST1_ANT_SSU<-rep(nfixed, nrow(z))}
	if (ppstype==4) z$ST1_ANT_SSU<-sstotal*((sqrt(z$ST1_VWITHIN+z$ST1_VAMONG)*Nk)/sum(sqrt(z$ST1_VWITHIN+z$ST1_VAMONG)*Nk))
	}

#Add in number of PSU in each strata
	z$ST1_NO_PSU<-c(table(x[,strata1])[])

	if(mkvalue==1){
		z$Vratio<-z[,Vwithin]/z[,Vamong]
		z$Cratio<-z[,PSUcost]/z[,Parcost]
		z$ST1_POPT<-sqrt(z$Vratio*z$Cratio)
		z$ST1_SEL_PSU<-mk.calc(z$ST1_ANT_SSU,z$ST1_POPT,z$ST1_NO_PSU)
		}
	if (mkvalue==0 & stratayesno==0){z$ST1_SEL_PSU<-z[,mkvar]}
	if (mkvalue==0 & stratayesno==1){
		temp<-cbind(x$strata.id, x[,mkvar])
		z$ST1_SEL_PSU<-unique(temp)[,2]}

#calculate inclus prob 
	x$ST1_SEL_PSU<-z$ST1_SEL_PSU[x$strata.id]
	ST1_PROB<-matrix(NA, max(x$strata.id),ncol=1)
	for (i in 1:max(x$strata.id)){ST1_PROB[i]<-list(inclusionprobabilities(x[,size1][which(x$strata.id==i)], z$ST1_SEL_PSU[i]))}
	x$ST1_PROB<-c(ST1_PROB,recursive=TRUE)
	x$strata.pop<-z[,size2][x$strata.id] #Nk
	x$strata.n<-z$ST1_ANT_SSU[x$strata.id]

#rearrange dataset if demographics used
	if(demoyesno=="0"){
		dem<-x
		names(dem)[demo1]<-paste("domain",1:length(demo1),sep="")
		dem<-reshape(dem,direction="long",varying=demo1,sep="",v.names="DOMAIN_SIZE_PSU",times=names(x[,demo1]))
		dem$DOMAIN_ID<-as.factor(dem$time)
		dem$strata.pop<-z[,size2][dem$strata.id] #Nk
		demographSums<-xtabs(DOMAIN_SIZE_PSU~strata.id+DOMAIN_ID,data=dem) #
		dem$DOMAIN_SIZE_STR<-demographSums[cbind(dem$strata.id,dem$DOMAIN_ID)] #Nk+l 
		dem$strata.n<-z$ST1_ANT_SSU[dem$strata.id]
		dem$time<-dem$id<-NULL
		row.names(dem)<-NULL
		PSU1<-match(names(x)[PSU1], names(dem))
		size1<-match(names(x)[size1], names(dem))
		strata1<-match(names(x)[strata1], names(dem))

# new domain allocation domainalloc = 0 for fixed equal allocation, 1 for proportional          ### new bit ###
		demostrat <- data.frame(demographSums)						### new bit ###
		if (domainalloc == 0) {								### new bit ###
			demostrat$prop <- 1/length(demo1)} else {				### new bit ###
			demostrat$prop <- demostrat$Freq / z[demostrat$strata.id, size2]	### new bit ###
			demostrat$Freq<-NULL							### new bit ###
		}
	}

#calculate n for those with inc==1 first 
	big<-x[which(x$ST1_PROB==1),]
	if (nrow(big)>0 & demoyesno=="0") {
		message <- paste("Sample sizes could not be calculated because the following PSU(s) have inclusion probabilities equal to 1: ", list(unique(big[,PSU1])))
		tkmessageBox(message=message)
		return()

	} else if (nrow(big)>0 & demoyesno=="1") {
		big$ST2_ANT_SSU<-big$strata.n*(big[,size1]/big$strata.pop)
		demfreq<-data.frame(xtabs(big$ST2_ANT_SSU~big$strata.id))
		z$minus.n[z$strata.id%in%demfreq[,1]]<-demfreq[,2]
		z$n.new<-ifelse(is.na(z$minus.n),z$ST1_ANT_SSU,z$ST1_ANT_SSU-z$minus.n)
		bigmk<-data.frame(table(big$strata.id))
		z$minus.mk[z$strata.id%in%bigmk[,1]]<-bigmk[,2]
		problem2<-ifelse(z$minus.mk>z$ST1_SEL_PSU,1,0)
		if(sum(problem2,na.rm=TRUE)>0) {
			tkmessageBox(message="Sample sizes have not been calculated because one (or more) stratum has too many PSU with inclusion probability of 1")
			tkfocus(CommanderWindow())
			}
		z$m.new<-ifelse(is.na(z$minus.mk),z$ST1_SEL_PSU,z$ST1_SEL_PSU-z$minus.mk)

#Join dem/x and big and add in other n's
		x<-merge(x,big,all=TRUE)
		z$DOMAINS<-1
		x$ST2_ANT_SSU<-ifelse(is.na(x$ST2_ANT_SSU),z$n.new[x$strata.id]/z$m.new[x$strata.id],x$ST2_ANT_SSU) #nkj
		x$ST2_PROB<-x$ST2_ANT_SSU/x[,size1]
#clean up
		z$minus.n<-z$n.new<-z$minus.mk<-z$m.new<-z$strata.pop.new<-NULL
		message("Sample sizes have been calculated, however one or more PSUs have inclusion probabilities equal to 1")


############# removed part including demo (add back in later)
#		if (demoyesno=="0"){
#			big<-dem[dem[,PSU1]%in%big[,PSU1],]
#			big <- merge(big, demostrat, by=c("strata.id", "DOMAIN_ID"))         ### new bit ###
#			big <- big[,c(3:ncol(big),1,2)]                                      ### new bit ###
#			big$ST2_ANT_SSU<-big$strata.n*(big[,size1]/big$strata.pop)*big$prop  ### new bit ###
#		} else {big$ST2_ANT_SSU<-big$strata.n*(big[,size1]/big$strata.pop)}

#re-adjust strata n
#		demfreq<-data.frame(xtabs(big$ST2_ANT_SSU~big$strata.id))
#		z$minus.n[z$strata.id%in%demfreq[,1]]<-demfreq[,2]
#		z$n.new<-ifelse(is.na(z$minus.n),z$ST1_ANT_SSU,z$ST1_ANT_SSU-z$minus.n)

#adjust m
#		if(demoyesno=="0"){
#			bigmk<-data.frame(table(big$strata.id)/length(demo1))  ##Need to change here demostrat$prop?
#			} else {bigmk<-data.frame(table(big$strata.id))}
#		z$minus.mk[z$strata.id%in%bigmk[,1]]<-bigmk[,2]
#		problem2<-ifelse(z$minus.mk>z$ST1_SEL_PSU,1,0)
#		if(sum(problem2,na.rm=TRUE)>0) {
#			tkmessageBox(message="Sample sizes have not been calculated because one (or more) stratum has too many PSU with inclusion probability of 1")
#			tkfocus(CommanderWindow())
#			}
#		z$m.new<-ifelse(is.na(z$minus.mk),z$ST1_SEL_PSU,z$ST1_SEL_PSU-z$minus.mk)
#	
#Join dem/x and big and add in other n's
#		if (demoyesno==0) {
#			x<-merge(dem,big,all=TRUE)
#			z$strata.pop.new<-z[,size2]-xtabs(x[,size1]~x[,strata1], subset=!is.na(x$ST2_ANT_SSU))/length(demo1)	
#			x$demo.n.new<-z$n.new[x$strata.id]/length(demo1) 							###Need to change here
#			x$demo.pop.remove<-(xtabs(x$DOMAIN_SIZE_PSU~x[,strata1]+x$DOMAIN_ID, subset=!is.na(x$ST2_ANT_SSU)))[cbind(x$strata.id,x$DOMAIN_ID)]
#			x$DOMAIN_SIZE_STR.new<-x$DOMAIN_SIZE_STR-x$demo.pop.remove
#			z$DOMAINS<-length(demo1) 							
#			x$ST2_ANT_SSU<-ifelse(is.na(x$ST2_ANT_SSU),(x$demo.n.new/z$m.new[x$strata.id])*((as.double(z$strata.pop.new[x$strata.id])*x$DOMAIN_SIZE_PSU)/(x$DOMAIN_SIZE_STR.new*x[,size1])),x$ST2_ANT_SSU)
#			x$ST2_PROB<-x$ST2_ANT_SSU/x$DOMAIN_SIZE_PSU
#			} else{
#			x<-merge(x,big,all=TRUE)
#			z$DOMAINS<-1
#			x$ST2_ANT_SSU<-ifelse(is.na(x$ST2_ANT_SSU),z$n.new[x$strata.id]/z$m.new[x$strata.id],x$ST2_ANT_SSU) #nkj
#			x$ST2_PROB<-x$ST2_ANT_SSU/x[,size1]
#			}
#
#Clean up
#		z$minus.n<-z$n.new<-z$minus.mk<-z$m.new<-z$strata.pop.new<-NULL
#		message("Sample sizes have been calculated, however one or more PSUs have inclusion probabilities equal to 1")
########################		
	
#Back to if there are no PSU with inc=1
	} else if (nrow(big)==0) {
		if (demoyesno==0) {
			x<-merge(dem, demostrat, by=c("strata.id", "DOMAIN_ID"))   ### new bit ###
			x <- x[,c(3:ncol(x),1,2)]                                  ### new bit ###
			x$demo.n<-z$ST1_ANT_SSU[x$strata.id]*x$prop    		   ### new bit ###
			x$ST2_ANT_SSU<-(x$demo.n/x$ST1_SEL_PSU)*((as.double(x$strata.pop)*x$DOMAIN_SIZE_PSU)/(x$DOMAIN_SIZE_STR*x[,size1]))
			x$ST2_PROB<-ifelse(x$ST2_ANT_SSU==0, 0, x$ST2_ANT_SSU/x$DOMAIN_SIZE_PSU)  
			z$DOMAINS<-length(demo1)                      
			}
		if (demoyesno==1) {
			x$ST2_ANT_SSU<-x$strata.n/x$ST1_SEL_PSU    #nkj
			x$ST2_PROB<-ifelse(x$ST2_ANT_SSU==0, 0, x$ST2_ANT_SSU/x[,size1]) 
			z$DOMAINS<-1                     
			} 
		}

#Add in extras and General Clean-up
	x$ALL_PROB<-x$ST1_PROB*x$ST2_PROB
	x$SAMPLING_WEIGHT<-1/x$ALL_PROB
	assign("sstotal", round(sum(z$ST1_ANT_SSU)), inherits=TRUE)
	if (mkvalue==1) z$ST1_COST <-(z[,PSUcost]*z$ST1_SEL_PSU)+(z[,Parcost]*z$ST1_ANT_SSU)
	if (mkvalue==1) z$STR_VARIANCE<-(z[,Vamong]/z$ST1_SEL_PSU)+(z[,Vwithin]/z$ST1_ANT_SSU)
	if (mkvalue==1 & demoyesno==1) {
		assign("TotalVar", sum((z$ww^2)*z$STR_VARIANCE), inherits=TRUE)
		assign("COEF", (sqrt(TotalVar)/T.Var)*100, inherits=TRUE)
		z$COEF<-(sqrt(z$STR_VARIANCE)/z$mean.Var)*100
		VarAmongStrat<-sum((z$ST1_VWITHIN+z$ST1_VAMONG)*z$ww)
		VarWithinStrat<-sum(z$ww*((z$mean.Var-T.Var)^2))
		assign("Deff", TotalVar/((VarAmongStrat+VarWithinStrat)/sum(z$ST1_ANT_SSU)), inherits=TRUE)
	}
	if (ppstype ==1 | ppstype==2) { z$ST1_CV <- 1 } else { z$ST1_CV <- 2 }
	z$strata.id<-x$strata.id<-x$strata.pop<-x$strata.n<-x$demo.n<-x$row.names<-z$ww<-z$strata.pop.new<-x$demo.n.new<-x$demo.pop.remove<-x$DOMAIN_SIZE_STR.new<-NULL
	x<-x[order(x[,match(PSUname, names(x))]),]
	row.names(x)<-c(seq(1:nrow(x)))


#return 2 datasets	
	utfall<-list(as.data.frame(z),as.data.frame(x))
	names(utfall)<-c("StrataSS","PSUSS")
	return(utfall)
}
