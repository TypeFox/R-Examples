# R program for LD_CI macro
#
#   Input:   
#               y: a vector of variable of interest
#               group: a vector of factor variable
#               time: a vector of the time variable
#               subject: a vector of independent subjects
#               alpha: significance level for the confidence interval
#               group.name: name of group vector
#               time.name: name of time vector
#               description: description of the output. Default is set to TRUE (show description)
#		time.order: a vector of time levels specifying the order.
#		group.order: a vector of group levels specifying the order.
#   Output:
#             confidence interval and bias estimation in a matrix form
#  
ld.ci<-function(y, time, subject, group=NULL, alpha=0.05, time.name="Time", group.name="Group", 
description=TRUE, time.order=NULL, group.order=NULL, rounds=4, plot.CI = TRUE, order.warning=TRUE)
{
#        For model description see Brunner et al. (2002)
#    
#        Author: Karthinathan Thangavelu (kthanga@gwdg.de)
#                     Department of Medical Statistics, Goettingen, Germany
#
#         Version:  01-01
#         Date: February 18, 2003
#
#
#        Editied by: Kimihiro Noguchi
#         Version:  01-02
#         Date: August 27, 2009
#
#        Editied by: Kimihiro Noguchi
#         Version:  01-03
#         Date: December 24, 2009
#
#    Key Variables:
#                timenames: names (or numbers) of the time
#                time: time factor
#                subject: subject factor
#                group: group factor
#                facs: levels of the group factor
#                A: number of levels of whole plot factor
#                T: number of levels of time factor
#                D: variable of interest
#                RD: ranks of the observations
#                N: total number of subject
#                Lamda: indicator of whether the data are present (0=missing, 1=observed)
#                Omat: Matrix of the observations
#                Rmat: Matrix of the ranks by different factors, where missing observations=0
#                Lamdamat: Matrix of the Lamda indicator by different times
#                sort1: sorted time factor levels
#                NN: total number of observations
#                RMeans: rank means of each unique factor combinations and each factor
#                Ni: number of observations at each level

     	   var<-y
	   if(is.null(var)||is.null(time)||is.null(subject)) 
		stop("At least one of the input parameters (y, time, or subject) is not found.")
	   
           sublen<-length(subject)
	   varlen<-length(var)
	   timlen<-length(time)
	   
	   if((sublen!=varlen)||(sublen!=timlen))
		stop("At least one of the input parameters (y, time, or subject) has a different length.")

	   if(!is.null(group))
	   {
		grolen<-length(group)
		if((sublen!=varlen)||(sublen!=grolen))
			stop("Length of the subject does not match the length of group")
	   }	

	if(is.null(group))
	{
 		group<-rep(1,length(subject))
	}

	### initialize variables

	glevel <- unique(group)
	tlevel <- unique(time)
	slevel <- unique(subject)
	t <- length(tlevel)
	s <- length(slevel)
	a <- length(glevel)

	if((t*s)!=length(var))
	stop("Number of levels of subject (",s, ") times number of levels of time (",t,") 
	is not equal to the total number of observations (",length(var),").",sep="")

	#    time order vector

	if(!is.null(time.order))
	{
		tlevel <- time.order
		tlevel2 <- unique(time)

		if(length(tlevel)!=length(tlevel2))	# if the levels of the order is different from the one in the data
		stop("Length of the time.order vector (",length(tlevel), ") 
		is not equal to the levels of time vector (",length(tlevel2),").",sep="")

		if(mean(sort(tlevel)==sort(tlevel2))!=1)     # if the elements in the time.order is different from the time levels
		stop("Elements in the time.order vector is different from the levels specified in the time vector.",sep="")		
	}

	#    group order vector

	 if(!is.null(group.order))
	 {
		glevel <- group.order
		glevel2 <- unique(group)

		if(length(glevel)!=length(glevel2))	# if the levels of the order is different from the one in the data
		stop("Length of the group.order vector (",length(glevel), ") 
		is not equal to the levels of group vector (",length(glevel2),").",sep="")

		if(mean(sort(glevel)==sort(glevel2))!=1)     # if the elements in the group.order is different from the group levels
		stop("Elements in the group.order vector is different from the levels specified in the group vector.",sep="")		
	 }

	#    sort data

	sortvector<-double(length(var))
	newtime<-double(length(var))
	newsubject<-double(length(var))
	newgroup<-double(length(var))

	for(i in 1:length(var))
        {
           row<-which(subject[i]==slevel)
	   col<-which(time[i]==tlevel)
	   newsubject[i]<-row
	   newtime[i]<-col
	   newgroup[i]<-which(group[i]==glevel)
	   sortvector[((col-1)*s+row)]<-i
        }   
	
	subject<-newsubject[sortvector]
	var<-var[sortvector]
	time<-newtime[sortvector]
	group<-newgroup[sortvector]

	#    sort again by group, and assign new subject numbers to subjects

	grouptemp<-order(group[1:s])
	groupplus<-(rep(c(0:(t-1)),e=s))*s
	groupsort<-(rep(grouptemp,t))+groupplus
	
	subject<-rep(c(1:s),t)
	var<-var[groupsort]
	time<-time[groupsort]
	group<-group[groupsort]

	#    organize data

	timenames<-unique(time)
	time<-factor(time)
	subject<-factor(subject)
	group<-factor(group)
	facs<-levels(group)
	A<-nlevels(group)
	T<-nlevels(time)
	D<-var
	RD <-  rank(D)
	N <- length(unique(subject)) 
	Lamda <-  1-as.numeric(is.na(D))

	### the program is not executed if there is a missing observation	

	if(is.element(0, Lamda)) 
	{
	cat("\n There are some missing in the data. This program does not ")
		cat("\n accept missing values. Please verify and try again.\n")
		return(" Status : Stopped execution of program due to error in input. ")
	}

	### transform the data vectors into matrix

	Omat <- matrix(NA, N, T)
	Rmat <-  matrix(NA, N, T)
	Lamdamat <-  matrix(NA, N, T)
	sort1<-unique(time)

	for(i in 1:T)
	{
		tempD<-D[which(time==sort1[i])] 
 		tempRD<-RD[which(time==sort1[i])]
 		tempLM<-Lamda[which(time==sort1[i])]
 		tempSJ<-subject[which(time==sort1[i])]
 		orderSJ<-order(tempSJ)
 		tempD<-tempD[orderSJ]
 		tempRD<-tempRD[orderSJ]
 		tempLM<-tempLM[orderSJ]

 		if(i==1)
 		{
  			tempFC<-group[which(time==sort1[i])]
  			FAC<-tempFC[orderSJ]
			FAC<-factor(FAC)
 		}

 		Omat[,i]<-tempD
 		Rmat[,i]<-tempRD
 		Lamdamat[,i]<-tempLM
	}

	D<-Omat
	NN <-  sum(Lamda)
	Rmat <-  Rmat * Lamdamat 
	RMeans <-  rep(0, A*T)
	Ni <-  apply(Lamdamat, 2, sum)
	Rmat.fac <-  cbind(FAC, Rmat)
	colnames(Rmat.fac) <-  NULL
	Lamdamat.fac <-  cbind(FAC, Lamdamat) 
	colnames(Lamdamat.fac) <-  NULL

	### calculate rank means

	fn.fact.manip <-  function(fullRmat, n, n.a, n.t, A.ni) 
	{
		res.AT <-  list(0)

		for(i in 1:n.a) 
			res.AT[[i]] <-  fullRmat[(fullRmat[,1]==i),-1]
	
		return(res.AT)
	}

	Ra.list <-  fn.fact.manip(Rmat.fac, N, A, T, as.vector(summary(FAC)))
	R1a.list <-  fn.fact.manip(Lamdamat.fac, N, A, T, as.vector(summary(FAC)))

	for(i in 1:A) 
		RMeans[((i-1)*T + 1): ((i-1)*T + T)] <-  apply(Ra.list[[i]], 2, sum) / apply(
		R1a.list[[i]], 2, sum) 

	RTE <- (RMeans - 0.5) / NN 
	GROUP <-  rep(glevel, e=T)
	TIME <-  rep(tlevel,A)
	Nobs <-  rep(0, A*T)

	for(i in 1:A) 
		Nobs[((i-1)*T + 1): ((i-1)*T + T)] <-  apply(R1a.list[[i]], 2, sum) 

	PRes1 <- cbind(RMeans, Nobs, RTE)
###########################################################




	rd.PRes1 <- round(PRes1, Inf)
	rd.PRes1 <- cbind(TIME, rd.PRes1)
	D.fac <-  cbind(FAC, D)

	fn.compli.sig <-  function(fullDmat, n, nn, n.a, n.t, FAC, A.ni) 
	{
		res.list <-  list(0)
		res.R.is <-  list(0)
		res.R.mis <-  list(0)
		res.bias <-  rep(0, 0)
		MAXI <-  max(fullDmat)

		for(i in 1:n.a) 
		{
			res.R.is[[i]] <-  list(0)
			res.R.mis[[i]] <-  list(0)
			temp.Dmat <-  fullDmat[(fullDmat[,1]==i), -1]
			temp.mat <-  matrix(0, A.ni[i], n.t)
			temp.Rmat <-  matrix(0, A.ni[i], n.t)

			for(j in 1:A.ni[i]) 
				temp.mat[j,] <-  rank(temp.Dmat[j,])

			temp.mat <-  temp.mat * matrix(1-as.numeric(is.na(c(temp.mat))), A.ni[i], n.t)
			temp.Rmat <-  matrix(rank(temp.Dmat), A.ni[i], n.t) * matrix(1-as.numeric(is.na(c(
			temp.Dmat))), A.ni[i], n.t)

			for(s in 1:n.t) 
			{
				res.R.is[[i]][[s]] <-  rep(0, 0)
				res.R.mis[[i]][[s]] <-  matrix(0, n, n.t+1)

				temp.Dmat <-  fullDmat

				t.r.is <-  rank(temp.Dmat[(temp.Dmat[,1]==i), -1][,s])
				t.l.is <-  1-as.numeric(is.na(c(temp.Dmat[(temp.Dmat[,1]==i), -1][,s])))
				t.r.is <-  t.r.is * t.l.is

				temp.Dmat[(temp.Dmat[,1]==i), -1][,s] <-  MAXI + 1
				t1.Rmat <-  matrix(rank(temp.Dmat[,-1]), n, n.t)
				t1.Rmat <-  t1.Rmat * matrix(1-as.numeric(is.na(c(temp.Dmat[,-1]))), n, n.t)
				t1.Rmat[(temp.Dmat[,1]==i), s] <-  0
				res.R.is[[i]][[s]] <-  t.r.is
				res.R.mis[[i]][[s]] <-  cbind(FAC, t1.Rmat)
				colnames(res.R.mis[[i]][[s]]) <-  NULL
				temp.bias <-  (1/A.ni[i])*(sum(temp.mat[,s])) - 0.5
				temp.bias <-  temp.bias - ((1/A.ni[i]) * ((1/A.ni[i]) * sum(temp.Rmat[,s]) - 0.5))
				temp.bias <-  (A.ni[i]/(nn * (A.ni[i] - 1))) * temp.bias 
				res.bias <-  c(res.bias, temp.bias)
			}
		}

		res.list[[1]] <-  res.R.is
		res.list[[2]] <-  res.R.mis
		res.list[[3]] <-  res.bias

		return(res.list)
	}


	A.ni <-  as.vector(summary(FAC))
	Da.list <-  fn.compli.sig(D.fac, N, NN, A, T, FAC, A.ni)
	Y.list <-  list(0)

	for(i in 1:A) 
	{
		Y.list[[i]] <-  matrix(0, A.ni[i], T)

		for(k in 1:A.ni[i]) 
		{
			for(s in 1:T) 
			{
				Y.list[[i]][k, s] <-  (1/NN) * ((2*Rmat.fac[(Rmat.fac[,1]==i),-1][k, s]) - 
				Da.list[[1]][[i]][[s]][k] - sum(Rmat.fac[(Rmat.fac[,1]==i),-1][k, ]) + sum(
				Da.list[[2]][[i]][[s]][(Da.list[[2]][[i]][[s]][,1]==i), -1][k,]))
			}
		}
	}

	Z.list <-  list(0)
	tau.u.is.sq <-  list(0)

	for(i in 1:A) 
	{
		Z.list[[i]] <-  list(0)
		tau.u.is.sq[[i]] <-  list(0)

		for(s in 1:T) 
		{
			Z.list[[i]][[s]] <-  list(0)
			tau.u.is.sq[[i]][[s]] <-  rep(0, A)

			for(u in 1:A) 
			{
				Z.list[[i]][[s]][[u]] <-  rep(0, A.ni[u])

				if(u != i) 
				{
					for(k in 1:A.ni[u]) 
					{
						Z.list[[i]][[s]][[u]][k] <-  (1/A.ni[i]) * 
						(sum(Rmat.fac[(Rmat.fac[,1]==u),-1][k, ]) - 
						sum(Da.list[[2]][[i]][[s]][(Da.list[[2]][[i]][[s]][,1]==u), -1][k,]))
					}

					tau.u.is.sq[[i]][[s]][u] <-  (1/(A.ni[u] - 1)) * sum((Z.list[[i]][[s]][[u]] 
					- mean(Z.list[[i]][[s]][[u]])) * (Z.list[[i]][[s]][[u]] - 
					mean(Z.list[[i]][[s]][[u]])))
				}
			}
		}
	}

	lam.is.sq <-  matrix(0, A, T)

	for(i in 1:A) 
	{
		for(s in 1:T) 
		{
			lam.is.sq[i, s] <-  (1/(A.ni[i] - 1)) * (sum((Y.list[[i]][,s] - mean(Y.list[[i]][,s])) * 
			(Y.list[[i]][,s] - mean(Y.list[[i]][,s]))))
		}
	}

	sig.sq <-  matrix(0, A, T)

	for(i in 1:A) 
	{
		for(s in 1:T) 
		{
			sig.sq[i, s] <-  ((N/A.ni[i]) * lam.is.sq[i, s]) + ((1/(N*T*T)) * sum(A.ni * 
			tau.u.is.sq[[i]][[s]]))
		}
	}

	P.vec <-  t(matrix(RTE, T, A))
	pL <-  matrix(0, A, T)
	pU <-  matrix(0, A, T)

	### calculate coonfidence intervals

	for(i in 1:A) 
	{
		for(s in 1:T) 
		{
			pg <-  ((NN * P.vec[i, s]) - A.ni[i]/2) / (NN-A.ni[i])
			sd <-  (NN * sqrt(sig.sq[i, s])) / (NN-A.ni[i])
			t1 <-  log(pg/(1-pg))
			t2 <-  (sd/ (pg * (1 - pg) * sqrt(N))) * qnorm((1-(alpha/2)))
			pgL <-  t1 - t2
			pgU <-  t1 + t2
			pL[i, s] <-  (A.ni[i]/(2*NN)) + ((NN - A.ni[i]) / NN) * (exp(pgL) / (1 + exp(pgL)))
			pU[i, s] <-  (A.ni[i]/(2*NN)) + ((NN - A.ni[i]) / NN) * (exp(pgU) / (1 + exp(pgU)))

		}
	}

	CI <-  cbind(matrix(t(pL), A*T, 1), matrix(t(pU), A*T, 1))

#-------------------------Create a nice output table---------------------------#
frank.RMeans <- round(RMeans,rounds)
frank.Nobs <- Nobs
frank.RTE <- round(RTE,rounds)
frank.Bias <- round(c(Da.list[[3]]),rounds)
frank.Var <-round(c(t(sig.sq)),rounds)
frank.Lower <- round(matrix(t(pL), A*T, 1),rounds)
frank.Upper <- round(matrix(t(pU), A*T, 1),rounds)

frank.Time<-TIME
frank.Group<-paste(group.name,GROUP,sep="")


	rd.PRes1 <- cbind(rd.PRes1, round(Da.list[[3]], 4), round(c(t(sig.sq)), 4), round(CI, 4))
	colnames(rd.PRes1) <-  c(time.name, "RankMeans", "Nobs", "RTE", "Bias", "Variance", "Lower_bound", "Upper_bound")
	rownames(rd.PRes1)<-paste(group.name,GROUP,sep="")

	### output descriptions

	if(description==TRUE)
	{
		cat(" Bias Estimation and Confidence Intervals for Relative Effects \n")
		cat("\n Total number of observations :", (N*T))
		cat("\n Total number of Subject :", N)
		cat("\n Significance level (alpha) : ", alpha)
		cat("\n\n Class level information ")
		cat("\n ----------------------- ")
		cat("\n Levels of sub-plot factor time : ", T)
		cat("\n Levels of whole-plot factor group : ", A)
		cat("\n\n Abbreviations ")
		cat("\n ----------------------- ")
           	cat("\n RankMeans = Rank means\n")
		cat("\n Nobs = Number of observations")
		cat("\n RTE = Relative treatment effect")
		cat("\n Bias = Bias estimation")
		cat("\n Variance = Variance estimation")
		cat("\n Lower_bound = Lower bound of the confidence interval")
		cat("\n Upper_bound = Upper bound of the confidence interval")
           	cat("\n Note: The description output above will disappear by setting description=FALSE in the input. See the help file for details.")
		cat("\n\n")
	}
	
	if(order.warning==TRUE)
	{
		cat(" LD CI Calculations ")
		cat("\n ----------------------- \n")
           	cat(" Order of the time and group levels.\n") 
           	cat(" Time level:  " , paste(tlevel),"\n")
           	cat(" Group level:  " , paste(glevel),"\n")
           	cat(" The order may be specified in time.order or group.order (does not affect the calculation).\n\n")
	}

#abline(h = 0, col = "red", lty = 1, lwd = 2)
 #for (i in 1:nc) {points(rep(k[i],2),c(Lower.SCI[i],Upper.SCI[i]),type="l")}
#upper<-"_"
#lower <- "_"
#points(x = k, y = Lower.SCI, pch = lower,cex=1.5)
#points(x = k, y = Upper.SCI, pch = upper,cex=1.5)
#axis(1, at = k, labels = rownames(C),font.axis=1.5,cex.axis=1.5)

#}


Frank <-data.frame(Group=frank.Group, Time = frank.Time, Nobs=frank.Nobs, RankMeans=frank.RMeans, RTE=frank.RTE,Bias=frank.Bias, Variance=frank.Var, Lower=frank.Lower, Upper = frank.Upper )


Frank
}

