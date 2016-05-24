# R program for F1_LD_F2 macro
#
#   Input:   
#               y: a vector of variable of interest
#               group: a vector of group variable
#               time1: a vector of the time1 variable
#               time2: a vector of the time2 variable
#               subject: a vector of independent subjects
#               time1.name: time1 factor name. Default is set to TimeC
#               time2.name: time2 factor name. Default is set to TimeT
#               group.name: whole plot factor names. Default is set to GroupA
#               description: description of the output. Default is set to TRUE (show description)
#		time1.order: a vector of time1 levels specifying the order.
#		time2.order: a vector of time2 levels specifying the order.
#		group.order: a vector of group levels specifying the order.
#   Output:
#             list of relative treatment effects, test results, covariance matrix
#  
f1.ld.f2 <- function(y, time1, time2, group, subject, time1.name="Time1",
time2.name="Time2", group.name="Group", description=TRUE, time1.order=NULL,
time2.order=NULL, group.order=NULL,plot.RTE=TRUE,show.covariance=FALSE,
order.warning=TRUE)
{
#        For model description see Brunner et al. (2002)
#    
#        Author: Karthinathan Thangavelu (kthanga@gwdg.de)
#                     Department of Medical Statistics, Goettingen, Germany
#
#         Version:  01-01
#         Date: March 1, 2003
#
#
#        Editied by: Kimihiro Noguchi
#         Version:  01-02
#         Date: August 25, 2009
#
#        Editied by: Kimihiro Noguchi
#         Version:  01-03
#         Date: December 24, 2009
#
#    Key Variables:
#                FAC: whole plot factor
#                time1: time1 factor
#                time2: time2 factor
#                subject: subject factor
#                D: variable of interest
#                Lamda: indicator of whether the data are present (0=missing, 1=observed)
#                levs: whole plot factor levels
#                A: number of levels of whole plot factor
#                N: total number of subject
#                C: number of levels of time1 factor
#                T: number of levels of time2 factor
#                CTcount: number of different combinations of the time factor levels
#                RD: ranks of the observations
#                uvector: number assigned to each unique combination of the factor levels
#                sort1: sorted time1 factor levels
#                sort2: sorted time2 factor levels
#                Rmat: Matrix of the observations by different factors, where missing observations=0
#                Lamdamat: Matrix of the Lamda indicator by different times
#                NN: total number of observations
#                RMeans: rank means of each unique factor combinations and each factor
#                Ni: number of observations at each level

#    check whether the input variables are entered correctly

	   var<-y
	   if(is.null(var)||is.null(time1)||is.null(time2)||is.null(group)||is.null(subject)) 
		stop("At least one of the input parameters (y, time1, time2, group, or subject) is not found.")
	   
           sublen<-length(subject)
	   varlen<-length(var)
	   tim1len<-length(time1)
	   tim2len<-length(time2)
	   grolen<-length(group)
	   
	   if((sublen!=varlen)||(sublen!=tim1len)||(sublen!=tim2len)||(sublen!=grolen))
		stop("At least one of the input parameters (y, time1, time2, group, or subject) has a different length.")


	library(MASS)

	sort1<-unique(time1)
	sort2<-unique(time2)
	sortg<-unique(group)
	sorts<-unique(subject)
	C <- length(sort1)
	T <- length(sort2)
	A <- length(sortg)
	N <- length(sorts)
	CTcount <- C*T
	uvector<-double(length(var))

	if((C*T*N)!=length(var))
	stop("Number of levels of subject (",N, ") times number of levels of time1 (",C,") times number of levels of time2 (",T,") 
	is not equal to the total number of observations (",length(var),").",sep="")

#    time and group order vectors

	   if(!is.null(time1.order))
	   {
		sort1 <- time1.order
		sort12 <- unique(time1)

		if(length(sort1)!=length(sort12))	# if the levels of the order is different from the one in the data
		stop("Length of the time1.order vector (",length(sort1), ") 
		is not equal to the levels of time1 vector (",length(sort12),").",sep="")

		if(mean(sort(sort1)==sort(sort12))!=1)     # if the elements in the time1.order is different from the time1 levels
		stop("Elements in the time1.order vector is different from the levels specified in the time1 vector.",sep="")		
	   }

	   if(!is.null(time2.order))
	   {
		sort2 <- time2.order
		sort22 <- unique(time2)

		if(length(sort2)!=length(sort22))	# if the levels of the order is different from the one in the data
		stop("Length of the time2.order vector (",length(sort2), ") 
		is not equal to the levels of time2 vector (",length(sort22),").",sep="")

		if(mean(sort(sort2)==sort(sort22))!=1)     # if the elements in the time2.order is different from the time levels
		stop("Elements in the time2.order vector is different from the levels specified in the time2 vector.",sep="")		
	   }

	   if(!is.null(group.order))
	   {
		sortg <- group.order
		sortg2 <- unique(group)

		if(length(sortg)!=length(sortg2))	# if the levels of the order is different from the one in the data
		stop("Length of the group.order vector (",length(sortg), ") 
		is not equal to the levels of group vector (",length(sortg2),").",sep="")

		if(mean(sort(sortg)==sort(sortg2))!=1)     # if the elements in the group.order is different from the group levels
		stop("Elements in the group.order vector is different from the levels specified in the group vector.",sep="")		
	   }

#    sort data

	newtime1<-double(length(var))
	newtime2<-double(length(var))
 
	for(i in 1:length(var))
	{
          uvector[i]<-T*which(sort1==time1[i])+which(sort2==time2[i])-T
	  newtime1[i]<-which(sort1==time1[i])
	  newtime2[i]<-which(sort2==time2[i])
        }  
 
	sortvector<-double(length(var))
	newsubject<-double(length(var))
	newgroup<-double(length(var))

	for(i in 1:length(var))
        {
          row<-which(subject[i]==sorts)
	  col<-uvector[i]
	  newsubject[i]<-row
	  newgroup[i]<-which(group[i]==sortg)
	  sortvector[((col-1)*N+row)]<-i
        }   

#    relabel the time and subject vectors (to deal with factor variables)

	subject<-newsubject[sortvector]
	var<-var[sortvector]
	time1<-newtime1[sortvector]
	time2<-newtime2[sortvector]
	group<-newgroup[sortvector]

#    sort again by group, and assign new subject numbers to subjects

	grouptemp<-order(group[1:N])
	groupplus<-(rep(c(0:(CTcount-1)),e=N))*N
	groupsort<-(rep(grouptemp,CTcount))+groupplus
	
	subject<-rep(c(1:N),CTcount)
	var<-var[groupsort]
	time1<-time1[groupsort]
	time2<-time2[groupsort]
	group<-group[groupsort]

#    pass the original sort variables to origsort, and assign new ones

	origsort1<-sort1
	origsort2<-sort2
	origsortg<-sortg
	origsorts<-sorts

	sort1<-unique(time1)
	sort2<-unique(time2)
	sortg<-unique(group)
	sorts<-unique(subject)

#    organize variables

	factor<-group
	factor.name<-group.name
	FAC<-factor(factor)
	time1<-factor(time1)
	time2<-factor(time2)
	subject<-factor(subject)
	D<-var
	Lamda <- 1-as.numeric(is.na(D)) 
	levs <- origsortg
	A <- nlevels(FAC)
	N <- nlevels(subject)
	C <- nlevels(time1)
	T <- nlevels(time2)
	CTcount <- C*T
	RD <- rank(D)
	uvector <- double(N)
	Rmat <- matrix(NA, N, CTcount)
	Lamdamat <- matrix(NA, N, CTcount)

	for(i in 1:(N*C*T))
	{
	 	uvector[i] <- T*which(sort1==time1[i])+which(sort2==time2[i])-T
 		Rmat[subject[i],uvector[i]]<-RD[i]
 		Lamdamat[subject[i],uvector[i]]<-Lamda[i]
	}

	NN <- sum(Lamda) 
	Rmat <- Rmat * Lamdamat 
	RMeans <- rep(0, (A + C + T + (A*C) + (C*T) + (T*A) + (A*C*T)))

	### This order is important ie, A, C, T, (A*C), (C*T), (T*A), (A*C*T) 

	Ni <- apply(Lamdamat, 2, sum) 

	fn.collapse.day <- function(mat, n, n.c, n.t) 
	{
		res <- matrix(0, nrow=n.t*n, ncol=n.c)
		for(i in 1:n.c) 
		{
			for(j in 1:n.t) 
			{
				res[(((j-1)*n+1):((j-1)*n+n)), i] <- mat[, ((i-1)*n.t + j)]
			}
		}
		return(res)
	}

	Rr <- fn.collapse.day(Rmat, N, C, T)
	R1r <- fn.collapse.day(Lamdamat, N, C, T)

	fn.collapse.time <- function(mat, n, n.c, n.t) 
	{
		res <- matrix(0, nrow=n.c*n, ncol=n.t)
		for(i in 1:n.t) 
		{
			for(j in 1:n.c) 
			{
				res[(((j-1)*n+1):((j-1)*n+n)), i] <- mat[, ((j-1)*n.t + i)]
			}
		}
		return(res)
	}
	
	Rs <- fn.collapse.time(Rmat, N, C, T)
	R1s <- fn.collapse.time(Lamdamat, N, C, T)
	
	n.a.vec<-tapply(factor,factor,length)/CTcount
	
	fn.fact.manip <- function(fullRmat, n, n.a.vec, n.c, n.t) 
	{
 	 res.list <- list(0)
 	 n.a<-length(n.a.vec)
 	 res.A <- list(0)
 	 res.AC <- list(0)
 	 res.AT <- list(0)
 	 res.ACT <- list(0)
 	 
 	 n.a.start<-c(1,(cumsum(n.a.vec)+1))[-(length(n.a.vec)+1)]
 	 n.a.end<-cumsum(n.a.vec)
	 
	 for(i in 1:n.a) 
	 {
	  n.a.len<-n.a.vec[i]
	  temp.Rmat <- fullRmat[((n.a.start[i]):(n.a.end[i])), ]
	  res.A[[i]] <- c(temp.Rmat)
	  res.AC[[i]] <- fn.collapse.day(temp.Rmat, n.a.len, n.c, n.t)
	  res.AT[[i]] <- fn.collapse.time(temp.Rmat, n.a.len, n.c, n.t)
	  res.ACT[[i]] <- temp.Rmat
	 }
	 
	 res.list[[1]] <- res.A
	 res.list[[2]] <- res.AC
	 res.list[[3]] <- res.AT
	 res.list[[4]] <- res.ACT
	 return(res.list)
	}

	Ra.list <- fn.fact.manip(Rmat, N, n.a.vec, C, T)
	R1a.list <- fn.fact.manip(Lamdamat, N, n.a.vec, C, T)

	RMeans[1 : A] <- unlist(lapply(Ra.list[[1]], sum)) / unlist(lapply(R1a.list[[1]], sum))

	RMeans[(A + 1) : (A + C)] <- (apply(Rr, 2, sum) / apply(R1r, 2, sum))
	RMeans[(A + C + 1) : (A + C + T)] <- (apply(Rs, 2, sum) / apply(R1s, 2, sum))

	RMeans[(A + C + T + 1) : (A + C + T + C*A)] <- 
	unlist(lapply(Ra.list[[2]],apply,MARGIN=2,sum)) / unlist(lapply(R1a.list[[2]],apply,MARGIN=2,sum))


	RMeans[(A + C + T + C*A + 1) : (A + C + T + C*A + C*T)] <- 
	(apply(Rmat, 2, sum) / apply(Lamdamat, 2, sum))

	RMeans[(A + C + T + C*A + C*T + 1) : (A + C + T + C*A + C*T + T*A)] <- 
	unlist(lapply(Ra.list[[3]],apply,MARGIN=2,sum)) / unlist(lapply(R1a.list[[3]],apply,MARGIN=2,sum))

	RMeans[(A + C + T + C*A + C*T + T*A + 1) : (A + C + T + C*A + C*T + T*A + A*C*T)] <- 
	unlist(lapply(Ra.list[[4]],apply,MARGIN=2,sum)) / unlist(lapply(R1a.list[[4]],apply,MARGIN=2,sum))

	RTE <- (RMeans - 0.5) / NN 

	time1.vec <- c(paste(time1.name, origsort1, sep=""))
	time2.vec <- c(paste(time2.name, origsort2, sep=""))

	fn.nice.out <- function(A, C, T) 
	{
		SOURCE <- rep(0,0)

		for(i in 1:A)
			SOURCE <- c(SOURCE, paste(factor.name, levs[i], sep=""))
		for(i in 1:C)
			SOURCE <- c(SOURCE, paste(time1.vec[i], sep=""))
		for(i in 1:T)
			SOURCE <- c(SOURCE, paste(time2.vec[i], sep=""))
		for(i in 1:A)
			for(j in 1:C)
				SOURCE <- c(SOURCE, paste(factor.name, levs[i], ":", time1.vec[j], sep=""))
		for(i in 1:C)
			for(j in 1:T)
				SOURCE <- c(SOURCE, paste(time1.vec[i], ":", time2.vec[j], sep=""))
		for(i in 1:A)
			for(j in 1:T)
				SOURCE <- c(SOURCE, paste(factor.name, levs[i], ":", time2.vec[j], sep=""))
		for(i in 1:A)
			for(j in 1:C)
				for(k in 1:T)
					SOURCE <- c(SOURCE, paste(factor.name, levs[i], ":", time1.vec[j], ":", time2.vec[k], sep=""))

		return(SOURCE)
	}

	SOURCE <- fn.nice.out(A, C, T)

	Nobs <- rep(0, (A + C + T + (A*C) + (C*T) + (T*A) + (A*C*T)))
	Nobs[1 : A] <- unlist(lapply(R1a.list[[1]], sum))
	Nobs[(A + 1) : (A + C)] <- (apply(R1r, 2, sum))
	Nobs[(A + C + 1) : (A + C + T)] <- (apply(R1s, 2, sum))

	Nobs[(A + C + T + 1) : (A + C + T + C*A)] <- 
	unlist(lapply(R1a.list[[2]],apply,MARGIN=2,sum))

	Nobs[(A + C + T + C*A + 1) : (A + C + T + C*A + C*T)] <- (apply(Lamdamat, 2, sum))

	Nobs[(A + C + T + C*A + C*T + 1) : (A + C + T + C*A + C*T + T*A)] <- 
	unlist(lapply(R1a.list[[3]],apply,MARGIN=2,sum))

	Nobs[(A + C + T + C*A + C*T + T*A + 1) : (A + C + T + C*A + C*T + T*A + A*C*T)] <- 
	unlist(lapply(R1a.list[[4]],apply,MARGIN=2,sum))

	PRes1 <- data.frame(RankMeans=RMeans, Nobs, RTE) 
	rd.PRes1 <- round(PRes1, Inf)
	rownames(rd.PRes1)<-SOURCE

	   model.name<-"F1 LD F2 Model"
     	if(description==TRUE)
     	{
	   cat("\n Total number of observations : ", NN)
	   cat("\n Total Number of subjects ", N)
	   cat("\n Total Number of missing observations : ", (N*CTcount - NN), "\n")
           cat("\n Class level information ")
           cat("\n ----------------------- ")
	   cat("\n Levels of", factor.name, "(whole-plot factor group) : ", A)
	   cat("\n Levels of", time1.name, "(sub-plot factor time1) : ", C)
	   cat("\n Levels of", time2.name, "(sub-plot factor time2) : ", T,"\n")
           cat("\n Abbreviations ")
           cat("\n ----------------------- \n")
           cat(" RankMeans = Rank means\n")
           cat(" Nobs = Number of observations\n")
           cat(" RTE = Relative treatment effect\n")
           cat(" Wald.test = Wald-type test statistic\n")
           cat(" ANOVA.test = ANOVA-type test statistic\n")
           cat(" ANOVA.test.mod.Box = modified ANOVA-type test statistic with Box approximation\n")
           cat(" covariance = Covariance matrix","\n")
           cat(" Note: The description output above will disappear by setting description=FALSE in the input. See the help file for details.","\n\n")
      	}

	if(order.warning==TRUE)
	{
           cat(" F1 LD F2 Model ") 
           cat("\n ----------------------- \n")
           cat(" Check that the order of the time1, time2, and group levels are correct.\n") 
           cat(" Time1 level:  " , paste(origsort1),"\n")
           cat(" Time2 level:  " , paste(origsort2),"\n")
           cat(" Group level:  " , paste(origsortg),"\n")
           cat(" If the order is not correct, specify the correct order in time1.order, time2.order, or group.order.\n\n")
	}

	fn.P.mat <- function(arg1) 
	{
		I <- diag(1, arg1, arg1)
		J <- matrix((1/arg1), arg1, arg1)
		return(I - J)
	}


	PA <- fn.P.mat(A)
	PC <- fn.P.mat(C)
	PT <- fn.P.mat(T)

	A1 <- matrix((1/A), 1, A)
	C1 <- matrix((1/C), 1, C)
	T1 <- matrix((1/T), 1, T)

	CA <- kronecker(PA, kronecker(C1, T1))
	CC <- kronecker(A1, kronecker(PC, T1))
	CT <- kronecker(A1, kronecker(C1, PT))
	CAC <- kronecker(PA, kronecker(PC, T1))
	CCT <- kronecker(A1, kronecker(PC, PT))
	CAT <- kronecker(PA, kronecker(C1, PT))
	CACT <- kronecker(PA, kronecker(PC, PT))

	fn.covr<-function(N,d,NN,Rmat,DatRMeans,Lamdamat, n.a.vec, C, T)
	{
		V<-matrix(0,d,d);
 	 	n.a.start<-c(1,(cumsum(n.a.vec)+1))[-(length(n.a.vec)+1)]
 	 	n.a.end<-cumsum(n.a.vec)
		
		fn.covr.block.mats <- function(N,d,NN,Ni,Rmat,DatRMeans,Lamdamat) 
		{
			V<-matrix(0,d,d)
			for(s in 1:d) 
			{
	   			for(sdash in 1:d) 
				{
	      				if(s==sdash) 
					{	    
		    				temp<-(Rmat[,s]-DatRMeans[s])*(Rmat[,s]-DatRMeans[s])
		    				V[s,sdash]<-V[s,sdash]+N*(Lamdamat[,s]%*%temp)/(NN^2*Ni[s]*(Ni[s]-1))
					}
	
	      				if(s!=sdash)
	       				{
		   				temp<-(Rmat[,s]-DatRMeans[s])*(Rmat[,sdash]-DatRMeans[sdash])
		   				temp1<-Lamdamat[,s]*Lamdamat[,sdash]
		   				ks<-(Ni[s]-1)*(Ni[sdash]-1)+Lamdamat[,s]%*%Lamdamat[,sdash]-1
		   				V[s,sdash]<-V[s,sdash]+N*(temp1%*%temp)/(NN^2*ks)
	       				} 	
	   			}
			}
	
			return(V);  
		}


		for(i in 1:(length(n.a.vec))) 
		{

			Ni <- apply(Lamdamat[((n.a.start[i]):(n.a.end[i])), (1: (C*T))], 2, sum)
			temp.mat <- fn.covr.block.mats((n.a.vec[i]), C*T, NN, Ni, Rmat[((n.a.start[i]):(n.a.end[i])), 
			(1: (C*T))],  DatRMeans[(((i-1)*(C*T) + 1) : ((i-1)*
			(C*T) + C*T))], Lamdamat[((n.a.start[i]):(n.a.end[i])), (1: (C*T))])

			V[(((i-1)*(C*T) + 1) : ((i-1)*(C*T) + C*T)), (((i-1)*(C*T) + 1) : ((i-1)*
			(C*T) + C*T))]  <-  (sum(n.a.vec)/(n.a.vec[i]))*temp.mat
		}

		return(V);

	}


	V <- fn.covr(N, A*C*T, NN, Rmat, RMeans[(A + C + T + C*A + C*T + T*A + 1) : 
	(A + C + T + C*A + C*T + T*A + A*C*T)],  Lamdamat, n.a.vec, C, T)

	SING.COV <- FALSE
	if(qr(V)$rank < (A*C*T)) SING.COV <- TRUE

	pvec <- RTE[(A + C + T + C*A + C*T + T*A + 1) : (A + C + T + C*A + C*T + T*A + A*C*T)]

	if(((C > 1) && (T > 1))) 
	{

		### Wald type statistics computed here

		WA <- N*t(CA%*%pvec)%*%ginv(CA%*%V%*%t(CA))%*%(CA%*%pvec)
		WC <- N*t(CC%*%pvec)%*%ginv(CC%*%V%*%t(CC))%*%(CC%*%pvec)
		WT <- N*t(CT%*%pvec)%*%ginv(CT%*%V%*%t(CT))%*%(CT%*%pvec)
		WAC <- N*t(CAC%*%pvec)%*%ginv(CAC%*%V%*%t(CAC))%*%(CAC%*%pvec)
		WCT <- N*t(CCT%*%pvec)%*%ginv(CCT%*%V%*%t(CCT))%*%(CCT%*%pvec)
		WAT <- N*t(CAT%*%pvec)%*%ginv(CAT%*%V%*%t(CAT))%*%(CAT%*%pvec)
		WACT <- N*t(CACT%*%pvec)%*%ginv(CACT%*%V%*%t(CACT))%*%(CACT%*%pvec)

		dfWA <- qr(CA%*%V%*%t(CA))$rank
		dfWC <- qr(CC%*%V%*%t(CC))$rank
		dfWT <- qr(CT%*%V%*%t(CT))$rank
		dfWAC <- qr(CAC%*%V%*%t(CAC))$rank
		dfWCT <- qr(CCT%*%V%*%t(CCT))$rank
		dfWAT <- qr(CAT%*%V%*%t(CAT))$rank	
		dfWACT <- qr(CACT%*%V%*%t(CACT))$rank

		if(!is.na(WA) && WA > 0) pWA <- pchisq(WA, dfWA,lower.tail=FALSE)
		else pWA <- NA
		if(!is.na(WC) && WC > 0) pWC <- pchisq(WC, dfWC,lower.tail=FALSE)
		else pWC <- NA
		if(!is.na(WT) && WT > 0) pWT <- pchisq(WT, dfWT,lower.tail=FALSE)
		else pWT <- NA
		if(!is.na(WAC) && WAC > 0) pWAC <- pchisq(WAC, dfWAC,lower.tail=FALSE)
		else pWAC <- NA
		if(!is.na(WCT) && WCT > 0) pWCT <- pchisq(WCT, dfWCT,lower.tail=FALSE)
		else pWCT <- NA
		if(!is.na(WAT) && WAT > 0) pWAT <- pchisq(WAT, dfWAT,lower.tail=FALSE)
		else pWAT <- NA
		if(!is.na(WACT) && WACT > 0) pWACT <- pchisq(WACT, dfWACT,lower.tail=FALSE)
		else pWACT <- NA

		W <- rbind(WA, WC, WT, WAC, WCT, WAT, WACT)
		pW <- rbind(pWA, pWC, pWT, pWAC, pWCT, pWAT, pWACT)
		dfW <- rbind(dfWA, dfWC, dfWT, dfWAC, dfWCT, dfWAT, dfWACT)

		WaldType <- data.frame(W, dfW, pW)
		rd.WaldType <- round(WaldType, Inf)
		Wdesc <- rbind(factor.name, time1.name, time2.name, paste(factor.name,":", time1.name, sep="") , 
		paste(time2.name, ":", time1.name, sep=""), paste(factor.name, ":", time2.name, sep=""), paste(factor.name, 
		":", time1.name, ":", time2.name, sep=""))
		colnames(rd.WaldType) <- c("Statistic", "df","p-value")
		rownames(rd.WaldType) <- Wdesc

		fn.tr <- function(mat) 
		{
			return(sum(diag(mat)))
		}

		

		RTE.B <- RTE[(A + C + T + C*A + C*T + T*A + 1) : (A + C + T + C*A + C*T + T*A + A*C*T)]

		### Box type statistics computed here

		BtA <- t(CA)%*%(ginv(CA%*%t(CA)))%*%CA
		BtC <- t(CC)%*%(ginv(CC%*%t(CC)))%*%CC
		BtT <- t(CT)%*%(ginv(CT%*%t(CT)))%*%CT
		BtAC <- t(CAC)%*%(ginv(CAC%*%t(CAC)))%*%CAC
		BtCT <- t(CCT)%*%(ginv(CCT%*%t(CCT)))%*%CCT
		BtAT <- t(CAT)%*%(ginv(CAT%*%t(CAT)))%*%CAT
		BtACT <- t(CACT)%*%(ginv(CACT%*%t(CACT)))%*%CACT

		TVA <- BtA%*%V
		BA <- (N/fn.tr(TVA)) * ((t(RTE.B)) %*% BtA %*% (RTE.B))
		BAf <- ((fn.tr(BtA%*%V))^2)/(fn.tr(BtA%*%V%*%BtA%*%V))
		dpr <- PA*diag(A)
		mat <- kronecker(diag(A), (1/(C*T))*t(matrix(1,(C*T),1)))
		va <- mat%*%V%*%t(mat)
		lambda <- solve(diag(n.a.vec) - diag(A))
		tem1 <- (fn.tr(dpr%*%va))^2
		tem2 <- fn.tr(dpr%*%dpr%*%va%*%va%*%lambda)
		BAf2 <- tem1/tem2

		if((!is.na(BA))&&(!is.na(BAf))&&(!is.na(BAf2))&&(BA > 0)&&(BAf > 0)) 
		{
   		 BAp <- pf(BA, BAf,Inf,lower.tail=FALSE)
   		 BApmod <- pf(BA, BAf,BAf2,lower.tail=FALSE)
		}
		else BAp <- NA

		TVC <- BtC%*%V
		BC <- (N/fn.tr(TVC)) * ((t(RTE.B)) %*% BtC %*% (RTE.B))
		BCf <- ((fn.tr(BtC%*%V))^2)/(fn.tr(BtC%*%V%*%BtC%*%V))
		if((!is.na(BC))&&(!is.na(BCf))&&(BC > 0)&&(BCf > 0)) BCp <- pf(BC, BCf, Inf,lower.tail=FALSE)
		else BCp <- NA

		TVT <- BtT%*%V
		BT <- (N/fn.tr(TVT)) * ((t(RTE.B)) %*% BtT %*% (RTE.B))
		BTf <- ((fn.tr(BtT%*%V))^2)/(fn.tr(BtT%*%V%*%BtT%*%V))
		if((!is.na(BT))&&(!is.na(BTf))&&(BT > 0)&&(BTf > 0)) BTp <- pf(BT, BTf, Inf,lower.tail=FALSE)
		else BTp <- NA

		TVAC <- BtAC%*%V
		BAC <- (N/fn.tr(TVAC)) * ((t(RTE.B)) %*% BtAC %*% (RTE.B))
		BACf <- ((fn.tr(BtAC%*%V))^2)/(fn.tr(BtAC%*%V%*%BtAC%*%V))
		if((!is.na(BAC))&&(!is.na(BACf))&&(BAC > 0)&&(BACf > 0)) BACp <- pf(BAC, BACf, Inf,lower.tail=FALSE)
		else BACp <- NA

		TVCT <- BtCT%*%V
		BCT <- (N/fn.tr(TVCT)) * ((t(RTE.B)) %*% BtCT %*% (RTE.B))
		BCTf <- ((fn.tr(BtCT%*%V))^2)/(fn.tr(BtCT%*%V%*%BtCT%*%V))
		if((!is.na(BCT))&&(!is.na(BCTf))&&(BCT > 0)&&(BCTf > 0)) BCTp <- pf(BCT, BCTf, Inf,lower.tail=FALSE)
		else BCTp <- NA

		TVAT <- BtAT%*%V
		BAT <- (N/fn.tr(TVAT)) * ((t(RTE.B)) %*% BtAT %*% (RTE.B))
		BATf <- ((fn.tr(BtAT%*%V))^2)/(fn.tr(BtAT%*%V%*%BtAT%*%V))
		if((!is.na(BAT))&&(!is.na(BATf))&&(BAT > 0)&&(BATf > 0)) BATp <- pf(BAT, BATf, Inf,lower.tail=FALSE)
		else BATp <- NA

		TVACT <- BtACT%*%V
		BACT <- (N/fn.tr(TVACT)) * ((t(RTE.B)) %*% BtACT %*% (RTE.B))
		BACTf <- ((fn.tr(BtACT%*%V))^2)/(fn.tr(BtACT%*%V%*%BtACT%*%V))
		if((!is.na(BACT))&&(!is.na(BACTf))&&(BACT > 0)&&(BACTf > 0)) BACTp <- pf(BACT, BACTf, Inf,lower.tail=FALSE)
		else BACTp <- NA


		B <- rbind(BA, BC, BT, BAC, BCT, BAT, BACT)
		pB <- rbind(BAp, BCp, BTp, BACp, BCTp, BATp, BACTp)
		dfB <- rbind(BAf, BCf, BTf, BACf, BCTf, BATf, BACTf)

		BoxType <- cbind(B, dfB, pB)
		rd.BoxType <- round(BoxType, Inf)
		Bdesc <- rbind(factor.name, time1.name, time2.name, paste(factor.name,":", time1.name, sep=""), 
		paste(time2.name, ":", time1.name, sep=""), paste(factor.name, ":", time2.name, sep=""), 
		paste(factor.name, ":", time1.name, ":", time2.name, sep=""))
		colnames(rd.BoxType) <- c("Statistic", "df", "p-value")
		rownames(rd.BoxType) <- Bdesc

		BoxTypeMod <- cbind(BA, BAf, BAf2, BApmod)
		rd.BoxTypeMod <- round(BoxTypeMod, Inf)
		colnames(rd.BoxTypeMod) <- c("Statistic", "df1", "df2", "p-value")
		rownames(rd.BoxTypeMod) <- factor.name
	}

	if(SING.COV) 
	{
		cat("\n Warning(s):\n")
		cat(" The covariance matrix is singular. \n")
	}

	if(!((C > 1) && (T > 1))) cat("Wald-type and Anova-type statistics cannot be computed, since either C or T is 1")
   if (show.covariance == FALSE) {
        V <- NULL
    }

	out.f1.ld.f2 <- list(RTE=rd.PRes1,Wald.test=rd.WaldType,ANOVA.test=rd.BoxType,ANOVA.test.mod.Box=rd.BoxTypeMod,covariance=V,model.name=model.name) 


   if (plot.RTE == TRUE) {
        id.rte <- A + T + C + A * T + C * T + A * C
        plot.rte <- rd.PRes1[, 3][(id.rte + 1):(id.rte + A * 
            C * T)]
        frank.group <- rep(0, 0)
        frank.time1 <- rep(0, 0)
        frank.time2 <- rep(0, 0)
        for (i in 1:A) {
            for (j in 1:C) {
                for (k in 1:T) {
                  frank.group <- c(frank.group, paste(levs[i]))
                  frank.time1 <- c(frank.time1, paste(origsort1[j]))
                  frank.time2 <- c(frank.time2, paste(origsort2[k]))
                }
            }
        }
        Frank <- data.frame(Group = frank.group, Time1 = frank.time1, 
            Time2 = frank.time2, RTE = plot.rte)
        plot.samples <- split(Frank, Frank$Group)
        lev.T1 <- levels(factor(plot.samples[[1]]$Time1))
        lev.T2 <- levels(factor(plot.samples[[1]]$Time2))
        lev.F <- levels(factor(Frank$Group))
        par(mfrow = c(1, A))
        for (hh in 1:A) {
            id.g <- which(names(plot.samples) == origsortg[hh])
            plot(1:T, plot.samples[[id.g]]$RTE[1:T], pch = 10, 
                type = "b", ylim = c(0, 1.1), xaxt = "n", xlab = "", 
                ylab = "", cex.lab = 1.5, xlim = c(0, T + 1), 
                lwd = 3)
            title(main = paste(group.name, origsortg[hh]), xlab = paste(time2.name))
            for (s in 1:C) {
                points(1:T, plot.samples[[id.g]]$RTE[plot.samples[[hh]]$Time1 == 
                  lev.T1[s]], col = s, type = "b", lwd = 3)
            }
            axis(1, at = 1:T, labels = origsort2)
            legend("top", col = c(1:C), paste(time1.name, lev.T1), 
                pch = c(rep(10, C)), lwd = c(rep(3, C)))
        }
    }

	




	return(out.f1.ld.f2)
}
