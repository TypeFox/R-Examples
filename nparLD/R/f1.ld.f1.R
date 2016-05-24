# R code for F1_LD_F1 macro
#
# Input:
#		y: a vector of variable of interest
#		group: a vector of group variable (factor level)
#		time : a vector of time variable
#		subject : a vector of independent subjects
#
# Optional Input:
#		w.pat: pattern matrix of order group level x time level
#		w.t : vector of order time level, pattern for interaction
#		w.g : vector of order group level, group pattern
#		time.name: name of the time vector. "Time" is set as default.
#		group.name: name of the time vector. "Group" is set as default.
#               description: description of the output. Default is set to TRUE (show description)
#		time.order: a vector of time levels specifying the order.
#		group.order: a vector of group levels specifying the order.
#
# Output:
#               list of relative treatment effects, test results, pattern result
#
f1.ld.f1 <- function(y, time, group, subject, w.pat=NULL, w.t=NULL, w.g=NULL, time.name="Time", group.name="Group", 
description=TRUE, time.order=NULL, group.order=NULL,plot.RTE=TRUE,show.covariance=FALSE, order.warning=TRUE)
{
#        For model description see Brunner et al. (2002)
#
#        Author: Mahbub Latif (mlatif@gwdg.de)
#                     Department of Medical Statistics, Goettingen, Germany
#
#         Version:  01-01
#         Date: February 18, 2003
#
#        Editied by: Kimihiro Noguchi
#         Version:  01-02
#         Date: August 18, 2009
#
#        Editied by: Kimihiro Noguchi
#         Version:  01-03
#         Date: December 24, 2009
#
#    Key Variables:
#                time: time factor
#                t: number of levels of time
#                a: number of levels of group
#                N: total number of observations
#                ind: indicator of whether there exists a missing observation (0=Yes,1=No)
#                N.na: total number of missing observations
#                subject: total number of subject
#                rscore: ranks of the variable of interest
#                rankmean: mean rank for each level of time
#                Nobs: total number of observations for each level of time
#                RTE: relative treatment effects

#    check whether the input variables are entered correctly

	   var<-y
	   if(is.null(var)||is.null(time)||is.null(group)||is.null(subject))
		stop("At least one of the input parameters (y, time, group, or subject) is not found.")

           sublen<-length(subject)
	   varlen<-length(var)
	   timlen<-length(time)
	   grolen<-length(group)

	   if((sublen!=varlen)||(sublen!=timlen)||(sublen!=grolen))
		stop("At least one of the input parameters (y, time, group, or subject) has a different length.")

#####################################################################################
# The following are the helper functions for the main function
# List of functions:
# rte: outputs the relative treatment effect
# case2x2: outputs statistics for 2 x 2 design
# wald.test: outputs Wald-type test statistics
# ANOVA.test: outputs ANOVA-type test statistics
# Simple.time.test: outputs test statistics for time effect
# pair.comp.test: outputs test statistics for paired comparison test statistics
# pattern.group: outputs test statistics for patterned alternatives for group effects
# df.p: calculates degrees of freedom for patterned alternatives
# one: changes matrix to a vector
# I: creates an identity matrix
# J: creates a unit matrix
# count.subj: counts the number of subjects in each level of group
# vi: calculates the variance equation (8.18)
# V: calculates the block diagonal covariance matrix
# mean.factor: calculates the mean of each factor
# df: calculates the degrees of freedom
# tr: calculates the trace of the matrix
######################################################################################

	# rte function for relative treatment effects
	rte <- function(group, time, indx, rscore)
	{
		a <- nlevels(group);
		t <- nlevels(time);
		tab <- t(matrix(mean.factor(rscore, group:time, indx), t, a))
       		rankmean.g <- as.vector(apply(tab, 1, mean));
        	rankmean.y <- as.vector(apply(tab, 2, mean));
		rankmean.gy <- mean.factor(rscore,group:time, indx)
		RankMean <- c(rankmean.g, rankmean.y, rankmean.gy);

		# no of observations per factors
		Nobs <- c(tapply(indx, group, sum), tapply(indx, time, sum), tapply(indx, group:time, sum))

		# rte
		RTE <- (1/sum(indx))*(RankMean - 0.5)

		# output

		out <- data.frame(RankMeans=RankMean, Nobs=Nobs, RTE=RTE)
        	levels(group) <- paste(group.name,glevel,sep="")
       		levels(time) <- paste(time.name,tlevel,sep="")
		row.names(out)<-c(levels(group), levels(time), levels(group:time))
		return(out)
	}

	## Case 2 x 2 as described in 8.1.2.
	case2x2 <- function(group, time, subj, rscore, ind)
	{
       		rscore <- rscore*ind
       		rscore.s <- split(rscore, group)
       		subj.s <- split(subj, group)
       		ind.s <- split(ind, group)
       		sigma2 <- rep(0,2)
       		Un <- 0
       		Un.const <- 0
       		v.den <- 0
       		UnT <- 0
       		UnT.c <- 0
       		vT <- 0
       		UnAT <- 0

       		for(i in 1:2)
       		{
          		junk <- t(sapply(split(as.vector(rscore.s[[i]]), as.vector(subj.s[[i]])), matrix))
          		junk.i <- t(sapply(split(as.vector(ind.s[[i]]), as.vector(subj.s[[i]])), matrix))
			junk<-as.matrix(junk)
			junk.i<-as.matrix(junk.i)
          		junk.m <- apply(junk,2,sum)/apply(junk.i,2,sum)
          		junk.s <- apply(junk,1,sum)

          		# Group effect
          		Un <- Un + sum(junk.m)*(-1)^(i+1)
          		sigma2[i] <- sum((junk.s - sum(junk.m))^2)/(nrow(junk)-1)
          		Un.const <- Un.const + sigma2[i]/nrow(junk)
          		v.den <- v.den+(sigma2[i]/nrow(junk))^2/(nrow(junk)-1)

          		# Time
          		junk.d <- apply(junk,1,diff)
          		tau2 <- sum((junk.d - diff(junk.m))^2)/(nrow(junk)-1)
          		UnT <- UnT - diff(junk.m)
          		UnT.c <- UnT.c + tau2/nrow(junk)
          		vT <- vT + (tau2/nrow(junk))^2/(nrow(junk)-1)

          		# Interaction
          		UnAT <- UnAT + diff(junk.m)*(-1)^i
       		}

       		# Group effect
       		Un <- Un/sqrt(Un.const)
       		v <- Un.const^2/v.den
		if(!is.na(Un)&&(v > 0))
       		{
			pGN <- (pnorm(abs(Un),lower.tail=FALSE))*2
       			pGT <- (pt(abs(Un), v,lower.tail=FALSE))*2
		}
		else
		{
			pGN <- NA
			pGT <- NA
		}
       		out <- data.frame(Statistics=Un, NN=pGN, DF=v, tt=pGT)

       		# Time effect
       		UnT <- UnT/sqrt(UnT.c)
       		vT <- UnT.c^2/vT
		if(!is.na(UnT)&&(vT > 0))
       		{
       			pTN <- (pnorm(abs(UnT),lower.tail=FALSE))*2
       			pTT <- (pt(abs(UnT),vT,lower.tail=FALSE))*2
		}
		else
		{
			pTN <- NA
			pTT <- NA
		}
       		out <- rbind(out, c(UnT,pTN,vT,pTT))

       		# Interaction
       		UnAT <- UnAT/sqrt(UnT.c)
		if(!is.na(UnAT)&&(vT > 0))
		{
       			pATN <- (pnorm(abs(UnAT),lower.tail=FALSE))*2
       			pATT <- (pt(abs(UnAT),vT,lower.tail=FALSE))*2
		}
		else
		{
			pATN <- NA
			pATT <- NA
		}
       		out<- rbind(out, c(UnAT,pATN,vT,pATT))
       		names(out) <- c("Statistic","p-value(N)","df","p-value(T)")
       		row.names(out) <- c(group.name, time.name, paste(group.name,":",time.name,sep=""))
       		return(list(case2x2=out))
	}

	# Wald test to test average group effect, average time effect, and global interaction effect
	wald.test <- function(group, time, subject, rscore, ind, ni)
	{
		n <- sum(ni);
		N <- sum(ind)
		a <- nlevels(group)
		t <- nlevels(time)
		V <- V(group, time, subject, rscore, ind, a, t, ni)$V
		R <- V(group, time, subject, rscore, ind, a, t, ni)$R

		# unconditional time mean
		tab <- t(matrix(mean.factor(rscore, group:time, ind), t, a))
        	t.mean <- as.vector(apply(tab, 2, mean));

		# Average group effect
		p <- (R - 0.5)/N; # RTE group x time
        	Pa <- I(a) - (1/a) * J(a)
		Pt <- I(t) - (1/t) * J(t)
		Pat <- kronecker(Pa, Pt)

		# second last equation of page 134
		cpg <- sqrt(n) * Pa %*% kronecker(I(a), (1/t)*t(one(t))) %*% p;

		# last equation of page 134
		Sigma <- kronecker(Pa,  (1/t)*t(one(t))) %*% V %*%  kronecker(Pa, (1/t)*one(t));

		# equation (8.10)
        	cvc <- Pa %*% Sigma %*% Pa
		Q.a <- t(cpg) %*% ginv(cvc) %*% cpg;

	        df.a <- tr(cvc%*%ginv(cvc))
		if(!is.na(Q.a) && (Q.a > 0)) pval.a <- round(pchisq(Q.a, df.a,lower.tail=FALSE),Inf)
		else pval.a <- NA;
		A <- c(W=Q.a, df=df.a, pval=pval.a);

		# Average time effect, eqn (8.9)
		S <- kronecker((1/a)*t(one(a)), I(t)) %*% V %*% kronecker((1/a)*one(a), I(t))

		cpt <- Pt %*% t.mean;

		# equation (8.20)
        	cvc <- Pt %*% S %*% Pt
		Q.t <- (n/N^2)*t(cpt) %*% ginv(cvc) %*% cpt;
        	df.t <- tr(cvc%*%ginv(cvc))
		if(!is.na(Q.t) && (Q.t > 0)) pval.t <- round(pchisq(Q.t, df.t,lower.tail=FALSE),Inf)
		else pval.t <- NA;
		T <- c(Q.t, df.t, pval.t);

		# Global interaction effect
        	Cat <- kronecker(Pa, Pt); # book notation of page 141

		# equation (8.26)
        	cvc <- Cat %*% V %*% t(Cat)
        	Q.at <- (n/N^2)*t(Cat %*% R) %*% ginv(cvc) %*% Cat %*% R;
        	df.at <- tr(cvc%*%ginv(cvc))
		if(!is.na(Q.at) && (Q.at > 0)) pval.at <- round(pchisq(Q.at, df.at,lower.tail=FALSE), Inf)
		else pval.at <- NA;
		AT <- c(Q.at, df.at, pval.at);

		# results
		out.w <- rbind(A, T, AT);
	        colnames(out.w) <- c("Statistic", "df", "p-value")
                rownames(out.w) <- c(group.name,time.name,paste(group.name,":",time.name,sep=""))
		out <- list(Wald.test=out.w);
	}

	# To test average group effect, average time effect, and global interaction effect
	anova.test <- function(group, time, subject, rscore, ind, a, ni)
	{
		group <- as.factor(group)
		time <- as.factor(time)
		t <- nlevels(time)

		t.mean <- apply(t(matrix(mean.factor(rscore, group:time, ind), t, a)),2,mean)

		n <- sum(ni);
		N <- sum(ind)

		V <- V(group, time, subject, rscore, ind, a, t, ni)$V
		R <- V(group, time, subject, rscore, ind, a, t, ni)$R

	     	 # Average group effect
 	     	p <- (R - 0.5)/N; # RTE corresponding group x time
 	   	Pa <- I(a) - (1/a) * J(a); # centering matrix
 	        Pt <- I(t) - (1/t) * J(t); # centering matrix
 	        Pat <- kronecker(Pa, Pt);

  	        # second last equation of page 134
 	        cpg <- sqrt(n) * Pa %*% kronecker(I(a), (1/t)*t(one(t))) %*% p;
		Sigma <- kronecker(Pa,  (1/t)*t(one(t))) %*% V %*%  kronecker(Pa, (1/t)*one(t));

		# average group effect
		# equation (8.11)
		F.a <- t(cpg) %*% cpg/sum(diag(Pa %*% Sigma));
		# equation (5.7)
		df1.a <- sum(diag(Pa %*% Sigma))^2/sum(diag((Pa %*% Sigma) %*% (Pa %*% Sigma)));
		if((!is.na(F.a))&&(!is.na(df1.a))&&(F.a > 0)&&(df1.a > 0)) pval.a<-pchisq(F.a*df1.a,df1.a,lower.tail=FALSE)
		else pval.a<-NA;
		A <- round(c(B=F.a, df=df1.a, pval=pval.a), Inf);

		# average time effect
  	      	# equation (8.9)
 	       S <- kronecker((1/a)*t(one(a)), I(t)) %*% V %*% kronecker((1/a)*one(a), I(t))

               cpt <- Pt %*% t.mean;

		# equation (8.21)
		F.t <- (n/N^2) * (t(cpt) %*% cpt)/sum(diag(Pt %*% S));
		df1.t <- sum(diag(Pt %*% S))^2/sum(diag(Pt %*% S %*% Pt %*% S));
		if((!is.na(F.t))&&(!is.na(df1.t))&&(F.t > 0)&&(df1.t > 0)) pval.t<-pchisq(F.t*df1.t,df1.t,lower.tail=FALSE)
		else pval.t<-NA;
		T <- round(c(F.t, df1.t, pval.t),Inf);

		# Global interaction effect
		F.at <- n * t(p) %*% Pat %*% p/sum(diag(Pat %*% V));
		df1.at <- sum(diag(Pat %*% V))^2/sum(diag(Pat %*% V %*% Pat %*% V));
		if((!is.na(F.at))&&(!is.na(df1.at))&&(F.at > 0)&&(df1.at > 0)) pval.at<-pchisq(F.at*df1.at,df1.at,lower.tail=FALSE)
		else pval.at<-NA;
		AT <- round(c(F.at, df1.at, pval.at), Inf);
		out.box <- rbind(A, T, AT);
	        colnames(out.box) <- c("Statistic", "df", "p-value")
                rownames(out.box) <- c(group.name,time.name,paste(group.name,":",time.name,sep=""))
		# modified Box-approximation
		df1 <- df(V, a, t, ni, ind)$df1
		df2 <- df(V, a, t, ni, ind)$df2
		if((!is.na(F.a)) && (!is.na(df1)) && (!is.na(df2)) && (F.a > 0) && (df1 > 0) && (df2 > 0)) pval.mb <- pf(F.a, df1, df2,lower.tail=FALSE)
		else pval.mb <- NA
		A <- rbind(round(c(B=F.a, df1=df1, df2=df2, pval=pval.mb), Inf));
	        colnames(A) <- c("Statistic", "df1", "df2","p-value")
                rownames(A) <- c(group.name)
		out <- list(ANOVA.test=out.box, ANOVA.test.mod.Box=A);

		return(out);
	}

	# Simple time test to test time effect
	simple.time.test <- function(name.group, a , t, ni, N, pat.dat, V1, R1)
	{
		# simple time effect
		Pt <- I(t) - (1/t) * J(t);
		wald <- as.data.frame(matrix(0, a, 3));
		rownames(wald)<-name.group
		anova <- as.data.frame(matrix(0, a, 3));
		rownames(anova)<-name.group
		normal <- as.data.frame(matrix(0, a, 4));
		rownames(normal)<-name.group
		names(wald) <- c("Statistic", "df", "p-value");
		names(anova) <- c("Statistic", "df", "p-value");
		names(normal) <- c("Statistic", "p-value(N)", "df", "p-value(T)");
		n <- sum(ni); k <-1;

		for(i in 1:a)
		{
			V <- V1[k:(i*t), k:(i*t)];
			R <- R1[k:(i*t),1];
			k <- i*t + 1;
			cp <- Pt %*% R;
			Q.s <- round((n/N^2) * t(cp) %*% ginv(Pt %*% V %*% Pt) %*% cp, Inf);
			F.s <- round((n/N^2) * t(cp) %*% cp/sum(diag(Pt %*% V)), Inf);
                	df <- tr((Pt %*% V %*% Pt)%*%ginv(Pt %*% V %*% Pt))
			df1 <- round(sum(diag(Pt %*% V))^2/sum(diag(Pt %*% V %*% Pt %*% V)), Inf);
			if((!is.na(Q.s)) && (!is.na(df)) && (Q.s > 0) && (df > 0)) pval <- round(pchisq(Q.s, df,lower.tail=FALSE), Inf)
			else pval <- NA;
			if((!is.na(F.s)) && (!is.na(df1)) && (F.s > 0) && (df1 > 0)) pval1<-round(pchisq(F.s*df1,df1,lower.tail=FALSE), Inf)
			else pval1 <- NA
			out <- c(Q.s, df, pval);
			out1 <- c(F.s, df1, pval1);
			wald[i,] <- out;
			anova[i,] <- out1;

			# pattern effect
			if(!is.null(pat.dat))
			{
				pi <- (R - 0.05)/N;
	     	     	        s2 <- (ni[i]/n) * as.numeric(pat.dat[i,]) %*% Pt %*% V %*% Pt %*% as.numeric(pat.dat[i,]);
        		        L <- round(sqrt(ni[i]/s2) * as.numeric(pat.dat[i,]) %*% Pt %*% pi, Inf);
				p.nor <- round((pnorm(L,lower.tail=FALSE)), Inf);
				df1 <- ni[i] -1;
				pval <- round((pt(L, df1,lower.tail=FALSE)), Inf);
				normal[i,] <- c(L, p.nor, df1, pval);
			}
		}

		sim.time.effect <- list(Wald.test.time=wald, ANOVA.test.time=anova);
		if(!is.null(pat.dat)) pat.time.effect <- list(pattern.time=normal)
		else pat.time.effect <- list(pattern.time=NULL)

		return(c(sim.time.effect, pat.time.effect));
	}

	# pairwise comparison test statistics
	pair.comp.test <- function(data, ni, w, lev.grp)
	{
		a <- nlevels(factor(data[,1]));
		t <- nlevels(factor(data[,2]));
		n <- sum(ni)
		N <- sum(data[,5])

		Pt <- I(t) - (1/t)*J(t)
		V11 <- V(data[,1], data[,2], data[,3], data[,4], data[,5], a, t, ni)$V;

		# arranging output
		out <- as.data.frame(matrix(0, 3*choose(a,2), 3));
		out.pat <- as.data.frame(matrix(0, choose(a,2), 4));
		names(out) <- c("Statistic", "df", "p-value");
		names(out.pat) <- c("Statistic", "p-value(N)", "df", "p-value(T)")
		Test <- rep(c(group.name,time.name,paste(group.name,":",time.name,sep="")), choose(a,2));
		Pairs <- rep(0, 3*choose(a ,2));


		k <- 1;
		ll <-1
		for(i in 1:(a-1))
		{
			for(j in (i+1):a)
			{
				data.p <- data[data[,1]==i | data[,1]==j,]
				ni.p <- matrix(c(ni[i],ni[j]), 2,1);

				nn <- sum(ni.p)
				NN <- sum(data.p[,5]);

				gr <- data.p[,1]
				tm <- data.p[,2]
				subj <- data.p[,3]
				rs <- data.p[,4]
				ind <- data.p[,5]

				V <- V(as.numeric(gr), tm, subj, rs, ind, 2, t, ni.p)$V
				R <- V(as.numeric(gr), tm, subj, rs, ind, 2, t, ni.p)$R

				out[k:(k+2),] <- anova.test(gr, tm, subj, rs, ind, 2, ni.p)$ANOVA.test
				Pairs[k:(k+2)] <- paste(group.name,lev.grp[i], ":",group.name,lev.grp[j],sep="")
				k <- k + 3;

				# pattern interactions
				if(!is.null(w))
				{
					w <- matrix(w, t,1)
					sign <- t(w)%*%Pt%*%(V[1:t,1:t]+V[(t+1):(2*t),(t+1):(2*t)])%*%Pt%*%w
					out.pat[ll,1] <- sqrt(nn/sign)*t(w-mean(w))%*%(R[1:t,] - R[(t+1):(2*t),])/NN
					out.pat[ll,2] <- pnorm(out.pat[ll,1],lower.tail=FALSE)
					posA <- matrix(0, a, 1);
					posA[i,]<- 1; posA[j,] <- 1;
					CC <- t(w)%*%Pt
					M <- kronecker(diag(c(posA)),CC)
					S <- M%*%V11%*%t(M)
					lambda <- solve(diag(c(ni))-I(a))
					out.pat[ll,3] <- tr(S)^2/tr(S*S*lambda)
					out.pat[ll,4] <- pt(out.pat[ll,1],out.pat[ll,3],lower.tail=FALSE)
					row.names(out.pat)[ll] <- paste(group.name,lev.grp[i],":",group.name,lev.grp[j],sep="");
					ll <- ll + 1
				}
			}
		}

		out <- cbind(Pairs, Test, out);
		if(!is.null(w)) pair.comp<- list(pair.comparison=out, pattern.pair.comparison=round(out.pat,Inf))
		else pair.comp <- list(pair.comparison=out, pattern.pair.comparison=NULL)
		return(pair.comp);
	}

	# patterned alternatives for group effects
	pattern.group <- function(group, time, subject, rscore, ind, a, t, ni, g.mean, w.g)
	{
		n <- sum(ni)
		N <- sum(ind)
		Pa <- I(a) - (1/a)*J(a)
		V <- V(group, time, subject, rscore, ind, a, t, ni)$V
		S <- kronecker(Pa, (1/t)*t(one(t)))%*%V%*%kronecker(Pa, (1/t)*one(t))
		w.g <- matrix(w.g, a, 1)
		lambda <- solve(diag(c(ni))-I(a))
		c <- Pa%*%w.g

		g <- (g.mean-.5)/N
		Kn <- sqrt(n)*t(w.g)%*%Pa%*%g
		sign <- t(c)%*%S%*%c

		df<- df.p(V, a, t, ni, ind, w.g)

		Ln <- Kn/sqrt(sign)
		pval.t <- pt(Ln, df,lower.tail=FALSE)
		pval.n <- pnorm(Ln,lower.tail=FALSE)
		out<-rbind(round(c(Ln=Ln, pval.N=pval.n, df=df, pval.t=pval.t),Inf))
	        colnames(out) <- c("Statistic", "p-value(N)", "df", "p-value(T)")
                rownames(out) <- c(group.name)
		return(list(pattern.group=out))
	}

	# degrees of freedom calculatin for patterned alternative
	df.p <- function(V, a, t, ni, ind, w)
	{
		Pa <- I(a) - (1/a)*J(a)
		lambda <- solve(diag(c(ni)) - I(a))
		c <- (diag(c(t(w)%*%Pa)))^2
		m <- kronecker(I(a), (1/t)*t(one(t)))
		s <- m%*%V%*%t(m)
		x <- c*s
		df <- (tr(x))^2/tr((x%*%x)*lambda)
		return(df)
	}

	# one vector
	one <- function(d) return(matrix(1, d, 1));

	# Identity matrix
	I <- function(d)
	{
		junk <- rep(1, length=d);
		junk <- diag(junk);
		return(junk);
	}

	# Unit matrix
	J <- function(d1, d2=d1) return(matrix(1,d1,d2));

	# count the number of subjects in each level of group
	count.subj <- function(group, subject)
	{
		group <- as.factor(group)
		table <- table(group, subject);
		n <- matrix(0,nlevels(group), 1);
		for(i in 1:nlevels(group)) n[i] <- length(colnames(table)[table[i,] > 0]);
		return(n);
	}

	# variance Vi equation (8.18)
	Vi <- function(data)
	{
		time <- factor(data$time)
		subj <- factor(data$subj)
		srank <- data$rscore
		indx <- data$indx

		t <- nlevels(time)
		n0 <- nlevels(subj)

		# matrix of order subject x time, R_ik - R_i
		srank.m <- matrix(tapply(srank, subj:time, sum), t, n0)
		ind.m <- matrix(tapply(indx, subj:time, sum), t, n0)

		junk <- srank.m*ind.m;
		t.junk <- t(ind.m);

		R <- apply(junk,1,sum)/apply(ind.m, 1,sum);
		junk1 <- t(junk -R)

		V <- matrix(0, t, t);

		for(i in 1:t)
		{
			for(j in 1:t)
			{
				ls <- sum(t.junk[,i]);
				if(i==j) V[i,i] <- sum(t.junk[,i]*junk1[,i]^2)/(ls*(ls-1))
				else
				{
					lss <- sum(t.junk[,i]*t.junk[,j]);
					kss <- (ls-1)*(sum(t.junk[,j])-1) + lss - 1;
					V[i,j] <- sum(t.junk[,i]*t.junk[,j]*junk1[,i]*junk1[,j])/kss;
				}
			}
		}

		out <- list(R=R,V=V);
		return(out);
	}

	# block diagonal covariance matrix
	V <- function(group, time, subj, rscore, indx, a, t, ni)
	{
		group <- as.factor(group)
		subj <- as.factor(subj)
		time <- as.factor(time)
		data <- data.frame(group, time, subj, rscore, indx)

		n <- sum(ni);
		N <- sum(indx);

		# split rscore by group
		split.d <- split(data, data$group)

	        # calculating covariance matrix for average group, time and global interaction effect
        	V <- matrix(0, a*t, a*t);
        	R <- matrix(0, a*t, 1);
        	i <- 1; k <- 1;
        	while(i <= a)
        	{
                	V[k:(i*t), k:(i*t)] <-  Vi(split.d[[i]])$V;
                	R[k:(i*t), 1] <- Vi(split.d[[i]])$R;
                	k <- i*t + 1;
                	i <- i+1;
        	}
		V <- V*n/N^2;
		out <- list(V=V, R=R);
	}

	# mean of each factor
	mean.factor <- function(y, fact, indx)
	{
		tab <- tapply(y, fact, sum)
		ind <- tapply(indx, fact, sum)
		return(c(tab/ind))
	}

	# degrees of freedom calculation
	df <- function(V, a, t, ni, ind)
	{
		Pa <- I(a) - (1/a)*J(a)
		Pt <- I(t) - (1/t)*J(t)
		c <- kronecker(Pa, (1/t)*t(one(t)))
		tt <- t(c)%*%ginv(c%*%t(c))%*%c
		tem1 <- tt%*%V;
		df1 <- (tr(tem1))^2/tr(tem1%*%tem1)
		dpr <- Pa*I(a)
		mat <- kronecker(I(a), (1/t)*t(one(t)))
		va <- mat%*%V%*%t(mat)
		lambda <- solve(diag(c(ni)) - I(a))
		tem1 <- (tr(dpr%*%va))^2
		tem2 <- tr(dpr%*%dpr%*%va%*%va%*%lambda)
		df2 <- tem1/tem2
		return(list(df1=df1, df2=df2))
	}

	# trace calculation
	tr <- function(x) return(sum(diag(x)))

# end of helper functions
######################################################################################

# main function

	library(MASS)
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

	group<-factor(group)
	time<-factor(time)
	subject<-factor(subject)

        score <- var
	N.na <- sum(is.na(score))
	ind <- 1 - is.na(score)
	N <- sum(ind)
	rscore <- rank(score)*ind

	data <- cbind(group, time, subject, rscore, ind) # changed on 16 August, 2004
	ni <- count.subj(group, subject)
	n <- sum(ni); # number of subjects in the experiment
	   model.name<-"F1 LD F1 Model"

        if(description==TRUE)
        {
           	cat(" Total number of observations: ",sum(ind),"\n")
           	cat(" Total number of subjects:  " , n,"\n")
          	cat(" Total number of missing observations: ",N.na,"\n")
           	cat("\n Class level information ")
           	cat("\n ----------------------- \n")
	   	cat(" Levels of", time.name, "(sub-plot factor time) : ", t,"\n")
	   	cat(" Levels of", group.name, "(whole-plot factor group) : ", a,"\n")
           	cat("\n Abbreviations ")
           	cat("\n ----------------------- \n")
           	cat(" RankMeans = Rank means\n")
           	cat(" Nobs = Number of observations\n")
           	cat(" RTE = Relative treatment effect\n")
           	cat(" case2x2 = tests for 2-by-2 design\n")
           	cat(" Wald.test = Wald-type test statistic\n")
           	cat(" ANOVA.test = ANOVA-type test statistic with Box approximation\n")
           	cat(" ANOVA.test.mod.Box = modified ANOVA-type test statistic with Box approximation\n")
           	cat(" Wald.test.time = Wald-type test statistic for simple time effect\n")
           	cat(" ANOVA.test.time = ANOVA-type test statistic for simple time effect\n")
           	cat(" N = Standard Normal Distribution N(0,1)\n")
           	cat(" T = Student's T distribution with respective degrees of freedom\n")
           	if(!is.null(w.pat))
           	{
            		pattern.string<-c(w.pat)
           	}
           	else
           	{
           		pattern.string<-"no pattern specified"
           	}
           	if(!is.null(w.t))
           	{
            		pattern.string.t<-w.t
           	}
           	else
           	{
            		pattern.string.t<-"no pattern specified"
           	}
           	if(!is.null(w.g))
           	{
            		pattern.string.g<-w.g
           	}
           	else
           	{
            		pattern.string.g<-"no pattern specified"
           	}

           	cat(" pattern.time (time effects) = Test against patterned alternatives in time using normal distribution (",pattern.string,")","\n")
          	cat(" pair.comparison = Tests for pairwise comparisions (without specifying a pattern)","\n")
          	cat(" pattern.pair.comparison = Test for pairwise comparisons with patterned alternatives in time (",pattern.string.t,")","\n")
           	cat(" pattern.group (group effects) = Test against patterned alternatives in group (",pattern.string.g,")","\n")
           	cat(" covariance = Covariance matrix","\n")
           	cat(" Note: The description output above will disappear by setting description=FALSE in the input. See the help file for details.","\n\n")
      }
	
	if(order.warning==TRUE)
	{
           	cat(" F1 LD F1 Model ")
           	cat("\n ----------------------- \n")
           	cat(" Check that the order of the time and group levels are correct.\n")
           	cat(" Time level:  " , paste(tlevel),"\n")
           	cat(" Group level:  " , paste(glevel),"\n")
           	cat(" If the order is not correct, specify the correct order in time.order or group.order.\n\n")
	}

	# unconditional group and time means
	tab <- t(matrix(mean.factor(rscore, group:time, ind), t, a))
        g.mean <- as.vector(apply(tab, 1, mean));
        t.mean <- as.vector(apply(tab, 2, mean));

	# covariance matrix
	V2 <- V(group, time, subject, rscore, ind, a, t, ni)$V;
	R <- V(group, time, subject, rscore, ind, a, t, ni)$R;

	SING.COV <- FALSE
	if(qr(V2)$rank < (t*a)) SING.COV <- TRUE
	if(SING.COV)
	{
		cat("\n Warning(s):\n")
		cat(" The covariance matrix is singular. \n")
	}

        sdat <- NULL
	rte <- list(RTE=rte(group, time, ind, rscore))
	rte.plot<-data.frame(rte)
	namen.plot<-rownames(rte.plot)[(a+1):(a+t)]
	namen.plot.g<-rownames(rte.plot)[1:a]
	#rte <- data.frame(rte(group, time, ind, rscore))
	### case2x2 is available only when there is no missing observation in the 2-by-2 design.
	### otherwise, it returns NULL.
        if(a==2 && t==2 && N.na==0)
	{
		out2 <- case2x2(group, time, subject, rscore, ind)
	}
        else
	{
		out2 <- list(case2x2=NULL)
	}
	   wald.test.t <- wald.test(group, time, subject, rscore, ind, ni);
	   anova.test.t <- anova.test(group, time, subject, rscore, ind, a, ni);
	   out2 <- c(out2, wald.test.t, anova.test.t)
	   simple.time.test.t <- simple.time.test(glevel, a, t, ni, N, w.pat, V2, R);
           pair.comp.t <- pair.comp.test(data, ni, w.t, glevel)
           out2 <- c(out2, simple.time.test.t,pair.comp.t)

        if(!is.null(w.g)) pattern.g <- pattern.group(group, time, subject, rscore, ind, a, t, ni, g.mean, w.g)
        else pattern.g <- NULL
            if (show.covariance == FALSE) {
        V2 <- NULL}
        out <- c(sdat, rte, out2, pattern.g, list(covariance=V2), model.name=model.name)
     if (plot.RTE == TRUE) {
        ptg <- expression(p[i][j])
        rte.plot <- data.frame(rte)[, 3][(a + t + 1):(a + t + 
            a * t)]
        kk.time <- 1:t
        plot(kk.time, rte.plot[kk.time], pch = 10, type = "b", 
            lwd = 3, ylim = c(0, 1.1), xaxt = "n", xlab = time.name, 
            ylab = ptg, cex.lab = 1.5)
        axis(1, at = kk.time, labels = namen.plot[kk.time])
        group.id.help <- sort(c(rep(1:a, t)))
        for (ss in 2:a) {
            points(kk.time, rte.plot[group.id.help == ss], pch = 10, 
                type = "b", lwd = 3, ylim = c(0, 1.1), xaxt = "n", 
                xlab = "", ylab = ptg, cex.lab = 1.5, col = ss)
        }
        legend("top", col = c(1:a), namen.plot.g, pch = c(rep(10, 
            a)), lwd = c(rep(3, a)))
        title(main = paste("Relative Effects"))
    }



        return(out)
}

