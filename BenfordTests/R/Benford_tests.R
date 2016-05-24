#
#    Do not delete!
#  File name		Benford_tests.R
#  Part of:		   BenfordTests (GNU R contributed package)
#  Author:			Dieter William Joenssen
#  Copyright:		Dieter William Joenssen
#  Email:			Dieter.Joenssen@TU-Ilmenau.de
#  Created:		   16 April 2013
#  Last Update: 	16 July 2015
#  Description:	R code for Package BenfordTests. Implemented functionions include following:
#                 Actual Tests:
#                 -chisq.benftest                  ~ Chi square test for Benford's law
#                 -ks.benftest                     ~ Kologormov-Smirnov test for Benford's law
#                 -mdist.benftest                  ~ Chebyshev distance based test for Benford's law
#                 -edist.benftest                  ~ Euclidean distance based test for Benford's law
#                 -usq.benftest                    ~ Freedman Watson U square test for Benford's law
#                 -meandigit.benftest              ~ mean digit test for Benford's law
#                 -jpsq.benftest                   ~ Pearson correlation test for Benford's law (removed moved ability to choose "spearman" and "kendall" for correlation)
#                 -signifd.analysis                ~ function to (graphically) analyze each digit individualy.
#                 -jointdigit.benftest             ~ function to test all digit frequencies jointly
#                 Supporting functions:
#                 -signifd                  	   ~ returns the specified number of first significant digits
#                 -signifd.seq                	   ~ sequence of all possible k first digits(i.e., k=1 -> 1:9)
#                 -qbenf                           ~ returns full cumulative probability function for Benford's distribution (first k digits)
#                 -pbenf                           ~ returns full probability function for Benford's distribution (first k digits)
#                 -rbenf                           ~ generates a random variable that satisfies Benford's law
#                 -simulateH0                      ~ calculates the H0-Distribution of all test statistics via simulation
## Tests for Benford's law
## Pearson's Chi squared test statistic
chisq.benftest<-function(x=NULL,digits=1,pvalmethod="asymptotic",pvalsims=10000)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("asymptotic", "simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}

#reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
  #get the amount of values that should be tested
   n<-length(first_digits)
   #the the observed frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   #calculate the relative frequencies
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   #calculate the chi square test statistic
   chi_square<-n*sum((rel_freq_of_digits-rel_freq_of_digits_H0)^2/rel_freq_of_digits_H0)

   #calc pval if using the asymptotic NULL-distribution
   if(pvalmethod==1)
   {
		pval<-1-pchisq(chi_square,df=length(signifd.seq(digits))-1)
   }
   if(pvalmethod==2)#calc pval if using the simulated NULL-distribution
   {
   #wrapper function for simulating the NULL distribution
   dist_chisquareH0<-simulateH0(teststatistic="chisq",n=n,digits=digits,pvalsims=pvalsims)
   
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated chi_square value
      pval<-1-sum(dist_chisquareH0<=chi_square)/length(dist_chisquareH0)
   }
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(chisq = chi_square), p.value = pval, method = "Chi-Square Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

##Kologormov Smirnov test (EDF type)
ks.benftest<-function(x=NULL,digits=1,pvalmethod="simulate",pvalsims=10000)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
     #get the amount of values that should be tested
   n<-length(first_digits)
    #the the observed (relative)frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   
   #calculate the deviations in the cumulative sums
   cum_sum_Ds<-cumsum(rel_freq_of_digits)-cumsum(rel_freq_of_digits_H0)
   #calculate the K-S-D-statistic
   K_S_D<-max(max(cum_sum_Ds),abs(min(cum_sum_Ds)))*sqrt(n)
   
   if(pvalmethod==1)#calc pval if using the simulated NULL-distribution
   {
   #wrapper function for simulating the NULL distribution
   dist_K_S_D_H0<-simulateH0(teststatistic="ks",n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated D value
      pval<-1-sum(dist_K_S_D_H0<=K_S_D)/length(dist_K_S_D_H0)
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(D = K_S_D), p.value = pval, method = "K-S Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

# Chebyshev Distance Test (crit values for one digit testing first by Morrow)
mdist.benftest<-function(x=NULL,digits=1,pvalmethod="simulate",pvalsims=10000)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
     #get the amount of values that should be tested
   n<-length(first_digits)
    #the the observed (relative)frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   
   #calculate the m_star statisitic
   #on a personal note, this isn't very different from the KSD.
   m_star<-sqrt(n)*max(abs(rel_freq_of_digits-rel_freq_of_digits_H0))
   
   if(pvalmethod==1)
   {
      #wrapper function for simulating the NULL distribution
      dist_m_star_H0<-simulateH0(teststatistic="mdist",n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated m_star value
      pval<-1-sum(dist_m_star_H0<=m_star)/length(dist_m_star_H0)
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(m_star = m_star), p.value = pval, method = "Chebyshev Distance Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

# Euclidean Distance Test(crit values for one digit testing first by Morrow)
edist.benftest<-function(x=NULL,digits=1,pvalmethod="simulate",pvalsims=10000)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
     #get the amount of values that should be tested
   n<-length(first_digits)
    #the the observed (relative)frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   
   #calculate the m_star statisitic
   #on a personal note, this isn't very different from the m_star statistic.
   d_star<-sqrt(n)*sqrt(sum((rel_freq_of_digits-rel_freq_of_digits_H0)^2))
   
   if(pvalmethod==1)
   {
   #wrapper function for simulating the NULL distribution
   dist_d_star_H0<-simulateH0(teststatistic="edist",n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated d_star value
      pval<-1-sum(dist_d_star_H0<=d_star)/length(dist_d_star_H0)
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(d_star = d_star), p.value = pval, method = "Euclidean Distance Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

# Freedman's modification of Watsons U^2 for the Benford distribution (originally 1 digit)
usq.benftest<-function(x=NULL,digits=1,pvalmethod="simulate",pvalsims=10000)
{
 #some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
     #get the amount of values that should be tested
   n<-length(first_digits)
    #the the observed (relative)frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   
   #calculate deviations betwen the cumulative sums
   cum_sum_Ds<-cumsum(rel_freq_of_digits-rel_freq_of_digits_H0)
   #calculate the U^2 test statistic
   U_square<-(n/length(rel_freq_of_digits))*(sum(cum_sum_Ds^2)-((sum(cum_sum_Ds)^2)/length(rel_freq_of_digits)))
   
   if(pvalmethod==1)
   {
      #wrapper function for simulating the NULL distribution
      dist_U_square_H0<-simulateH0(teststatistic="usq",n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated U_square value
      pval<-1-sum(dist_U_square_H0<=U_square)/length(dist_U_square_H0)
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(U_square = U_square), p.value = pval, method = "Freedman-Watson U-squared Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

#Normed mean deviation test for Benfords distribution first proposed as descriptive test statistic by Judge and Schechter
meandigit.benftest<-function(x=NULL,digits=1,pvalmethod="asymptotic",pvalsims=10000)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("asymptotic", "simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #get specified number of first digits and number of numbers
   first_digits<-signifd(x,digits)
   n<-length(first_digits)
   #get empirical mean digit value
   mu_emp<-mean(first_digits)
   #get expected mean digit value under NULL
   mu_bed<-sum(signifd.seq(digits)*pbenf(digits))
   #get variance of mean digit value under NULL
   var_bed<-sum(((signifd.seq(digits)-mu_bed)^2)*pbenf(digits))
   #normalize to get a_star
   a_star<-abs(mu_emp-mu_bed)/(max(signifd.seq(digits))-mu_bed)
   
   #calc pval if using the asymptotic NULL-distribution
   if(pvalmethod==1)
   {
      pval<-(1-pnorm(a_star,mean=0,sd=sqrt(var_bed/n)/(9-mu_bed)))*2
   }
   if(pvalmethod==2)
   {
         #wrapper function for simulating the NULL distribution
         dist_a_star_H0<-simulateH0(teststatistic="meandigit",n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated a_star value
      pval<-1-sum(dist_a_star_H0<=a_star)/length(dist_a_star_H0)
	  #if this were a two-sided test, the p_value would be adjusted as follows:
	  #if(pval>.5){pval<- (1- pval)*2}
	  #else{pval<- pval*2}
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(a_star = a_star), p.value = pval, method = "Judge-Schechter Normed Deviation Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}


# Shapiro-Francia type (correlation based) test for Benford's distribution first proposed by Joenssen (2013)
jpsq.benftest<-function(x=NULL,digits=1,pvalmethod="simulate",pvalsims=10000)
{
 #some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
     #get the amount of values that should be tested
   n<-length(first_digits)
    #the the observed (relative)frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   
   #calculate the Jstat statistic
   J_stat_squ<-cor(rel_freq_of_digits,rel_freq_of_digits_H0)
   #square the J_stat_statistic and adjust the sign
   J_stat_squ<-sign(J_stat_squ)*(J_stat_squ^2)
   if(pvalmethod==1)
   {
   #wrapper function for simulating the NULL distribution
      dist_J_stat_H0<- simulateH0(teststatistic="jpsq",n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated J_stat_square value
      pval<-sum(dist_J_stat_H0<=J_stat_squ)/length(dist_J_stat_H0)
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(J_stat_squ = J_stat_squ), p.value = pval, method = "JP-Square Correlation Statistic Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

# Hotelling type test for Benford's distribution first proposed by Joenssen (2015)
jointdigit.benftest<-function(x = NULL, digits = 1, eigenvalues="all", tol = 1e-15, pvalmethod = "asymptotic", pvalsims = 10000)
{
   #some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("asymptotic"))#, "simulate"
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   #Might need this in the future
   decompose=TRUE
   
   #reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
   #get the amount of values that should be tested
   n<-length(first_digits)
   #the the observed frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   #calculate the relative frequencies
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   
   #calculate covariance matrix under the NULL
   covariance_matirx<-outer(rel_freq_of_digits_H0,rel_freq_of_digits_H0,"*")*-1
   diag(covariance_matirx)<-rel_freq_of_digits_H0*(1-rel_freq_of_digits_H0)
   #ignore multiplication by n b/c only eigenvectors are desired
   #covariance_matirx<-covariance_matirx*n
   if(decompose)
   {
      eigenval_vect<-eigen(covariance_matirx,symmetric = TRUE)
      eigenval_vect_result<-eigenval_vect
      # identify which eigenvalues = 0
      eigen_to_keep<-abs(eigenval_vect$values)>tol
      #toss out the eigenvalues... = 0;
      eigenval_vect$values<-eigenval_vect$values[eigen_to_keep]
      eigenval_vect$vectors<-eigenval_vect$vectors[,eigen_to_keep]
      
      #determine which nonzero eigenvalues to keep
      if(length(eigenvalues)>0)
      {
         if(is.character(eigenvalues))
         {
            if(length(eigenvalues)==1)
            {
               eigenvalues <- pmatch(tolower(eigenvalues), c("all","kaiser"))
               if(eigenvalues == 1)
               {
                  eigen_to_keep<-1:length(eigenval_vect$values)
               }
               if(eigenvalues == 2)
               {
                  eigen_to_keep<-which(eigenval_vect$values>=mean(eigenval_vect$values))
               }
            }
            else
            {stop("Error: 'is.character(eigenvalues) && length(eigenvalues)!=1', use only one string!")}
         }
         else
         {
            if(is.numeric(eigenvalues)&all(eigenvalues>=0,na.rm = TRUE))
            {
               eigen_to_keep<-eigenvalues[!is.na(eigenvalues)]
               eigen_to_keep<-eigen_to_keep[eigen_to_keep<=length(eigenval_vect$values)]
               if(length(eigen_to_keep)<=0)
               {stop("Error: No eigenvalues remain.")}
            }
            else
            {stop("Error: non string value for eigenvalues must numeric vector of eigenvalue indexes! No negative indexing allowed.")}
         }
      }else{stop("Error: 'length(eigenvalues)<=0'!")}
      
      
      #reduce to selected eigenvalues;
      eigenval_vect$values<-eigenval_vect$values[eigen_to_keep]
      eigenval_vect$vectors<-eigenval_vect$vectors[,eigen_to_keep]
      
      principle_components<-rel_freq_of_digits%*%eigenval_vect$vectors
      true_components_means<-rel_freq_of_digits_H0%*%eigenval_vect$vectors
      
      if(length(eigenval_vect$values)==1)
      {
         hotelling_T<-(n/eigenval_vect$values)*((principle_components-true_components_means)^2)
      }else{
         hotelling_T<-n*(principle_components-true_components_means)%*%solve(diag(eigenval_vect$values))%*%t(principle_components-true_components_means)
      }
      deg_free<-length(principle_components)
   }else{
      hotelling_T<-n*(rel_freq_of_digits-rel_freq_of_digits_H0)%*%solve(covariance_matirx)%*%t(rel_freq_of_digits-rel_freq_of_digits_H0)
      deg_free<-length(rel_freq_of_digits)
   }
   
   #calc pval if using the asymptotic NULL-distribution
   if(pvalmethod==1)
   {
      #pval<-pchisq(q = hotelling_T,df = length(principle_components))
      pval<-1-pchisq(q = hotelling_T,df = deg_free)
   }
   if(pvalmethod==2)#calc pval if using the simulated NULL-distribution
   {
      #Reserved for when implemented
      #wrapper function for simulating the NULL distribution
      #dist_chisquareH0<-simulateH0(teststatistic="chisq",n=n,digits=digits,pvalsims=pvalsims)
      
      #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated chi_square value
      #pval<-1-sum(dist_chisquareH0<=chi_square)/length(dist_chisquareH0)
   }
   
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(Tsquare = hotelling_T), p.value = pval, method = "Joint Digits Test", 
                data.name = deparse(substitute(x)),eigenvalues_tested=eigen_to_keep,eigen_val_vect=eigenval_vect_result)
   class(RVAL) <- "htest"
   return(RVAL)
}



### Aditional functions provided
## returns first "digits" significant digits of numerical vector x
signifd<-function(x=NULL, digits=1)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x needs to be numeric.")}
   #calculate the first significant digits
   x<-abs(x)
   return(trunc((10^((floor(log10(x))*-1)+digits-1))*x))   
}

##returns the sequence of all possible leading digits for "digits" leading digits
#ie 1-> 1:9; 2-> 10:99; 3-> 100:999 etc.
signifd.seq<-function(digits=1)
{return(seq(from=10^(digits-1),to=(10^(digits))-1))}

# returns complete cumulative distribution function of Benford distribution for the given amount of significant digits
qbenf<-function(digits=1)
{
	return(cumsum(pbenf(digits)))
}

# returns complete probability distribution function of Benford distribution for the given amount of significant digits
pbenf<-function(digits=1)
{
   pbenf_for_seq<-function(leaddigit=10)
   {
      return(log10(1+(1/leaddigit)))
   }
   benf_table<-table(signifd.seq(digits))-1
   benf_table<-benf_table+sapply(signifd.seq(digits),FUN=pbenf_for_seq)
   
   return(benf_table)
}

#returns a n-long sample numbers satisfying Benford's law
rbenf<-function(n)
{
   return(10^(runif(n)))
}

#simulates the H0-Distribution of various tests offered in BenfordTests
simulateH0<-function(teststatistic="chisq",n=10,digits=1,pvalsims=10)
{
   teststatistic<-match.arg(arg = teststatistic, choices = c("chisq","edist","jpsq","ks","mdist","meandigit","usq"), several.ok = FALSE)
   if(teststatistic=="chisq")
   {
      H0_chi_square<-rep(0,pvalsims)
      H0_chi_square<- .C("compute_H0_chi_square", H0_chi_square = as.double(H0_chi_square), digits = as.integer(digits),
                         pbenf = as.double(pbenf(digits)),qbenf=as.double(qbenf(digits)),n = as.integer(n),
                         n_sim=as.integer(pvalsims))$H0_chi_square
      return(H0_chi_square)
   }
   if(teststatistic=="edist")
   {
      H0_dstar <- rep(0, pvalsims)
      H0_dstar <- .C("compute_H0_dstar", H0_dstar = as.double(H0_dstar), 
                     digits = as.integer(digits), pbenf = as.double(pbenf(digits)), 
                     qbenf = as.double(qbenf(digits)), n = as.integer(n), 
                     n_sim = as.integer(pvalsims))$H0_dstar
      return(H0_dstar)
   }
   if(teststatistic=="jpsq")
   {
      H0_J_stat <- rep(0, pvalsims)
      H0_J_stat <- .C("compute_H0_J_stat", H0_J_stat = as.double(H0_J_stat), 
                      digits = as.integer(digits), pbenf = as.double(pbenf(digits)), 
                      qbenf = as.double(qbenf(digits)), n = as.integer(n), 
                      n_sim = as.integer(pvalsims))$H0_J_stat
      return(H0_J_stat)
   }
   if(teststatistic=="ks")
   {
      H0_KSD <- rep(0, pvalsims)
      H0_KSD <- .C("compute_H0_KSD", H0_KSD = as.double(H0_KSD), 
                   digits = as.integer(digits), pbenf = as.double(pbenf(digits)), 
                   qbenf = as.double(qbenf(digits)), n = as.integer(n), 
                   n_sim = as.integer(pvalsims))$H0_KSD
      return(H0_KSD)
   }
   if(teststatistic=="mdist")
   {
      H0_mstar <- rep(0, pvalsims)
      H0_mstar <- .C("compute_H0_mstar", H0_mstar = as.double(H0_mstar), 
                     digits = as.integer(digits), pbenf = as.double(pbenf(digits)), 
                     qbenf = as.double(qbenf(digits)), n = as.integer(n), 
                     n_sim = as.integer(pvalsims))$H0_mstar
      return(H0_mstar)
   }
   if(teststatistic=="meandigit")
   {
      H0_astar <- rep(0, pvalsims)
      H0_astar <- .C("compute_H0_astar", H0_astar = as.double(H0_astar), 
                     digits = as.integer(digits), pbenf = as.double(pbenf(digits)), 
                     qbenf = as.double(qbenf(digits)), n = as.integer(n), 
                     n_sim = as.integer(pvalsims))$H0_astar
      return(H0_astar)
   }
   if(teststatistic=="usq")
   {
      H0_U_square <- rep(0, pvalsims)
      H0_U_square <- .C("compute_H0_U_square", H0_U_square = as.double(H0_U_square), 
                        digits = as.integer(digits), pbenf = as.double(pbenf(digits)), 
                        qbenf = as.double(qbenf(digits)), n = as.integer(n), 
                        n_sim = as.integer(pvalsims))$H0_U_square
      return(H0_U_square)
   }
   
   
}

#a method for graphically analyzing the first digits. Also gives pvalues for investigating each individual digit
signifd.analysis<-function(x=NULL,digits=1,graphical_analysis=TRUE,freq=FALSE,alphas=20,tick_col="red",ci_col="darkgreen",ci_lines=c(.05))
{
   if(length(alphas)==1)
   {
      if(alphas>1)
      {
         alphas=seq(from=0,to=.5,length.out=alphas+2)[-c(1,alphas+2)]
      }
   }
   n<-length(x)
   first_digits<-signifd(x, digits)
   pdf_benf<-pbenf(digits)
   
   
   freq_of_digits <- table(c(first_digits, signifd.seq(digits))) - 1
   E_vals<-pdf_benf*n
   Var_vals<-pdf_benf*n*(1-pdf_benf)
   Cov_vals<-outer(pdf_benf,pdf_benf)*-1*n
   diag(Cov_vals)<-Var_vals
   
   pval<-rep(0,length(pdf_benf))
   
   for(i in 1:length(pval))
   {
      pval[i]<-pnorm(q=freq_of_digits[i],mean=E_vals[i],sd=sqrt(Var_vals[i]))
      if(pval[i]>.5)
      {pval[i]<- (1- pval[i])*2}else{pval[i]<- pval[i]*2}
   }
   
   if(graphical_analysis)
   {
      mids<-seq(from=0,to=1,length.out=length(E_vals)+2)
      ci_line_length<-(mids[2]-mids[1])*(2/5)
      mids<-mids[-c(1,length(mids))]
      ## number formating function
      numformat <- function(val,trailing=4) { sub("^(-?)0.", "\\1.", sprintf(paste(sep="","%.",trailing,"f"), val)) }
      trailing<-0
      
      
      ci_cols<-colorRampPalette(colors=c("white",ci_col),interpolate="linear")(length(alphas)+1)[-1]
      ci_cols<-c(ci_cols,rev(ci_cols))
      
      alphas<-c(alphas/2,0.5,rev(1-(alphas/2)))
      cis<-sapply(alphas,FUN=qnorm,mean=E_vals,sd=sqrt(Var_vals))
      CIs=t(cis)
      colnames(CIs)<-signifd.seq(digits)
      rownames(CIs)<-alphas
      if(!freq)
      {
         cis<-cis/n
         freq_of_digits<-freq_of_digits/n
         
         #if(sum(cis>1)>0){warning("n is small, normal approximation may not be accurate!\nSome confidence intervals > 1.",call.=FALSE)}
         cis[cis>1]<-1
         trailing<-4
      }
      results<-list(summary=rbind(freq=freq_of_digits,pvals=pval),CIs=CIs)
      
      #if(sum(cis<0)>0){warning("n is small, normal approximation may not be accurate!\n Some confidence intervals <0.",call.=FALSE)}
      cis[cis<0]<-0
      
      lr_mid<-cbind(mids-ci_line_length,mids+ci_line_length)
      
      
      plot(x=0,y=0,xlim=c(0,1),ylim=c(0,max(cis)*1.3),type="n",axes=FALSE,xlab="summary",ylab="")
      for(i in 1:dim(cis)[1])
      {
         for(j in 1:(dim(cis)[2]-1))
         {
            polygon(x=lr_mid[i,c(1,1,2,2)],y=cis[i,c(j,j+1,j+1,j)],col=ci_cols[j],border=FALSE)
         }
      }
      dim_cis<-dim(cis)
      dim(cis)<-NULL
      
      
      posy<-seq(from=0,to=max(cis)*1.3,length.out=10)
      axis(side=2,at=round(posy-5*(10^-(digits+1)),digits),las=1)
      
      if(any(ci_lines!=FALSE))
      {
         if((!is.logical(ci_lines))&(all(ci_lines<1)&all(ci_lines>0)))
         {
            ci_lines<-c(ci_lines/2,0.5,rev(1-(ci_lines/2)))
            
            cis<-sapply(ci_lines,FUN=qnorm,mean=E_vals,sd=sqrt(Var_vals))
            if(!freq)
            {cis<-cis/n}
            CIs=t(cis)
            colnames(CIs)<-signifd.seq(digits)
            rownames(CIs)<-ci_lines
            results$CIs<-CIs
            dim_cis<-dim(cis)
            dim(cis)<-NULL
            if(!freq)
            {cis[cis>1]<-1}
            cis[cis<0]<-0
            j<-1
            for(i in 1:length(cis))
            {
               lines(lr_mid[j,],rep(cis[i],2))
               if(j==dim(lr_mid)[1])
               {j<-1}
               else
               {j<-j+1}
            }
         }
         else
         {
            j<-1
            for(i in 1:length(cis))
            {
               lines(lr_mid[j,],rep(cis[i],2))
               if(j==dim(lr_mid)[1])
               {j<-1}
               else
               {j<-j+1}
            }
         }
      }
      
      
      points(mids,freq_of_digits,col=tick_col,pch=3)
      if(digits==1)
      {
         mtext(c("digit: ",names(E_vals)),side=1,line=0,at=c(-1*ci_line_length,mids))
         if(freq){mtext(c("freq:  ",numformat(freq_of_digits,trailing)),side=1,line=1,at=c(-1*ci_line_length,mids))}
         else{mtext(c("rel freq:  ",numformat(freq_of_digits,trailing)),side=1,line=1,at=c(-1*ci_line_length,mids))}
         
         mtext(c("pval:  ",numformat(pval)),side=1,line=2,at=c(-1*ci_line_length,mids))
      }
      dim(cis)<-dim_cis
      abline(h=0)
   }
   
   if(!graphical_analysis)
   {
      if(any(ci_lines!=FALSE)&(!is.logical(ci_lines))&(all(ci_lines<1)&all(ci_lines>0)))
      {alphas<-c(ci_lines/2,0.5,rev(1-(ci_lines/2)))}
      else
      {alphas<-c(alphas/2,0.5,rev(1-(alphas/2)))}
      cis<-sapply(alphas,FUN=qnorm,mean=E_vals,sd=sqrt(Var_vals))
      if(!freq)
      {
         cis<-cis/n
         freq_of_digits<-freq_of_digits/n
         
         #if(sum(cis>1)>0){warning("n is small, normal approximation may not be accurate!\nSome confidence intervals > 1.",call.=FALSE)}
         #cis[cis>1]<-1
      }
      #if(sum(cis<0)>0){warning("n is small, normal approximation may not be accurate!\n Some confidence intervals <0.",call.=FALSE)}
      #cis[cis<0]<-0
      CIs=t(cis)
      colnames(CIs)<-signifd.seq(digits)
      rownames(CIs)<-alphas
      results<-list(summary=rbind(freq=freq_of_digits,pvals=pval),CIs=CIs)
      return(results)
   }
   return(results)
}