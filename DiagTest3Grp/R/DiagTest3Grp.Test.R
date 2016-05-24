##########################################################################################
#Under normality assumption
#Omnibus Chi-square test compare >=3 markers to examine if all have statistically equivalent summary measures
#When exactly 2 markers are compared, z-test on equivalence in summary measures between two markers
#
#########################################################################################

DiagTest3Grp.Test <- function(dat,paired=FALSE,type=c("VUS","Youden"),p=0,q=0,mu=0,conf.level=0.95,alternative=c("two.sided","less","greater"))
  {
    
    ###Inputs:
     #(1) dat: For unpaired data (see paired argument), dat is a list of K components, each corresponding to a marker. Each component of the list should be a data frame of n observation by 2 columns: group must be at 1st column with levels of D-, D0 and D+ corresponding to increasing disease severity and numeric marker measurements at 2nd column ; for paired data, dat should be a data frame with observations at row and group/markers at column (group must be the 1st column with levels D-, D0. D+, followed by markers at 2nd to the last column
    
     #(2)paired: whether all markers are measured on the same set (paired=TRUE) or not (paired=FALSE)
     #(3)type: which summary measure (VUS or Youden) will be used to indicate the discriminative ability of the three-group diagnostic tests. Allows unique argument matching,i.e., can use "V" , "Y" as abbreviation 
     #(4)p, q: if testing on "VUS" and partial VUS satisfying specificity >=p and sensitivity >=q is desired, default=0 for full VUS.
     #(5):conf.level: default 95% confidence interval on point estimate (abs diff between two markers's summary measures)
     ##(6)alternative:fixed as  "two.sided" for testing >=3 markers.For testing two markers, can be "two.sided/less/greater" (marker1 at first column-marker2 at second column)
    
    ###outputs:
    ###(1)statistic: for independent samples, statistic is the calculated Chi-square-statistic with df=K-1;for paired samples, statistic is the calculated z-statistic following N(0,1)
    ###(2)p.value: calculated p-value
    ###(3)estimate: the markers' summary measures
    ###(4)conf.int: the 95% CI on the absoluate difference between the two markers' summary measures, only returned for two marker (independent samples or paired samples) , for >=3 markers, NA
    ###(5)method: "VUS" or "Youden", which summary measure was tested
    ###(6)
    BetweenMarker.cov <- function(marker1, marker2,n.minus,n0,n.plus,rho.minus,rho0,rho.plus,type="VUS")          
      {
           ###Calculate the covariance between the two markers' summary measures: cov(theta_1 and theta_2)
        
           ##input: (1)marker1, marker2: each a data frame (for marker 1, 2 respectively) containing:  a, b, c, d,(mostly a column of res0$patialDeriv outputted above) and V.a, V.b, V.c, V.d
           ##        (2)n.minus, n0,n.plus: the number of obs in each group that both markers are available
            
           #calculate the i<=4, j>4 between marker terms for eq (8) on P686 of Xiong et al on Biometrical,
           ###output: covariance
            
        mu.s.cov <- function(n,rho,s1,s2)
          {
            
            ##the covariance between the the two normal mean parameters in the same Diagnosis group for marker 1&2, e.g. in the D- group, cov(mu_{-1},mu_{-2})=rho_{-12}*s_{-1}*s{-2}/n_{-}
            mu.cov <- rho*s1*s2/n
                
            ##the covariance between the mu parameter in the same Diagnosis group for marker 1&2,e.g. in the D- group, cov(s_{-1},s_{-2})=s_{-1}*s{-2}*rho_{-12}^2*(1-rho_{-12}^2)^3 / [2*n*(3*rho_{-12}^4-rho_{-12}^6-3*rho_{-12}^2+1)]

            s.cov <- s1*s2*rho^2*(1-rho^2)^3 /(2*n*(3*rho^4-rho^6-3*rho^2+1))
            return(c(mu.cov=mu.cov,s.cov=s.cov))
          }
              
        if(type=="VUS")
          {
            ##the same param but two markers (index 1, 2 for marker 1, 2)
            a1.a2.cov <- marker1$a*marker2$a*(rho0^2/(2*n0)+rho.minus^2/(2*n.minus))
            b1.b2.cov <- rho.minus^2*marker1$b*marker2$b/(2*n.minus)+rho0*marker1$a*marker2$a/n0+rho.minus/n.minus##
            c1.c2.cov <- marker1$c*marker2$c*(rho0^2/(2*n0)+rho.plus^2/(2*n.plus))
            d1.d2.cov <- rho.plus^2*marker1$d*marker2$d/(2*n.plus)+rho0*marker1$c*marker2$c/n0+rho.plus/n.plus##
            
            ##diff para
            a1.b2.cov <- rho.minus^2*marker1$a*marker2$b/(2*n.minus)
            a1.c2.cov <- rho0^2*marker1$a*marker2$c/(2*n0)
            b1.d2.cov <- rho0*marker1$a*marker2$c/n0
            c1.d2.cov <- rho.plus^2*marker1$c*marker2$d/(2*n.plus)

            ###
            b1.a2.cov <- rho.minus^2*marker1$b*marker2$a/(2*n.minus)
            c1.a2.cov <- rho0^2*marker1$c*marker2$a/(2*n0)
            d1.b2.cov <- rho0*marker1$c*marker2$a/n0
            d1.c2.cov <- rho.plus^2*marker1$d*marker2$c/(2*n.plus)
              
            cov0 <- marker1$V.a*marker2$V.a*a1.a2.cov+marker1$V.b*marker2$V.b*b1.b2.cov+marker1$V.c*marker2$V.c*c1.c2.cov+marker1$V.d*marker2$V.d*d1.d2.cov+marker1$V.a*marker2$V.b*a1.b2.cov+marker1$V.a*marker2$V.c*a1.c2.cov+marker1$V.b*marker2$V.d*b1.d2.cov+marker1$V.c*marker2$V.d*c1.d2.cov +marker1$V.b*marker2$V.a*b1.a2.cov+marker1$V.c*marker2$V.a*c1.a2.cov+marker1$V.d*marker2$V.b*d1.b2.cov+marker1$V.d*marker2$V.c*d1.c2.cov
            
          }
        else
          {
            res.minus <- mu.s.cov(n.minus,rho.minus,marker1$s.minus,marker2$s.minus)
            res0 <- mu.s.cov(n0,rho0,marker1$s0,marker2$s0)
            res.plus <- mu.s.cov(n.plus,rho.plus,marker1$s.plus,marker2$s.plus)
            
            mu.minus.cov <- res.minus["mu.cov"]
            s.minus.cov <- res.minus["s.cov"]
            
            mu0.cov <- res0["mu.cov"]
            s0.cov <- res0["s.cov"]
            
            mu.plus.cov <- res.plus["mu.cov"]
            s.plus.cov <- res.plus["s.cov"]
            
            cov0 <- marker1$Y.mu.minus*marker2$Y.mu.minus*mu.minus.cov+marker1$Y.mu0*marker2$Y.mu0*mu0.cov+marker1$Y.mu.plus*marker2$Y.mu.plus*mu.plus.cov+marker1$Y.s.minus*marker2$Y.s.minus*s.minus.cov+marker1$Y.s0*marker2$Y.s0*s0.cov+marker1$Y.s.plus*marker2$Y.s.plus*s.plus.cov
          }
        return(cov0)
      }


    ##################Start ###################

    type <- match.arg(type)##allows incomplete input in the argument, e.g., "V" for "VUS" and "Y" for "Youden"
    
    if(!type%in%c("VUS","Youden")) stop("type must be VUS or Youden!!")
    
    dat.str <- class(dat)
    if(dat.str=="list" & paired) stop("Paired data (where markers measured on the same setof subjects) should be organized as a data frame with samples at row and group& marker measurements at column, not a list!")
    if(dat.str=="data.frame" & !paired) stop("unpaired data (where markers measured on different sets subjects) should be organized as a list (with each component for a marker), not a data frame!")

    ####check data structure
    if(dat.str=="list") K <- length(dat)
    if(dat.str=="data.frame") K <- ncol(dat)-1##the first column is group membership
    if(K<2) stop("must have >=2 markers!")

    alternative <- match.arg(alternative)
    data.name <- paste("Test ", deparse(substitute(dat))," on ", type,sep="")
    
    ###Construct the contrast matrix of K*(K-1) at row r, the r-th position is 1 and the (r+1) position is -1
    A <- matrix(0,K-1,K)
    diag(A) <- 1
    for(ii in 1:(K-1)) A[ii,ii+1] <- -1
       
    ###independent samples
    if(!paired)
      {
        ###Obtain point estimate on the summary measure and its associated variance
        markerID <- names(dat)
        res0 <- lapply(dat,function(marker.dat)
                       {
                         marker.list <- split(marker.dat[,2],marker.dat[,1])##make sure group at 1st column and marker measurements at 2nd column
                         
                         xx <- na.exclude(marker.list$'D-')
                         yy <- na.exclude(marker.list$D0)
                         zz <- na.exclude(marker.list$'D+')                 
                         
                         if(type=="VUS")
                           {
                             vus.res <- VUS(x=xx,y=yy,z=zz,method="Normal",p=p,q=q)
                             data.frame(estimate=vus.res$estimate,variance=vus.res$variance,vus.res$partialDeriv)
                           }
                         else 
                           {
                             youden.res <- Youden3Grp(x=xx,y=yy,z=zz,method="Normal")
                             data.frame(estimate=youden.res$estimate,variance=youden.res$variance,youden.res$partialDeriv)
                           }                                                                                               
                       })
        #browser()
        res0 <- do.call("rbind",res0)
        theta <- res0$estimate
        variance <- res0$variance
        
        ###construct the covariance matrix on theta estimates with variance on estimate at the diagonal and all other entries are 0
        Sigma <- diag(variance)
      }

    ###paired samples, need to calculate the covariance
    if(paired)
      {
        
       markerID <- names(dat)[-1]
       ###Obtain point estimate on the summary measure and its associated variance
        
        res0 <- apply(dat[,-1],2,function(marker.dat)
                       {
                         marker.list <- split(marker.dat,dat[,1])##make sure group at 1st column and marker measurements at 2nd column
                         
                         xx <- na.exclude(marker.list$'D-')
                         yy <- na.exclude(marker.list$D0)
                         zz <- na.exclude(marker.list$'D+')                 
                         
                         if(type=="VUS")
                           {
                             vus.res <- VUS(x=xx,y=yy,z=zz,method="Normal",p=p,q=q)
                             #partialDeriv <- vus.res$partialDeriv
                             data.frame(estimate=vus.res$estimate,variance=vus.res$variance,vus.res$partialDeriv)
                           }
                         else 
                           {
                             youden.res <- Youden3Grp(x=xx,y=yy,z=zz,method="Normal")
                             data.frame(estimate=youden.res$estimate,variance=youden.res$variance,youden.res$partialDeriv)
                           }                                                                                               
                       })
        #browser()##make sure the output is a data frame with the names
        res0 <- do.call("rbind",res0)
        theta <- res0$estimate
        variance <- res0$variance
        
        partialDeriv <- subset(res0,select=-c(estimate,variance))##a data frame of K row (each corresponding to a marker) and columns are (1) partial derivative w.r.t a, b, c, d: V.a, V.b, V.c, V.d (2) variance, covariance: a.var, b.var,c.var,d.var,a.b.cov,a.c.cov,b.d.cov,c.d.cov

        ###obtain within each diagnosis group correlation among  markers
        
        cor.minus <- cor(dat[dat[,1]=="D-",-1],use="pairwise.complete.obs")##a matrix of K*K, cor within D- group
        cor0 <- cor(dat[dat[,1]=="D0",-1],use="pairwise.complete.obs")
        cor.plus <- cor(dat[dat[,1]=="D+",-1],use="pairwise.complete.obs")
        
        ###construct the covariance matrix by placing the variance on theta estimates at the diagnonal term
        Sigma <- diag(variance)
        
        ####put at the entry (i,j) of Sigma the covariance of marker i and marker j     
        
        if(K==2)
          {
            ###get sample size with complete obs on each pair of marker and have group label
            temp.dat <- na.exclude(dat[,c(1,2,3)])
            n.minus <- sum(temp.dat[,1]=="D-")
            n0 <- sum(temp.dat[,1]=="D0")
            n.plus <- sum(temp.dat[,1]=="D+")

            #browser()
            current.cov <-  BetweenMarker.cov(marker1=partialDeriv[1,],  marker2=partialDeriv[2,],n.minus=n.minus,n0=n0,n.plus=n.plus,rho.minus=cor.minus[1,2],rho0=cor0[1,2],rho.plus=cor.plus[1,2],type=type)
            Sigma[1,2] <- Sigma[2,1] <- current.cov
            
          }else{
            for (i1 in 1:(K-1))
              {
                for(i2 in (i1+1):K)
                  {
                    temp.dat <- na.exclude(dat[,c(1,i1,i2)])
                    n.minus <- sum(temp.dat[,1]=="D-")
                    n0 <- sum(temp.dat[,1]=="D0")
                    n.plus <- sum(temp.dat[,1]=="D+")
                    current.cov <- BetweenMarker.cov(marker1=partialDeriv[i1,],  marker2=partialDeriv[i2,],n.minus=n.minus,n0=n0,n.plus=n.plus,rho.minus=cor.minus[i1,i2],rho0=cor0[i1,i2],rho.plus=cor.plus[i1,i2],type=type)
                    Sigma[i1,i2] <- Sigma[i2,i1] <- current.cov
                  }
              }
          }

      }##end of paired 
    #print(Sigma)
    ###calculate A*theta
    A.theta <- A%*%matrix(theta,ncol=1)
    ###Calcualte the Chi-square stat and associated P-value
    A.Sigma.A <- A%*%Sigma%*%t(A)
    inv.mat <- solve(A.Sigma.A)
    #print(inv.mat)
    chisq.stat <- t(A.theta)%*%inv.mat%*%A.theta
    names(chisq.stat) <- "Chi-square"
    
    df <- K-1
    names(df) <- "df"
    chisq.p.value <- pchisq(chisq.stat,df=K-1,lower.tail=FALSE)


    names(mu) <- paste("diff in ",type,sep="")
    
    ####For two markers, do z-test
    estimate <- theta
    names(estimate) <- paste(paste(type," of ",sep=""),markerID,sep="")
    if(K==2)
      {
        diff <- theta[1]-theta[2]
        diff.se <- sqrt(as.numeric(A.Sigma.A))
        
        ###alternative=="two.sided"
        z.stat <- (diff-mu)/diff.se
        
        names(z.stat) <- "Z-stat"
        
        if(alternative=="two.sided")
          {
            z.p.value <- 2*pnorm(-abs(z.stat),lower.tail=TRUE)
            #ci <- CI.normal(diff,abs.diff.se,alpha=1-conf.level)
            alpha <- 1-conf.level
            cint <- qnorm(1-alpha/2)*diff.se
            cint <- diff+c(-cint,cint)
          }
        else if (alternative=="less")
          {
            z.p.value <- pnorm(z.stat,lower.tail=TRUE)
            cint <- c(-Inf,diff+qnorm(conf.level)*diff.se)
          }else if (alternative=="greater")
            {
              z.p.value <- pnorm(z.stat,lower.tail=FALSE)
              cint <- c(diff-qnorm(conf.level)*diff.se,Inf)
            }
        attr(cint,"conf.level") <- conf.level
        
        method <- "normal-test"
        
        mean.sd <- c(mean=0,sd=1)

        rval <- list(statistic=z.stat,p.value=z.p.value,parameter=mean.sd,estimate=estimate,conf.int=cint,method=method,type=type,null.value=mu,alternative=alternative,data.name=data.name,Sigma=Sigma)
        
      }else{
        method <- "Chi-Sqaure test"
        rval <- list(statistic=chisq.stat,p.value=chisq.p.value,parameter=df,estimate=estimate,conf.int=NA,method=method,type=type,null.value=mu,alternative=alternative,data.name=data.name,Sigma=Sigma)
      }
    class(rval) <- "htest"
    return(rval)
    
  }
