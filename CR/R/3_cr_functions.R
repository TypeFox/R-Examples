
curerate<-function(rho = 0,numreps,cureobs,curerx,medobs,medrx,acrate,
              probrx,actime,futime,info,crits,alpha = .025){
# ********************************************************************
# This function performs the sample size and power calculations for the 
# set of all pairs (actime[i],futime[j]), i=1,...,numac, j=1,...,numfu.
# 
# Output - an object of class "CureRate".
#    

    # Error checking in the parameter entries
    #****************************************

    error_lst<-list()
    error_lst[[1]]<-"'numreps' - Number of replicates must be an integer > 0."
    error_lst[[2]]<-"'acrate' - Accrual rate must be a scalar > 0."
    error_lst[[3]]<-"'probrx' - Prob assigned to Rx must be a scalar in [0,1]."
    error_lst[[4]]<-"'alpha' - Significance level must be a scalar in (0,1)."
    error_lst[[5]]<-"'rho' - rho should be a numeric scalar; it specifies the G-rho test."
    error_lst[[6]]<-"'cureobs' - Obs cure rate must be a scalar in [0,1]." 
    error_lst[[7]]<-"'curerx' - Rx cure rate must be a scalar in [0,1]."
    error_lst[[8]]<-"'medobs' - Obs median survival must be a scalar > 0."
    error_lst[[9]]<-"'medrx' - Rx median survival must be a scalar > 0."
    error_lst[[10]]<-"'actime' - accrual durations must be numeric and distinct."
    error_lst[[11]]<-"'actime' - accrual duration must be > 0."
    error_lst[[12]]<-"'actime' - accrual durations must increase."
    error_lst[[13]]<-"'futime' - follow-up durations must be numeric and distinct."
    error_lst[[14]]<-"'futime' - follow-up duration must be > 0."
    error_lst[[15]]<-"'futime' - follow-up durations must increase."
    error_lst[[16]]<-"'info' - information times must be entered as a numeric vector."
    error_lst[[17]]<-"'info' - information times must be distinct."
    error_lst[[18]]<-"'info' - information times must be in (0,1]."
    error_lst[[19]]<-"'info' - information times must increase."
    error_lst[[20]]<-"'crits' - critical values must be entered as a numeric vector."
    error_lst[[21]]<-"'crits' - critical values must be distinct."
    error_lst[[22]]<-"'crits' - critical values must be > 0."
    error_lst[[23]]<-"'crits' - critical values must decrease."
    error_lst[[24]]<-"lengths of 'info' and 'crits' must agree."


    error.curerate <- function(rho = 0,numreps,cureobs,curerx,medobs,medrx,
            acrate,probrx,actime,futime,info,crits,alpha){
    ## Checks for the errors in the entered input parameters of curerate()

        error.lst <- list();
        if (!(is.numeric(numreps) && length(numreps)==1 && numreps>0)) return(1)
        if (!(is.numeric(acrate) && length(acrate)==1 && acrate>0)) return(2) 
        if (!(is.numeric(probrx) && length(probrx)==1 && probrx>=0 && probrx<=1)) return(3) 
        if (!(is.numeric(alpha) && length(alpha)==1 && alpha>0 && alpha<1)) return(4) 
        if (!(is.numeric(rho) && length(rho)==1) ) return(5)
        if (!(is.numeric(cureobs) && length(cureobs)==1 && cureobs>=0 && cureobs<=1)) return(6)
        if (!(is.numeric(curerx) && length(curerx)==1 && curerx>=0 && curerx<=1)) return(7)
        if (!(is.numeric(medobs) && length(medobs)==1 && medobs>0)) return(8) 
        if (!(is.numeric(medrx) && length(medrx)==1 && medrx>0)) return(9)
        numac <- length(actime)
        if (!(is.numeric(actime) && numac == length(unique(actime)))) return(10)
        o<-order(actime)
        if(actime[o][1]<=0) return(11)
if (sum(actime != actime[o]) != 0) return(12)
#        for(i in 1:numac){
 #           if(actime[i] > actime[o][i]) return(12)
  #      }
        numfu <-length(futime)    
        if (!(is.numeric(futime) && numfu == length(unique(futime)))) return(13)
        o=order(futime)
        if (futime[o][1]<=0) return(14)
if (sum(futime != futime[o]) != 0) return(15)
#        for(i in 1:numfu){
 #           if(futime[i] > futime[o][i]) return(15)
  #      }
        if (!is.numeric(info)) return(16)
        info<-as.vector(info)
        numinfo<-length(info)
        if(numinfo>length(unique(info))) return(17)
        o=order(info)
        if(info[o][1]<=0 || info[o][numinfo]>1) return(18)
if (sum(info != info[o]) != 0) return(19)
#        for(i in 1:numinfo){
 #           if(info[i] > info[o][i]) return(19)
  #      } 
        if(!is.numeric(crits)) return(20)
        crits<-as.vector(crits)
        numcrits<-length(crits)
        if(numcrits > length(unique(crits))) return(21)
        o=order(crits,decreasing = TRUE)
        if(crits[o][numcrits]<=0) return(22)
if (sum(crits != crits[o]) != 0) return(23)
#        for(i in 1:numcrits){
 #           if(crits[i] < crits[o][i]) return(23)
  #      }
        if(numinfo != numcrits) return(24)
        return(0)
    }## End of error.curerate()

    rc<-error.curerate(rho,numreps,cureobs,curerx,medobs,medrx,
            acrate,probrx,actime,futime,info,crits,alpha)
    if(rc>=1)
        stop(error_lst[[rc]])

    
    ##### Echoing inputs to file 


    survtest <- function(x,d,z,rho){
    # 
    #### ***** d and z MUST be coded 0 and 1 *****
    #
    # This function implements the G-rho family of Harrington and Fleming (1982), 
    # tests, with weights on each death of S(t)^rho, where S is the Kaplan-Meier 
    # estimate of survival. With rho = 0 this is the log-rank or Mantel-Haenszel 
    # test, and with rho = 1 it is equivalent to the Peto & Peto modification of 
    # the Gehan-Wilcoxon test, for DATA WITH TWO GROUPS AND NO TIES.
    # 
    # Input
    #     x[] - vector of event times with NO TIES;
    #     d[] - vector of censoring indicators must have values 0 and 1
    #           (d=0 - censored, d=1 failure);
    #     z[] - vector of group indicators must have values 0 and 1 
    #           (z=0 - control group, z=1 - treatment group);
    #     x[],d[], and z[] must be of same length, and 
    #     x[] does NOT have to be sorted in ascending order;
    #     rho - scalar, the parameter "rho" for the G-rho family of tests.
    #
    # Output 
    #    test - scalar, the value of the G-rho test statistic.
    #

        n=length(x)
        o=order(x)
        d<-d[o]
        z<-(1-z)[o]
        N1frac <- rev(cumsum(rev(z))) / n:1
        obs1 <- d * z 
        exp1 <- d * N1frac
        Srho=c(1,cumprod(1-1/n:2*d[-n])^rho)  
        obs1=Srho%*%obs1
        exp1=Srho*exp1
        var  <- (exp1 * (1-N1frac))%*%Srho
        ## factor (Nj-Oj)/(Nj-1) in var terms is not needed since Oj is 0 or 1
        exp1 <- sum(exp1)
        #score=obs1-exp1
        test = 0
        if(var > 0.0)  test = (obs1-exp1)/sqrt(var)
        return(test)

    } ## End of "survtest()"


    cr<-function(ac_time,fu_time,rho,numreps,cureobs,curerx,medobs,medrx,
            acrate,probrx,info,crits){
    # **********************************************************
    # This functions performs the main part of the program, the sample 
    # size and power calculations for a single pair (actime,futime).
    #

        lambda1<-log(2.)/medobs
        lambda2<-log(2.)/medrx
    
        numlook <- length(info)
        power = vector(mode="numeric",length=numlook)
        deaths = matrix(NA,nrow=numreps,ncol=numlook)
        timept = matrix(NA,nrow=numreps,ncol=numlook)
        reject = matrix(NA,nrow=numreps,ncol=numlook)
        numtest = matrix(NA,nrow=numreps,ncol=numlook)

        ## Compute expected number of observations
        numobs<- ac_time*acrate

        ## Compute expected number of deaths
        expected = (1-probrx)*(1. - cureobs)*acrate*
                   (ac_time-(medobs/log(2.))*
                   (exp(-(log(2.)/medobs)*fu_time)-
                   exp(-(log(2.)/medobs)*(fu_time+ac_time))))+
                   probrx*(1. - curerx)*acrate*
                   (ac_time-(medrx/log(2.))*
                   (exp(-(log(2.)/medrx)*fu_time)-
                   exp(-(log(2.)/medrx)*(fu_time+ac_time))))


        ## Compute the expected number of deaths at the interim looks
        ## based on the informational times
        expinfo = floor(expected * info + 0.5)

        ## Check whether or not the expected number of deaths increase
        ## from an interim look to the next one. If not, STOP!
        if(length(unique(expinfo)) != length(expinfo))
            stop("info times too close")

        for(n in 1:numreps){ 
            treat <- rbinom(numobs,1,probrx)
            lambda <- lambda1+(lambda2-lambda1)*treat               
            survtime <- rexp(numobs)/lambda 
            probcure <- cureobs+(curerx-cureobs)*treat
            cure <- rbinom(length(probcure),1,probcure)
            accrual <- runif(numobs) * ac_time

            ## FINDING time points for analysis based on information
            x = accrual + survtime
       
            ## sorting survival + accrual times
            o = order(x)
            x = x[o]

            ## finding event points for analysis based on information        
            notcured  = which(!cure[o])
            numnotcured = length(notcured)

            if(expinfo[numlook] <= numnotcured)  events = expinfo
            else events = c(expinfo[expinfo < numnotcured],numnotcured)
            nlook1 = length(events)

            ## finding time points for analysis
            time = x[notcured[events]]

            deaths[n,1:nlook1] = events
            timept[n,1:nlook1] = time
            numtest[n,1:nlook1] = 1

            for(m in 1:nlook1){
                ## establishing censoring and in sample indictors
                survcens =  (!cure) * ((survtime+accrual) <= time[m])
                insample = (accrual <= time[m])

                ## extracting data available at analysis time point
                d = survcens[insample]
                x = survtime[insample]*d + (time[m] - accrual[insample])*(1-d)
                z = treat[insample]

                ## test survival curve differences
                test = survtest(x,d,z,rho) 
                if(test >= crits[m]){
                    reject[n,m] = 1
                    break
                }
            } ## End of "for(m in 1:nlook1)"
        } ## End of "for(n in 1:numreps)"

        numtestSums = colSums(numtest,na.rm=TRUE)
        rejectSums = colSums(reject,na.rm=TRUE)
        power = rejectSums/numtestSums
        beta <- sum(power)
        deaths_sum <- colSums(deaths,na.rm=TRUE)/numtestSums
        timept_sum <- colSums(timept,na.rm=TRUE)/numtestSums

        ## Output of R-fc cr()
        rlst_cr <- list(numobs=numobs,numreps=numreps, timept=timept_sum,
                      deaths=deaths_sum,power=power,beta=beta,reject=rejectSums,
                      numtest=numtestSums,expected=expected,events=events,
                      time=time,len_insample=length(insample),len_x=length(x),
                      len_d=length(d),len_z=length(z),x=x,d=d,z=z,test=test)

        return(rlst_cr)

    }## End of R-fc "cr()"


    ## set the test name
    if(rho==0) test_name <- "Test Type: Logrank Test"
    if(rho==1) test_name <- "Test Type: Pete-Peto Wilcoxon Test"
    if(rho==-1) test_name <- "Test Type: Gray-Tsiatis Test"
    if(rho != 0 && rho != 1 && rho != -1) 
        test_name <- paste("Test Type: G-rho, rho=",rho)

    numlook <- length(info)
    numac <- length(actime)
    numfu <- length(futime)

    crout=new("CureRate",
        cureobs=cureobs,medobs=medobs,curerx=curerx,medrx=medrx, 
        info=info,crits=crits,alpha=alpha,
        actime=actime,futime=futime,rho=rho,
        acrate=acrate,probrx=probrx,numreps=as.integer(numreps),
        numobs=matrix(0, nrow=numac, ncol=numfu),
        testname=test_name,
        timept=array(0,dim=c(numac,numfu,numlook)),
        deaths=array(0,dim=c(numac,numfu,numlook)),
        power=array(0,dim=c(numac,numfu,numlook)),
        beta=matrix(0,nrow=numac, ncol=numfu),
        printflag=1L
    )

    for(i in 1:numac){
        for(j in 1:numfu){
            rlst_cr <- cr(actime[i],futime[j],rho,numreps,cureobs,curerx,
                          medobs,medrx,acrate,probrx,info,crits)
            crout@numobs[i,j] <- rlst_cr$numobs
            crout@timept[i,j,] <- rlst_cr$timept
            crout@deaths[i,j,] <- rlst_cr$deaths
            crout@power[i,j,] <- rlst_cr$power
            crout@beta[i,j] <- rlst_cr$beta
        }   
    }

    return(crout)

}## End of R-fc "curerate()"


showcr<-function(cr,full.results=FALSE, indac=0, indfu=0){
# Customize the display of an object of class "CureRate".

    if(full.results==TRUE) cr@printflag=as.integer(2)
    else{
        if(indac==0 || indfu==0) cr@printflag=as.integer(1)
        else{
            cr@indac=as.integer(indac)
            cr@indfu=as.integer(indfu)
            cr@printflag=as.integer(3)
        }
    }
    show(cr)

    return(invisible(cr))
}


