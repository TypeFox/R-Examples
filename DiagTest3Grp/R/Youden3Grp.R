#################################################################################################################
#The wrapper function to obtain point estimate and variance on Youden index estimation
# by calling functions Youden3Grp.PointEst() and Youden3Grp.Variance.Normal() and Youden3Grp.Variance.Bootstrap()
#
###################################################################################################################

Youden3Grp <- function(x,y,z, method=c("Normal","TN","EMP","KS","KS-SJ"),randomStart.N=1,optim.method=NULL,t.minus.start=NULL,t.plus.start=NULL,lam.minus=1/3,lam0=1/3,lam.plus=1/3,typeIerror=0.05,margin=0.05,NBOOT=10,seed.seq=NULL,alpha=0.05,FisherZ=FALSE,...)
  {

     ###Inputs:
    ###x,y,z  diagnostic test marker measurements in the D- (healthy), D0,D+(disease) groups,
    ###NOTE:!!!! x, y, z must have increasing mean, i.e., mean(x)<mean(y)<mean(z)
    ###method: "Normal"/"TN"/"EMP"/"KS"/"KS-SJ"
      ###For EMP/KS, we have to do numerical optimization for estimates of Youden index and cutoff points. Thus it will depend on initial values when data have a wide range
      ####If using TN, need all values to be positive,if original input data are not all positive, a small positive constant will be added, which will not
      ####               affect the youden index estimation but will on optimal cutoff, we adjust this back at the end for cut off estimation
    
    ###optim.method: Optimization methods have these options,method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),by default l-BFGS-B"
    ###                and BFGS needs gradient function provided, In KS, the gradient is easy to be provided, thus may be better use BFGS by providing gr=KS.gradient()
    
    ####we do random start (multiple starting t.minus and t.plus for EMP and KS method in case of bumpy goal function
    
    ####randomStart.N=1, only 1 starting value for t.minus and t.plus , set randomStart.N=10 if want to get the t.minus and t.plus at the maximum value
    ####start.t.minus,start.t.plus:starting values for t.minus and t.plus for use in the optimization steps for emp and KS methods,
    ###                                              it missing or NULL, will use normal solutions or multiple start (if randomStart.N>1) points
    ####...: can take the other arguments as in optim() function

    ## lam.minus, lam0,lam.plus: for sample size calculation, the expected proportion of samples in the D-, D0 and D+ group, which can be equal or not
    ##NBOOT:# of bootstrap iterations
    ###seed.seq: the random seed for bootstrapping procedure to obtain bootstrapping variance and percentile confidence interval if methods other than "Normal" is adopted
    ###alpha: significance level to obtain CI
    ###FisherZ: see VUS.R
    
     ######2. Ouput: a DiagTest3Grp object, consisting of a list of components
    ####
    ####(1)type: a character value, "VUS" for VUS and "Youden" for the extended Youden index, indicating which summary measure is outputted
    ####(2)method: a character value. For VUS, method can be "Normal" or "NonPar" (nonparametric) ; for Youden index, choices are Normal/TN/EMP/KS/KS-SJ, indicating which method is used to estimate the summary measure
    ####(3)dat: a list. Three components are "x","y","z", recording the marker measurements under D-, D0 and D+ respectively
    ####(4)dat.summary: a data frame with 3 rows (D-, D0,D+) and 3 columns (number of observations,mean, SD).
    ####(5)estimate: a numerical value. Point estimate for the summary measure
    ####(6)variance: variance on the summary measure estimate. For method="Normal", output normal variance; for other methods output variance from bootstrapping    
    ####(7)CI: a named numeric vector of length 2. confidence interval,with name "2.5%", 97.5%" if significance level is set to be 5%. For both VUS and the Youden index, when method="Normal", the CI is normal CI while bootstrap method was used under other methods.
    ####(8)cut.point: a named numeric vector of length 2. optimal cut-points with name "t.minus" for lower optimal cut point and "t.plus" for upper optimal point
    ####(9)classify.prob: a named numeric of 3 values. Estimates on the three group correct classification probabilities:Sp=specificity=Pr(x<=t-|D-), Se=sensitivity=Pr(z>=t+|D+), Sm=SPr(t-<=y<=t+|D0). For VUS, it's empirical estimation. For Youden index, depending on method adopted for the Youden index estimate, the three probabilities will be estimated using specified method.
    ####(10)sampleSize: the sample size to estimate the summary measure within given margin of error and type-I error rate
    ####(11)alpha: the significance level for the CI computation
    ####(12)typeIerror: type-I error rate
    ####(13)margin: the margin of errors (precision) to calculate the sample size,see SampleSize.VUS and SampleSize.Youden3Grp
        
    #####REMOVE (14) and (15)
    ####(14)transform.dat: a list similar to dat but with an additional component lambda, recording the estimated lambda parameter in Box-cox transformation . Specifically for method TN applied to obtain Youden index, recording the data transformed according with Box-Cox transformation.
    ####(15)transform.dat.summary: similar to dat.summary but on the Box-Cox transformed data as described in item transform.dat.

    method <- match.arg(method)
    
    ######Youden index, cut-point and sample size point estimate
    
    pointEst <- Youden3Grp.PointEst(x=x,y=y,z=z, method=method,randomStart.N=randomStart.N,optim.method=optim.method,t.minus.start=t.minus.start,t.plus.start=t.plus.start,lam.minus=lam.minus,lam0=lam0,lam.plus=lam.plus,typeIerror=typeIerror,margin=margin,FisherZ=FisherZ,...)
    
    ######Variance Computation
    if(method=="Normal")
      {
        temp.res <- Youden3Grp.Variance.Normal(x=x,y=y,z=z,alpha=alpha)
        partialDeriv <- temp.res$partialDeriv##the partial derivatives and others to be used for test markers 
      }
    else if(method %in%c("TN","EMP","KS","KS-SJ"))
      {
        temp.res <- Youden3Grp.Variance.Bootstrap(x=x,y=y,z=z,method=method,seed.seq=seed.seq,randomStart.N=randomStart.N,NBOOT=NBOOT,alpha=alpha,t.minus.start=t.minus.start,t.plus.start=t.plus.start)
        partialDeriv <- NA
      }
    else stop("methods should be Normal/TN/EMP/KS/KS-SJ")

    if(!FisherZ) CI <- temp.res$youden.CI else CI <- temp.res$youden.z.CI
    
    out0 <- list(type="Youden",method=method,dat=pointEst$dat,dat.summary=pointEst$dat.summary,estimate=ifelse(FisherZ,pointEst$est$youden.z,pointEst$est$youden),variance=ifelse(FisherZ,temp.res$var.youden.z,temp.res$var.youden),CI=CI,cut.point=c(t.minus=pointEst$est$t.minus,t.plus=pointEst$est$t.plus),classify.prob=c(Sp=pointEst$est$Sp,Sm=pointEst$est$Sm,Se=pointEst$est$Se),sampleSize=pointEst$est$sampleSize,alpha=alpha,typeIerror=typeIerror,margin=margin,partialDeriv=partialDeriv)####cut.point in original data scale thoug t.minus.TN and t.plus.TN record the optimal cut-points after Box-Cox transformation 
    ###note:var.t.minus, var.t.plus can be provided too but will not be outputted here for the DiagTest3Grp class.
    
    class(out0) <- "DiagTest3Grp"
           
    return(out0)

  }

