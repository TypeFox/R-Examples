VUS <- function(x,y,z,method=c("Normal","NonPar"),p=0,q=0,alpha=0.05,NBOOT=100,subdivisions=50000,lam.minus=1/3,lam0=1/3,lam.plus=1/3,typeIerror=0.05,margin=0.05,FisherZ=FALSE,optimalCut=TRUE,cut.seq=NULL,optimize=FALSE,...)
  {
    ######################################################################################################################################################
    #######A wrapper function provides: VUS (volume under ROC surface), variance and (1-alpha)*100% under both normal distribution assumption and 
    #######                       nonparametric VUS estimate and associated bootstrap variance and (1-alpha)*100% basic quatile
    ######1. Input
    ###(1)x,y,z:vectors are the test marker measurements for the 3 ordinal group (D-/D0/D+ or healthy/medium diseased/diseased)
    ### NOTE:For estimate under normal, mean values of the three groups should increase with severity of disease from D- to D0 and D+,
    ####       if reverse, simply take negated x,y,z as inputs instead
    ####(2)method: "Normal"
    ###(3)p,q:used for VUS under normal assumption. p=minimum desired specificity for the D- group, q=minimum desired sensitivity for D+ group, 0<=p,q<1
    ###(4)alpha:for the nonparametric estimate, provide (1-alpha)% boostrap basic quantile CI
    ###(5)NBOOT:for the nonparametric estimate, # of bootstrap samples to draw to obtain bootstrap variance and CI
    ###(6) optimalCut=T: the functin will return the optimal cut point from VUS analyses
    ###(7) cut.seq: the sequence of values from which the optimal cut-point will be selected from, by default=NULL, will use the unique values of the collection of x,y,z
    ####(8) FisherZ
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

    ##########################################################################################################################################################
    
    #if(!(method%in%c("Normal","NonPar"))) stop("The method choice should be  Normal/NonPar!")

    method <- match.arg(method)
    
    if(method=="Normal")
      {
        temp.res <- Normal.VUS(x=x,y=y,z=z,p=p,q=q,alpha=alpha,subdivisions=subdivisions,lam.minus=1/3,lam0=1/3,lam.plus=1/3,typeIerror=0.05,margin=0.05,FisherZ=FisherZ,...)
        partialDeriv <- temp.res$partialDeriv##a data.frame
      }else if(method=="NonPar")
      {
        temp.res <- NonParametric.VUS(x=x,y=y,z=z,FisherZ=FisherZ)##estimate VUS
        boot.res <- NonParametric.VUS.var(x=x,y=y,z=z,alpha=alpha,NBOOT=NBOOT,FisherZ=FisherZ)
        temp.res$variance <-boot.res$var0
        temp.res$CI <- boot.res$CI
        temp.res$sampleSize <- NA
        partialDeriv <- NA
      }else stop("ERROR: The method choice should only be Normal/NonPar!")

    if(optimalCut)
      {
        cut.point <- VUS.CutPoint(x=x,y=y,z=z,cut.seq=cut.seq)
        ###correct classification probability
        sp.sm.se <- Sp.Sm.Se(x=x,y=y,z=z,t.minus=cut.point[1],t.plus=cut.point[2])##Sp=specificity Pr(x<=t-|D-),##Se=sensitivity Pr(z>=t+|D+) ## Sm=Pr(t-<=y<=t+|D0)  
      }else
    {
      cut.point <- NA
      sp.sm.se <- NA
    }

   
    out0 <- list(type="VUS",method=method,dat=temp.res$dat,dat.summary=temp.res$dat.summary,estimate=temp.res$estimate,variance=temp.res$variance,CI=temp.res$CI,cut.point=cut.point,classify.prob=sp.sm.se,sampleSize=temp.res$sampleSize,alpha=alpha,typeIerror=typeIerror,margin=margin,partialDeriv=partialDeriv)
    
    class(out0) <- "DiagTest3Grp"
    
    return(out0)
  }

