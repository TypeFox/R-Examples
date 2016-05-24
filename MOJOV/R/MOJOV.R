MOJOV.read<-function(x=NULL,genofile=NULL,phenofile=NULL,indelfile=NULL,
                     header=T,column=1,...)
{
    if(is.null(genofile)&class(x)!="MOJOV")
        stop("You must input a correct genotype file name and path.")
    if(is.null(phenofile))
        stop("You must input a correct phenotype file name and path.")
    if(class(x)=="MOJOV")
    {
        genoFile=genoFile(x)
        genoIID=genoIID(x)
        genoChr=genoChr(x)
        genoPosi=genoPosi(x)
        genoVariant=genoVariant(x)
        genoH=genoH(x)
        genoGene=genoGene(x)
    }
    else
    {
        x1<-read.csv(genofile,header=T,...)
        names(x1)<-c("IID","chromosome","Position","SNP","Mut","Gene")
        if(!is.null(indelfile)){
            if(file.access(indelfile)>=0){
                x2<-read.csv(indelfile,header=T,...)
                genofile<-c(genofile,indelfile)
                names(x2)<-c("IID","chromosome","Position","SNP","Mut","Gene")
                x<-rbind(x1,x2)
            }else{
                warning("Your indel file has some error, it doesn't work.")
            }
        }else{
            x=x1
        }
        g<-x[,5]
        gsub("het","Het",g)
        gsub("hom","Hom",g)
        geno<-ifelse(g=="het",1,2)
        geno<-as.numeric(geno)
        genoFile=genofile
        genoIID=as.character(x[,1])
        genoChr=as.character(x[,2])
        genoPosi=as.numeric(x[,3])
        genoVariant=as.character(x[,4])
        genoH=geno
        genoGene=as.character(x[,6])
    }
    y<-read.csv(phenofile,header=header,...)
    MOJOVObj<-new("MOJOV",genoFile=genoFile,genoIID=genoIID,
                  genoChr=genoChr,genoPosi=genoPosi,
                  genoVariant=genoVariant,genoH=genoH,
                  genoGene=genoGene,
                  phenoFile=phenofile,phenoIID=as.character(y[,1]),
                  phenoGender=as.numeric(y[,3]),phenoAge=as.numeric(y[,2]),
                  phenoHeight=as.numeric(y[,4]),phenoWeight=as.numeric(y[,5]),
                  phenoLabel=names(y)[5+column],phenoData=as.numeric(y[,5+column]),
                  phenoColumn=column,phenoIndNum=dim(y)[1])
    return(MOJOVObj)
}


MOJOV.analysis<-function(x=NULL,MAF=0.01,ROI="scan",savefile=NULL,
                         codeMethod=c("Proportion","Indicator","ChuanhuaXin"),
                         weightMethod=NULL,
                         testMethod=c("FTest","WaldTest","LRT","Sandwich","all"),...)
{
    if(class(x)!="MOJOV")
        stop("Must specify a MOJOV object.")
    codeMethod<-match.arg(codeMethod)
    weightMethod<-match.arg(weightMethod)
    testMethod<-match.arg(testMethod)
    
    genoIID=genoIID(x)
    genoVariant=genoVariant(x)
    genoH=genoH(x)
    genoGene=genoGene(x)
    genoChr=genoChr(x)
    genoPosi=genoPosi(x)
    phenoIID=phenoIID(x)
    phenoIndNum=phenoIndNum(x)
    adjustData=adjustData(x)    
    varTot=unique(genoVariant)
    varMAF=MAF
       
    if(is.null(ROI))
    {
        tmp<-MOJOV.genoMatrix(genoIID=genoIID,genoVariant=genoVariant,
                              genoH=genoH,phenoIID=phenoIID,MAF=varMAF,...)
        varFreq<-apply(tmp,2,sum)/(phenoIndNum*2)
        tmp<-MOJOV.genoVector(x=tmp,y=adjustData,codeMethod=codeMethod,
                              weightMethod=weightMethod,...)
        stat<-MOJOV.linearRegAnalysis(x=tmp,y=adjustData,
                                      testMethod=testMethod,...)
        resultPvalue=unlist(stat)
        cat("P-value:",stat)
        regionType="Whole"
        regionChr=unique(genoChr)
        regionPStart=numeric()
        regionPStop=numeric()
        regionGene=genoVariant
        regionFile=character()
    }
    else
    {
        if(ROI%in%genoGene)
        {
            regionGene=ROI
            index<-(1:length(genoGene))[genoGene==regionGene]
            tmp<-MOJOV.genoMatrix(genoIID=genoIID[index],genoVariant=genoVariant[index],
                                  genoH=genoH[index],phenoIID=phenoIID,MAF=varMAF,...)
            varFreq<-apply(tmp,2,sum)/(phenoIndNum*2)
            tmp<-MOJOV.genoVector(x=tmp,y=adjustData,codeMethod=codeMethod,
                                  weightMethod=weightMethod,...)
            return(lm(adjustData~tmp))
#This step has been edited. It will return a lm object directly.
            stat<-MOJOV.linearRegAnalysis(x=tmp,y=adjustData,
                                          testMethod=testMethod,...)
            resultPvalue=unlist(stat)
            cat("P-value:",stat)
            regionType="Gene"
            regionChr=genoChr[index][1]
            regionPStart=min(genoPosi[index])
            regionPStop=max(genoPosi[index])
            regionFile=character()
        }
        else
        {
            if(file.exists(ROI))
            {
                regionFile=ROI
                r<-read.csv(regionFile,header=FALSE,...)
                if(dim(r)[2]!=4)
                    stop("ROI file has an incorrect format.")
                pb<-txtProgressBar(min=1,max=length(r[,2]),style=3)
                pbi=1                
                regionType="Position"
                regionChr=as.character(r[,2])
                regionPStart=r[,3]
                regionPStop=r[,4]
                regionGene=as.character(r[,1])
                if(length((regionPStop>regionPStart)==FALSE)>0)
                    stop("Stop should greater than start position for ROI file.")
                stat<-numeric()
                for ( ii in 1:length(regionChr))
                {
                    index<-genoChr==regionChr[ii]
                    IID_tmp<-genoIID[index]
                    Variant_tmp<-genoPosi[index]
                    H_tmp<-genoH[index]
                    Posi_tmp<-genoVariant[index]
                    
                    index<-(Posi_tmp>regionPStart[ii]&Posi_tmp<regionPStart[ii])
                    IID_tmp<-IID_tmp[index]
                    Variant_tmp<-Variant_tmp[index]
                    H_tmp<-H_tmp[index]
                    
                    index<-genoPosi[genoVariant==ii]
                    tmp<-MOJOV.genoMatrix(genoIID=IID_tmp,genoVariant=Variant_tmp,
                                          genoH=H_tmp,phenoIID=phenoIID,MAF=varMAF,...)
                    varFreq<-apply(tmp,2,sum)/(phenoIndNum*2)
                    tmp<-MOJOV.genoVector(x=tmp,y=adjustData,codeMethod=codeMethod,
                                          weightMethod=weightMethod,...)
                    stat<-rbind(stat,MOJOV.linearRegAnalysis(x=tmp,y=adjustData,
                                                             testMethod=testMethod,...))
                    setTxtProgressBar(pb,pbi)
                    pbi<-pbi+1
                }
                resultPvalue=stat
            }
            else
            {
                if(ROI=="scan")
                {
                    geneTmp=unique(genoGene)
                    pb<-txtProgressBar(min=1,max=length(geneTmp),style=3)
                    pbi=1 
                    stat<-numeric()
                    for( ROIT in geneTmp)
                    {
                        regionGene=ROIT
                        index<-(1:length(genoGene))[genoGene==regionGene]
                        tmp<-MOJOV.genoMatrix(genoIID=genoIID[index],genoVariant=genoVariant[index],
                                              genoH=genoH[index],phenoIID=phenoIID,MAF=varMAF,...)
                        varFreq<-apply(tmp,2,sum)/(phenoIndNum*2)
                        tmp<-MOJOV.genoVector(x=tmp,y=adjustData,codeMethod=codeMethod,
                                              weightMethod=weightMethod,...)
                        stat<-rbind(stat,MOJOV.linearRegAnalysis(x=tmp,y=adjustData,
                                                                 testMethod=testMethod,...))
                        setTxtProgressBar(pb,pbi)
                        pbi<-pbi+1
                    }
                    resultPvalue=as.vector(stat)
                    regionType=geneTmp
                    regionChr=character()
                    regionPStart=numeric()
                    regionPStop=numeric()
                    regionFile=character()
                }
                else
                    stop("ROI should be a gene symbol from your file or a ROI file path.")
            }
        }
    }
    
    varRare<-varTot[varFreq<=varMAF]
    if(!is.null(savefile))
    {
        stat<-matrix(stat,,9,byrow=T)
        write.csv(x=data.frame(stat,phenotype=rep(phenoLabel(x),dim(as.matrix(stat))[1])),file=savefile,row.names=FALSE)
    }
    
    if(is.null(weightMethod))
        weightMethod=character()
    MOJOVObj<-new("MOJOV",varMAF=varMAF,varFreq=varFreq,varRare=varRare,
                  regionFile=regionFile,regionType=regionType,
                  regionChr=regionChr,regionPStart=regionPStart,
                  regionGene=regionGene,analyCode=codeMethod,
                  analyWeighted=weightMethod,
                  resultMethod=testMethod,resultPvalue=resultPvalue,
                  adjustAuto=adjustAuto(x),adjustData=adjustData(x),
                  adjustGender=adjustGender(x),adjustPower=adjustPower(x),
                  adjustPowerPvalue=adjustPowerPvalue(x),
                  adjustTerms=adjustTerms(x),
                  adjustPvalue=adjustPvalue(x),
                  genoFile=genoFile(x),genoIID=genoIID(x),
                  genoChr=genoChr(x),genoPosi=genoPosi(x),
                  genoVariant=genoVariant(x),genoH=genoH(x),
                  genoGene=genoGene(x),phenoFile=phenoFile(x),
                  phenoIID=phenoIID(x),phenoGender=phenoGender(x),phenoAge=phenoAge(x),
                  phenoHeight=phenoHeight(x),phenoWeight=phenoWeight(x),
                  phenoLabel=phenoLabel(x),phenoData=phenoData(x),
                  phenoColumn=phenoColumn(x),phenoIndNum=phenoIndNum(x))
    return(MOJOVObj)
}

MOJOV.genoMatrix<-function(genoIID=NULL,genoVariant=NULL,
                           genoH=NULL,phenoIID=NULL,MAF=NULL)
{
    if(is.null(genoIID))
        stop("You should specify your genotype ID infomation.")
    if(is.null(genoVariant))
        stop("You should specify your genotype varaint lablels.")
    if(is.null(genoH))
        stop("You should specify your genotype varaint infomation.")
    if(is.null(phenoIID))
        stop("You should specify your phenotype ID infomation.")
    variantLabels<-unique(genoVariant)
    variantNum<-length(variantLabels)
    individualsNum<-length(phenoIID)
    I<-matrix(0,individualsNum,variantNum)
    for( i in 1:length(genoH))
    {
        col=match(genoVariant[i],variantLabels)
        row=match(genoIID[i],phenoIID)
        I[row,col]<-genoH[i]
    }
    return(I)
}
MOJOV.genoVector<-function(x=NULL,y=NULL,
                           codeMethod=c("Proportion","Indicator","ChuanhuaXin"),
                           weightMethod=c("ChuanhuaXin"))
{
  if(is.null(x)|is.null(y))
    stop("You need specify your data.")
  if(class(x)!="matrix")
    stop("The data x should be a matrix.")
  if(length(y)!=dim(x)[1])
    stop("The genotype matrix x is not conformable.")
  codeMethod<-match.arg(codeMethod)
  
  if(codeMethod=="ChuanhuaXin")
  {
      weight<-MOJOV.weight(x=x,y=y,weightMethod=weightMethod)
      X<-apply(x,1,function(x)sum(x*weight))
  }
  if(codeMethod=="Proportion")
  {
      X<-apply(x,1,sum)/(dim(x)[2]*2)
  }
  if(codeMethod=="Indicator")
  {
      X<-apply(x,1,sum)
      X<-unlist(lapply(X,function(x)ifelse(x>0,1,0)))
  }
  return(X)
}
MOJOV.linearRegAnalysis<-function(x=NULL,y=NULL,
                                  testMethod=c("FTest","WaldTest","LRT","Sandwich","all"),...)
{
  if(is.null(x)|is.null(y))
    stop("You need input correct data for x and y.")
  if(dim(as.matrix(x))[2]!=1|dim(as.matrix(y))[2]!=1)
    stop("The vector : X or Y is not one-dimension.")
  if(length(x)!=length(y))
    stop("X and Y is not same length.")
        
  testMethod<-match.arg(testMethod)
  
  if(testMethod=="FTest")
  {
      if(length(unique(x))==1)
          return(1)
      return(summary(lm(y~x))$coefficient[2,4])
  }
  if(testMethod=="WaldTest")
  {
      if(length(unique(x))==1)
          return(1)
      return(MOJOV.wald.test(x,y,...))
  }
  if(testMethod=="Sandwich")
  {
      if(length(unique(x))==1)
          return(1)
      return(MOJOV.saws(x,y,...))
  }
  if(testMethod=="SurveyWald")
  {
      if(length(unique(x))==1)
          return(1)
      return(MOJOV.regTermTest(x,y,...))
  }
  if(testMethod=="LRT")
  {
      if(length(unique(x))==1)
          return(1)
      return(MOJOV.regTermTest(x,y,method="LRT",...))
  }
  if(testMethod=="all")
  {
      if(length(unique(x))==1)
          return(rep(1,9))
      statistic.1<-MOJOV.wald.test(x,y,...)
      statistic.2<-MOJOV.saws(x,y,...)
      statistic.3<-MOJOV.regTermTest(x,y)
      statistic.4<-MOJOV.regTermTest(x,y,method="LRT")
      statistic<-unlist(c(wald.test.chi2=statistic.1$statistic[6]$result$chi2[1],
                          wald.test.p=statistic.1$statistic[6]$result$chi2[3],
                          saws.p1=statistic.2$statistic[9][1],
                          saws.p2=statistic.2$statistic[9][2],
                          regTT.wald.F=statistic.3$statistic[3],
                          regTT.wald.P=statistic.3$statistic[7],
                          regTT.LRT=statistic.4$statistic[3],
                          regTT.LRT=statistic.4$statistic[6],
                          pr=summary(statistic.1$fit)[12]$coefficient[2,4]
      ))
      return(statistic)
  }


}
MOJOV.phenotype<-function(x=NULL,auto=FALSE,gender=TRUE,power=1,Terms=1:3)
{
    if(class(x)!="MOJOV")
        stop("Must specify a MOJOV object.")
    if(!is.logical(auto))
        stop("auto should be logical.")
    if(!is.logical(gender))
        stop("gender should be logical.")
    adjustAuto=auto
    x_which<-phenoData(x)
    phenoAge=phenoAge(x)
    phenoHeight=phenoHeight(x)
    phenoWeight=phenoWeight(x)
    phenoGender=phenoGender(x)
    x_x<-matrix(c(phenoAge,phenoHeight,phenoWeight),,3)
    labels=c("Age","Height","Weight")
    
    if(adjustAuto==T)
    {
        if(is.numeric(power))
        {
            if(length(power)<=1) 
                power<-c(seq(-5,5,0.01))
            power<-power[power!=0]
            pvalue_tmp<-numeric()
            power_tmp<-numeric()
            for ( i in power )
            {
                pvalue_tmp<-c(pvalue_tmp,unlist(shapiro.test(x_which^i)[2]))
                power_tmp<-c(power_tmp,i)
            }
            pvalue_max<-pvalue_tmp[pvalue_tmp==max(pvalue_tmp)]
            power_max<-power_tmp[pvalue_tmp==max(pvalue_tmp)]
            x_which<-x_which^power_max
        }
        else
        {
            if(power!=FALSE)
                power<-c(seq(-5,5,0.01))
            if(power==F)
            {
                power_max=1
                pvalue_max=unlist(shapiro.test(x_which)[2])        
            }
        }
        
        fit0<-glm(x_which~x_x)
        p_vector<-summary(fit0)[12]$coefficients[(2:4),4]
        Terms<-p_vector<=0.05
        if(length(Terms[Terms==T]>0))
        {
            fit<-glm(x_which~x_x[,(1:3)[Terms]])
            adjustTerms=names(x_x)[1:3][Terms]
            adjustPvalue=unlist(summary(fit)[12]$coefficients[2:(length(adjustTerms)+1),4])
            x_res<-fit$residuals  
            MOJOVObj<-new("MOJOV",adjustAuto=adjustAuto,adjustData=x_res,
                          adjustGender=FALSE,adjustPower=power_max,
                          adjustPowerPvalue=pvalue_max,
                          adjustTerms=labels[Terms],
                          adjustPvalue=adjustPvalue,
                          genoFile=genoFile(x),genoIID=genoIID(x),
                          genoChr=genoChr(x),genoPosi=genoPosi(x),
                          genoVariant=genoVariant(x),genoH=genoH(x),
                          genoGene=genoGene(x),phenoFile=phenoFile(x),
                          phenoIID=phenoIID(x),phenoGender=phenoGender(x),phenoAge=phenoAge(x),
                          phenoHeight=phenoHeight(x),phenoWeight=phenoWeight(x),
                          phenoLabel=phenoLabel(x),phenoData=phenoData(x),
                          phenoColumn=phenoColumn(x),phenoIndNum=phenoIndNum(x))
            return(MOJOVObj)
        }
        else
        {                          
            MOJOVObj<-new("MOJOV",adjustAuto=adjustAuto,adjustData=x_which,
                          adjustGender=FALSE,adjustPower=power_max,
                          adjustPowerPvalue=pvalue_max,
                          adjustTerms="NA",adjustPvalue=Inf,
                          genoFile=genoFile(x),genoIID=genoIID(x),
                          genoChr=genoChr(x),genoPosi=genoPosi(x),
                          genoVariant=genoVariant(x),genoH=genoH(x),
                          genoGene=genoGene(x),phenoFile=phenoFile(x),
                          phenoIID=phenoIID(x),phenoGender=phenoGender(x),phenoAge=phenoAge(x),
                          phenoHeight=phenoHeight(x),phenoWeight=phenoWeight(x),
                          phenoLabel=phenoLabel(x),phenoData=phenoData(x),
                          phenoColumn=phenoColumn(x),phenoIndNum=phenoIndNum(x))
            return(MOJOVObj)
            
        }
        
        
    }
    else
    {
        if(mode(power)=="numeric"&length(power)==1)
        {
            x_which=x_which^power
        }else{stop("power should be a number.")}
        if(gender==TRUE)
        {
            #For primary sort
            x_ID<-1:length(x_which)      
            Gender_1<-phenoGender==1
            Gender_2<-phenoGender==2
            x_1<-x_which[Gender_1]
            x_ID_1<-x_ID[Gender_1]      
            x_2<-x_which[Gender_2]
            x_ID_2<-x_ID[Gender_2]
            
            if(mode(Terms)=="numeric")
            {
                x_fit_1<-glm(x_1~x_x[Gender_1,Terms])
                x_fit_2<-glm(x_2~x_x[Gender_2,Terms])
                x_res_1<-x_fit_1$residuals
                x_res_2<-x_fit_2$residuals
                adjustTerms=labels[Terms]
            }
            x_res<-c(x_res_1,x_res_2)[order(c(x_ID_1,x_ID_2))]
            fit<-list(x_fit_1,x_fit_2)
            adjustPvalue=c(unlist(summary(x_fit_1)[12]$coefficients[2:(length(adjustTerms)+1),4]),
                           unlist(summary(x_fit_2)[12]$coefficients[2:(length(adjustTerms)+1),4]))
            names(adjustPvalue)<-c(paste(adjustTerms,"Gender1",sep=""),paste(adjustTerms,"Gender2",sep=""))
        }
        if(gender==FALSE)
        {
            fit<-glm(x_which~x_x[,Terms])
            adjustPvalue=unlist(summary(fit)[12][2:(length(adjustTerms)+1),4])
            x_res<-fit$residuals
        }
        if(power==F)power=1
        adjustPowerPvalue=unlist(shapiro.test(x_which)[2])
        MOJOVObj<-new("MOJOV",adjustAuto=adjustAuto,adjustData=x_res,
                      adjustGender=gender,adjustPower=power,
                      adjustPowerPvalue=adjustPowerPvalue,
                      adjustTerms=labels[Terms],
                      adjustPvalue=adjustPvalue,
                      genoFile=genoFile(x),genoIID=genoIID(x),
                      genoChr=genoChr(x),genoPosi=genoPosi(x),
                      genoVariant=genoVariant(x),genoH=genoH(x),
                      genoGene=genoGene(x),phenoFile=phenoFile(x),
                      phenoIID=phenoIID(x),phenoGender=phenoGender(x),phenoAge=phenoAge(x),
                      phenoHeight=phenoHeight(x),phenoWeight=phenoWeight(x),
                      phenoLabel=phenoLabel(x),phenoData=phenoData(x),
                      phenoColumn=phenoColumn(x),phenoIndNum=phenoIndNum(x))
        return(MOJOVObj)
    }
    
    x_x<-as.matrix(x[,c(2,4,5)])
    if(auto==TRUE)
    {
        if(mode(power)=="numeric")
        {
            if(length(power)<=1) 
                power<-c(seq(-5,5,0.01))
        }
        else
        {
            if(mode(power)=="logical"&length(power)==1)
            {
                if(power!=FALSE)
                {
                    power<-c(seq(-5,5,0.01))
                }
            }
            stop("power should be TRUE or FALSE or numeric vector.")
        }
        if(mode(power)=="numeric")
        {
            power<-power[power!=0]
            pvalue_tmp<-numeric()
            power_tmp<-numeric()
            for ( i in power )
            {
                pvalue_tmp<-c(pvalue_tmp,unlist(shapiro.test(x_which^i)[2]))
                power_tmp<-c(power_tmp,i)
            }
            pvalue_max<-pvalue_tmp[pvalue_tmp==max(pvalue_tmp)]
            power_max<-power_tmp[pvalue_tmp==max(pvalue_tmp)]
            x_which<-x_which^power_max
        }
        
        fit0<-glm(x[,6]~x_x)
        p_vector<-summary(fit0)[12]$coefficients[(2:4),4]
        Terms<-p_vector<=0.05
        fit<-glm(x_which~x_x[,(1:3)[Terms]])
        x_res<-fit$residuals
        phenotype<-list(phenotype=x_res,fit=fit,fit0=fit0,
                        power=data.frame(power=power_max,pvalue=pvalue_max))
        class(phenotype)<-"MOJOV.ad.Phenotype"
        return(phenotype)
    }else{
        phenotype<-list(phenotype=x_res,fit=fit)
        class(phenotype)<-"MOJOV.ad.Phenotype"
        return(phenotype)
    }
    
}


MOJOV.regTermTest<-function(x,y,...)
{
  fit<-glm(y~x)
  statistic<-regTermTest(fit,"x",...)
  return(list(fit=fit,statistic=statistic))
}
MOJOV.saws<-function(x,y,...)
{
  out<-lmfitSaws(model.matrix(~x),y)
  statistic<-saws(out,...)
  fit<-glm(y~x)
  return(list(fit=fit,statistic=statistic))
}
#### MOJOV class definition ####
setClass(Class="MOJOV",
         representation = representation(
             genoFile="character",
             genoIID="character",
             genoChr="character",
             genoPosi="numeric",
             genoVariant="character",
             genoH="numeric",
             genoGene="character",
             phenoFile="character",
             phenoIID="character",
             phenoGender="numeric",
             phenoAge="numeric",
             phenoHeight="numeric",
             phenoWeight="numeric",
             phenoLabel="character",
             phenoData="numeric",
             phenoColumn="numeric",
             phenoIndNum="numeric",
             adjustAuto="logical",
             adjustData="numeric",
             adjustGender="logical",
             adjustPower="numeric",
             adjustPowerPvalue="numeric",
             adjustTerms="character",
             adjustPvalue="numeric",
             varMAF="numeric",
             varFreq="numeric",
             varTot="character",
             varRare="character",
             regionFile="character",
             regionType="character",
             regionChr="character",
             regionPStart="numeric",
             regionPStop="numeric",
             regionGene="character",
             analyCode="character",
             analyWeighted="character",
             resultMethod="character",
             resultPvalue="numeric"
             ),
         prototype=prototype()
)

#### Methods for extracting slots ###
# genotype file
setGeneric('genoFile', function(object) standardGeneric('genoFile'))
setMethod('genoFile','MOJOV', function(object) object@genoFile)

# genotype individuals ID
setGeneric('genoIID', function(object) standardGeneric('genoIID'))
setMethod('genoIID','MOJOV', function(object) object@genoIID)

# genotype chromosome
setGeneric('genoChr', function(object) standardGeneric('genoChr'))
setMethod('genoChr','MOJOV', function(object) object@genoChr)

# genotype position
setGeneric('genoPosi', function(object) standardGeneric('genoPosi'))
setMethod('genoPosi','MOJOV', function(object) object@genoPosi)

# genotype variants labels
setGeneric('genoVariant', function(object) standardGeneric('genoVariant'))
setMethod('genoVariant','MOJOV', function(object) object@genoVariant)

# genotype status
setGeneric('genoH', function(object) standardGeneric('genoH'))
setMethod('genoH','MOJOV', function(object) object@genoH)

# genotype gene labels
setGeneric('genoGene', function(object) standardGeneric('genoGene'))
setMethod('genoGene','MOJOV', function(object) object@genoGene)

# phenotype file
setGeneric('phenoFile', function(object) standardGeneric('phenoFile'))
setMethod('phenoFile','MOJOV', function(object) object@phenoFile)

# phenotype individuals ID
setGeneric('phenoIID', function(object) standardGeneric('phenoIID'))
setMethod('phenoIID','MOJOV', function(object) object@phenoIID)

# phenotype gender
setGeneric('phenoGender', function(object) standardGeneric('phenoGender'))
setMethod('phenoGender','MOJOV', function(object) object@phenoGender)

# phenotype age
setGeneric('phenoAge', function(object) standardGeneric('phenoAge'))
setMethod('phenoAge','MOJOV', function(object) object@phenoAge)

# phenotype height
setGeneric('phenoHeight', function(object) standardGeneric('phenoHeight'))
setMethod('phenoHeight','MOJOV', function(object) object@phenoHeight)

# phenotype weight
setGeneric('phenoWeight', function(object) standardGeneric('phenoWeight'))
setMethod('phenoWeight','MOJOV', function(object) object@phenoWeight)

# phenotype type label
setGeneric('phenoLabel', function(object) standardGeneric('phenoLabel'))
setMethod('phenoLabel','MOJOV', function(object) object@phenoLabel)

# phenotype data
setGeneric('phenoData', function(object) standardGeneric('phenoData'))
setMethod('phenoData','MOJOV', function(object) object@phenoData)

# which column of source phenotype file 
setGeneric('phenoColumn', function(object) standardGeneric('phenoColumn'))
setMethod('phenoColumn','MOJOV', function(object) object@phenoColumn)

# the number of individuals from phenotype file
setGeneric('phenoIndNum', function(object) standardGeneric('phenoIndNum'))
setMethod('phenoIndNum','MOJOV', function(object) object@phenoIndNum)

# if auto choose the best parameters for adjusting phenotype
setGeneric('adjustAuto', function(object) standardGeneric('adjustAuto'))
setMethod('adjustAuto','MOJOV', function(object) object@adjustAuto)

# the phenotype data after adjusted
setGeneric('adjustData', function(object) standardGeneric('adjustData'))
setMethod('adjustData','MOJOV', function(object) object@adjustData)

# if adjusting the phenotype based different gender
setGeneric('adjustGender', function(object) standardGeneric('adjustGender'))
setMethod('adjustGender','MOJOV', function(object) object@adjustGender)

# the power for adjusting
setGeneric('adjustPower', function(object) standardGeneric('adjustPower'))
setMethod('adjustPower','MOJOV', function(object) object@adjustPower)

# the pvalue from shapiro test for the power
setGeneric('adjustPowerPvalue', function(object) standardGeneric('adjustPowerPvalue'))
setMethod('adjustPowerPvalue','MOJOV', function(object) object@adjustPowerPvalue)

# the terms for adjusting
setGeneric('adjustTerms', function(object) standardGeneric('adjustTerms'))
setMethod('adjustTerms','MOJOV', function(object) object@adjustTerms)

# the pvalue from t test for linear regression
setGeneric('adjustPvalue', function(object) standardGeneric('adjustPvalue'))
setMethod('adjustPvalue','MOJOV', function(object) object@adjustPvalue)

# the file name of region of interest
setGeneric('regionFile', function(object) standardGeneric('regionFile'))
setMethod('regionFile','MOJOV', function(object) object@regionFile)

# the type of region position or gene
setGeneric('regionType', function(object) standardGeneric('regionType'))
setMethod('regionType','MOJOV', function(object) object@regionType)

# the chromosome of ROI
setGeneric('regionChr', function(object) standardGeneric('regionChr'))
setMethod('regionChr','MOJOV', function(object) object@regionChr)

# the position start of ROI
setGeneric('regionPStart', function(object) standardGeneric('regionPStart'))
setMethod('regionPStart','MOJOV', function(object) object@regionPStart)

# the position stop of ROI
setGeneric('regionPStop', function(object) standardGeneric('regionPStop'))
setMethod('regionPStop','MOJOV', function(object) object@regionPStop)

# the labels of genes for ROI
setGeneric('regionGene', function(object) standardGeneric('regionGene'))
setMethod('regionGene','MOJOV', function(object) object@regionGene)

# the model for analysis
setGeneric('analyCode', function(object) standardGeneric('analyCode'))
setMethod('analyCode','MOJOV', function(object) object@analyCode)

# the method of weighted for model
setGeneric('analyWeighted', function(object) standardGeneric('analyWeighted'))
setMethod('analyWeighted','MOJOV', function(object) object@analyWeighted)

# the significance test for linear regression
setGeneric('resultMethod', function(object) standardGeneric('resultMethod'))
setMethod('resultMethod','MOJOV', function(object) object@resultMethod)

# the significance P value for linear regression
setGeneric('resultPvalue', function(object) standardGeneric('resultPvalue'))
setMethod('resultPvalue','MOJOV', function(object) object@resultPvalue)

### show ###
setMethod('show', signature='MOJOV', definition=function(object) 
    {    
    symbols1=paste(c(rep("#",40),"\n"),collapse="") 
    symbols2=paste(c(rep("=",40),"\n"),collapse="") 
    region<-1:(min(3,length(object@phenoIID)))
    cat(symbols1)
    cat(" MOJOV is a R package for rare variant\n")
    cat("  Author:Conda E-mail:wkh983@sina.com\n")
    cat(symbols2)
    cat("Genotype Infomation\n")
    cat(symbols2)
    cat("File:")
    cat(object@genoFile)
    cat("\nIndividuals:")
    cat(object@genoIID[region])
    cat("... ...\nChromosome(s):")
    region<-1:(min(3,length(unique(object@genoChr))))
    cat(unique(object@genoChr)[region])
    cat("\nGene Labels:")
    region<-1:(min(3,length(unique(object@genoGene))))
    cat(unique(object@genoGene)[region])
    cat("\nVariant Labels:")
    region<-1:(min(2,length(unique(object@genoVariant))))
    cat(unique(object@genoVariant)[region])
    cat("... ...\nPositions:")
    region<-1:(min(3,length(object@genoPosi)))
    cat(unique(object@genoPosi)[region])
    cat(".\nGenotype Status:")
    region<-1:(min(10,length(object@genoH)))
    cat(object@genoH[region])
    cat("... \n\n")
    cat(symbols2)
    
    cat("Phenotype Infomation\n")
    cat(symbols2)
    cat("File:")
    cat(object@phenoFile)
    cat("\nColumn:")
    cat(object@phenoColumn)
    cat("\nPhenotype Label:")
    cat(object@phenoLabel)
    cat("\n\n")
    cat(symbols2)
    
    cat("Correction Infomation\n")
    if(length(object@adjustAuto)!=0)
    {
        if(object@adjustAuto==T)
            cat("Auto correction\n")
        else
            cat("Non-Auto correction\n")
        cat("Gender:")
        cat(object@adjustGender)
        cat(paste("\nExponent:",object@adjustPower,"  ",object@adjustPowerPvalue,sep=""))
        cat("\nTerms:\n")
        cat(object@adjustTerms)
        cat("\n")
        cat(object@adjustPvalue)
    }
    cat("\n\n")
    cat(symbols2)
    
    cat("Region of interest Infomation\n")
    cat(symbols2)
    if(length(object@regionType)>1)
    {
        cat("Type:")
        cat(object@regionType[1:min(3,length(object@regionType))])
        cat("...")
    }
    if(length(object@regionType)==1)
    {
        cat("Type:")
        cat(object@regionType)
        if(object@regionType=="Single")
        {
            cat("\nGene:")
            cat(object@regionGene)
        }
        if(object@regionType=="Gene")
        {
            cat("File:")
            cat(object@regionFile)
            cat("\nGene:")
            region<-1:(min(3,length(object@regionGene)))
            cat(object@regionGene[region])
        }
        if(object@regionType=="Position")
        {
            cat("File:")
            cat(object@regionFile)
            cat("\nChromosome:")
            region<-1:(min(3,length(object@regionChr)))
            cat(object@regionChr[region])
            cat("\nPosition Start:")
            cat(object@regionPStart[region])
            cat("\nPosition Stop:")
            cat(object@regionPStopt[region])
            cat("\nGene:")
            region<-1:(min(3,length(object@regionGene)))
            cat(object@regionGene[region])
        }
    }
    cat("\n\n")
    cat(symbols2)
    
    cat("Result after analysis\n")
    cat(symbols2)
    cat("Model:")
    cat(object@analyCode)
    cat("\nWeight Method:")
    cat(object@analyWeighted)
    cat("\nSignificance Test Methods:")
    cat(object@resultMethod)
    cat("\nSignificance Test Result:\n")
    region<-1:(min(3,length(object@resultPvalue)))
    cat(object@resultPvalue[region])
    cat("...\n")
    cat(symbols1)    
}
)

MOJOV.simulation<-function(cohortSize=500,nReps=2,theta=10,sites=NULL,affectNum=NULL,
                           MAF=0.01,totalMAF=0.05,lambda=1,sd=NULL,type = c("alpha","belta"),
                           sampleNum=100,outFile=NULL,plot=FALSE,
                           codeMethod=c("Proportion","Indicator","ChuanhuaXin"),
                           weightMethod=c("ChuanhuaXin"),
                           testMethod=c("FTest","WaldTest","LRT","Sandwich","all"),
                           save=NULL)
{
    type<-match.arg(type)
    if(!is.null(sites))
        theta<-paste(theta,"-s",sites)
    if(is.null(outFile))
        outFile<-"MOJOV.simulation.out"
    if(sampleNum>cohortSize)
        stop("The sample number should not greater than population.")
    cohortSize<-cohortSize*2
  commandCharacter<-paste("ms",cohortSize,nReps,"-t",theta,">",outFile)
  cat("ms is running...\nIt may take lots of time, please wait for it patiently...\n")
  system(commandCharacter)
  inFile<-file(outFile,"r")
  systemInfo<-readLines(inFile,n=2)
  pb<-txtProgressBar(min=0,max=nReps,style=3)
  
    #data collection
    result<-numeric()
    resultSummary<-numeric()
  
  for ( i in 1:nReps)
  {
    garbage<-readLines(inFile,n=2)
    if(!is.null(sites))
        garbage<-readLines(inFile,n=1)
    segsites<-as.numeric(unlist(strsplit(readLines(inFile,n=1),split=" "))[2])
    positions<-as.numeric(unlist(strsplit(readLines(inFile,n=1),split=" "))[-1])
    tmp<-readLines(inFile,n=cohortSize)
    tmp<-as.numeric(unlist(strsplit(tmp,split="")))
    x<-matrix(tmp,cohortSize,segsites,byrow=T)
    sitesLabel<-paste("SITE",1:segsites,sep="")
    
    #pairing together to form an individual at random
    sampleIndex<-sample(1:cohortSize,cohortSize)
    sampleIndex<-matrix(sampleIndex,,2)
    sample<-unlist(lapply(1:dim(sampleIndex)[1],function(i)apply(x[sampleIndex[i,],],2,sum)))
    sample<-matrix(sample,,segsites,byrow=T)
    dimnames(sample)<-list(paste("ID",(1:(cohortSize/2)),sep=""),sitesLabel)
    
    #Generating phenotype vector
    freq<-apply(sample,2,sum)/cohortSize
    resultSummary<-rbind(resultSummary,summary.freq(freq))

    if(is.null(sd))
        stdd<-sd(freq)*lambda
    else
        stdd<-sd
    if(type=="alpha")
    {
        index_tmp<-(1:segsites)[freq<=MAF]
        index<-numeric()
        if(is.null(affectNum))
        {
            total<-0
            for(ii in sample(index_tmp,length(index_tmp)))
            {
                total=total+freq[ii]
                index=c(index,ii)
                if(total>totalMAF)break
            }
        }
        else
        {
            if(affectNum<=1)
                affect<-floor(affectNum*length(index_tmp))
            else
                affect<-affectNum
            index<-sample(index_tmp,affect)
        }
        if(length(index)<2)
        {
            warning("The affect sites is less than 2.")
            next
        }
        rareSample<-sample[,index]
        meanPhenotype<-apply(rareSample,1,sum)
        phenotype<-unlist(lapply(meanPhenotype,function(x)rnorm(1,x,stdd)))
    }
    else
        phenotype<-rnorm((cohortSize/2),0,stdd)

    
    #sample at random
    index<-sample((1:(cohortSize/2)),sampleNum)
    x<-sample[index,]
    y<-phenotype[index]
    
    
    result<-rbind(result,unlist(MOJOV.linearRegAnalysis(x=MOJOV.genoVector(x,y,
                                                                           weightMethod=weightMethod,
                                                                           codeMethod=codeMethod),y=y,
                                                        testMethod=testMethod)))
#    return(list(sample=sample,rareSample=rareSample,phenotype=phenotype,
#                x=x,y=y,index=index,result=result))
    
    setTxtProgressBar(pb,i)
  }
    apply(as.matrix(result),2,mean)
    cat("\n")
    if(!is.null(save))
    write.csv(x=result,file=save,row.names=FALSE)
    close(inFile)
    if(is.null(outFile))
        file.remove("MOJOV.simulation.out")
    Summary=apply(as.matrix(resultSummary),2,mean)
    if(plot==TRUE)
        barchart(Summary,col="pink",horizontal=T,xlab="Mean",ylab="Frequency",main="Distribution of Frequency")
    list(TestResult=result,
         SummaryResult=resultSummary,
         Summary=Summary)

}
MOJOV.wald.test<-function(x,y,...)
{
  fit<-glm(y~x)
  statistic<-wald.test(b=coef(fit),Sigma=vcov(fit),Terms=2,...)
  return(list(fit=fit,statistic=statistic))
}
MOJOV.weight<-function(x=NULL,y=NULL,weightMethod="ChuanhuaXin")
{
    weightMethodVector<-c("ChuanhuaXin")
    if(length(weightMethod)!=1&is.character(weightMethod))
        stop("weightMethod should be one character.")
    if(!weightMethod%in%weightMethodVector)
        stop("weightMethod is not proper.")
    sitesNum<-dim(x)[2]
    individualsNum<-dim(x)[1]
    if(weightMethod == "ChuanhuaXin")
    {
        delta<-ifelse(y<mean(y)+sd(y)&y<mean(y)-sd(y),1,0)
        
        weight<-unlist(lapply(1:sitesNum,function(col.index)
        {
            p<-(sum(x[,col.index]*delta)+1)/(2*(sum(delta)+1))
            1/(2*individualsNum*p*(1-p))^(1/2)
        }))
        return(weight)
    }
    if(weightMethod == "")
    {
        weight<-rep(1,individualsNum)
        return(weight)
    }
}

MOJOV.export<-function(x=NULL,file="MOJOV.result.csv",p=0.001)
{
    if(class(x)!="MOJOV")
        stop("Must specify a MOJOV object.")
    regionType<-regionType(x)
    resultPvalue<-resultPvalue(x)
    resultMethod<-resultMethod(x)
    analyCode<-analyCode(x)
    genoFile<-genoFile(x)
    phenoFile<-phenoFile(x)
    phenoLabel<-phenoLabel(x)
    
    if(length(regionType)>1)
    {
        note1<-paste("##Genotype File :",genoFile)
        note2<-paste("##Phenotype File :",phenoFile)
        note3<-paste("##Phenotype Name :",phenoLabel)
        note4<-paste("##Code Method :",analyCode)
        note5<-paste("##Test Method :",resultMethod)
        note6<-paste("##Time Stamp :", Sys.time())
        note<-c("##This file is exported by MOJOV package.",
                note1,note2,note3,note4,note5,note6)

        write.table(x=data.frame(note),row.names=F,col.names=F,
                    file=file,quote=FALSE)
        dataSet<-data.frame(Gene=regionType,pvalue=resultPvalue)
        write.table(x=dataSet,row.names=F,append=T,file=file,sep=",")
    }
    return(x@regionType[x@resultPvalue<=p])
}

summary.freq<-function(x, bin=c(0,0.01,.05,.1,.2,.3,.4,.5))
{
    len<-length(bin)-1
    total<-length(x)
    labels<-character()
    result<-numeric()
    x[x>0.5]<-1-x[x>0.5]
    for(i in 1:len)
    {
        result<-c(result,length(x[x>bin[i]&x<=bin[i+1]]))
        labels<-c(labels,paste("[",bin[i],",",bin[i+1],"]",sep=""))
    }
    result<-result/total
    names(result)<-labels
    result
}

MOJOV.summary<-function(x=NULL,ROI="scan",bin=c(0,0.01,.05,.1,.2,.3,.4,.5),
                        MAF=0.05,...)
{
    if(class(x)!="MOJOV")
        stop("Must specify a MOJOV object.")
    genoIID=genoIID(x)
    genoChr=genoChr(x)
    genoPosi=genoPosi(x)
    genoVariant=genoVariant(x)
    genoH=genoH(x)
    genoGene=genoGene(x)
    phenoFile=phenoFile(x)
    phenoIID=phenoIID(x)
    phenoGender=phenoGender(x)
    phenoAge=phenoAge(x)
    phenoHeight=phenoHeight(x)
    phenoWeight=phenoWeight(x)
    phenoLabel=phenoLabel(x)
    phenoData=phenoData(x)
    phenoColumn=phenoColumn(x)
    phenoIndNum=phenoIndNum(x)
    
    if(ROI=="scan")
    {
        geneList=unique(genoGene)
        pb<-txtProgressBar(min=1,max=length(geneList),style=3)
        pbi=1 
        sumFreq<-numeric()
        geneVarNum<-numeric()
        geneRareNum<-numeric()
        for( i in geneList)
        {
            regionGene=i
            index<-(1:length(genoGene))[genoGene==regionGene]
            tmp<-MOJOV.genoMatrix(genoIID=genoIID[index],genoVariant=genoVariant[index],
                                  genoH=genoH[index],phenoIID=phenoIID,MAF=MAF,...)
            varFreq<-apply(tmp,2,sum)/(phenoIndNum*2)
            sumFreq<-rbind(sumFreq,summary.freq(varFreq,bin=bin))
            geneVarNum<-c(geneVarNum,length(varFreq))
            geneRareNum<-c(geneRareNum,length(varFreq[varFreq<=MAF]))
            setTxtProgressBar(pb,pbi)
            pbi<-pbi+1
        }
        
        bin=names(summary.freq(varFreq,bin=bin))
        
        #total information
        totalVariant<-sum(geneVarNum)
        totalRareVar<-sum(geneRareNum)
        totalGene<-length(geneList)

        #gene total
        geneRatio<-geneRareNum/geneVarNum
        geneDist<-apply(sumFreq,2,mean)
        
        #gene variants
        geneVarTmpValue<-geneVarNum[order(geneVarNum,decreasing=T)]
        geneVarTmpGene<-geneList[order(geneVarNum,decreasing=T)]
        if(geneVarTmpGene[1]=="NULL")
        {
            geneVarTmpValue<-geneVarTmpValue[-1]
            geneVarTmpGene<-geneVarTmpGene[-1]
        }
        geneVarMaxValue<-geneVarTmpValue[1]
        geneVarMaxGene<-geneVarTmpGene[1]
        geneVarMinValue<-geneVarTmpValue[totalGene]
        geneVarMinNum<-length(geneVarTmpValue[geneVarTmpValue==geneVarMinValue])
        geneVarMeanValue<-mean(geneVarTmpValue)
        
        #gene rare variants
        geneRareTmpValue<-geneRareNum[order(geneRareNum,decreasing=T)]
        geneRareTmpGene<-geneList[order(geneRareNum,decreasing=T)]
        if(geneRareTmpGene[1]=="NULL")
        {
            geneRareTmpValue<-geneRareTmpValue[-1]
            geneRareTmpGene<-geneRareTmpGene[-1]
        }
        geneRareMaxValue<-geneRareTmpValue[1]
        geneRareMaxGene<-geneRareTmpGene[1]
        geneRareMinValue<-geneRareTmpValue[length(geneRareTmpValue)]
        geneRareMinNum<-length(geneRareTmpValue[geneRareTmpValue==geneRareMinValue])
        geneRareMeanValue<-mean(geneRareTmpValue)
        
        #gene ratio
        geneRatioMaxValue<-max(geneRatio)
        geneRatioMinValue<-min(geneRatio)
        geneRatioMaxNum<-length(geneRatio[geneRatio==geneRatioMaxValue])
        geneRatioMinNum<-length(geneRatio[geneRatio==geneRatioMinValue])
        geneRatioMeanValue<-mean(geneRatio)
                
        MOJOV.SummaryObject<-new("MOJOV.Summary",
                                 ROI=ROI,
                                 MAF=MAF,
                                 bin=bin,
                                 
                                 #total information
                                 totalGene=totalGene,
                                 totalVariant=totalVariant,
                                 totalRareVar=totalRareVar,
                                 
                                 #gene total
                                 geneVarNum=geneVarNum,
                                 geneRareNum=geneRareNum,
                                 geneList=geneList,
                                 geneDist=geneDist,
                                 geneRatio=geneRatio,
                                 
                                 #gene variants
                                 geneVarMaxValue=geneVarMaxValue,
                                 geneVarMaxGene=geneVarMaxGene,
                                 geneVarMinValue=geneVarMinValue,
                                 geneVarMinNum=geneVarMinNum,
                                 geneVarMeanValue=geneVarMeanValue,
                                 
                                 #gene rare variants
                                 geneRareMaxValue=geneRareMaxValue,
                                 geneRareMaxGene=geneRareMaxGene,
                                 geneRareMinValue=geneRareMinValue,
                                 geneRareMinNum=geneRareMinNum,
                                 geneRareMeanValue=geneRareMeanValue,
                                 
                                 #gene ratio
                                 geneRatioMaxValue=geneRatioMaxValue,
                                 geneRatioMaxNum=geneRatioMaxNum,
                                 geneRatioMinValue=geneRatioMinValue,
                                 geneRatioMinNum=geneRatioMinNum,
                                 geneRatioMeanValue=geneRatioMeanValue
                                 
        )
        return(MOJOV.SummaryObject)
    }
    else
    {
        if(ROI%in%genoGene)
        {
            index<-(1:length(genoGene))[genoGene==regionGene]
            tmp<-MOJOV.genoMatrix(genoIID=genoIID[index],genoVariant=genoVariant[index],
                                  genoH=genoH[index],phenoIID=phenoIID,MAF=MAF,...)
            varFreq<-apply(tmp,2,sum)/(phenoIndNum*2)
            geneDist<-summary.freq(varFreq,bin=bin)
            geneVarNum<-length(varFreq)
            geneRareNum<-length(varFreq[varFreq<=MAF])
            geneRatio<-geneRareNum/geneVarNum
            bin<-names(geneDist(varFreq,bin=bin))
            
            MOJOV.SummaryObject<-new("MOJOV.Summary",
                                     ROI=ROI,
                                     MAF=MAF,
                                     bin=bin,
                                     
                                     #total information
                                     totalGene=numeric(),
                                     totalVariant=numeric(),
                                     totalRareVar=numeric(),
                                     
                                     #gene total
                                     geneVarNum=geneVarNum,
                                     geneRareNum=geneRareNum,
                                     geneList="character",
                                     geneDist=geneDist,
                                     geneRatio=geneRatio,
                                     
                                     #gene variants
                                     geneVarMaxValue=numeric(),
                                     geneVarMaxGene=numeric(),
                                     geneVarMinValue=numeric(),
                                     geneVarMinNum=numeric(),
                                     geneVarMeanValue=numeric(),
                                     
                                     #gene rare variants
                                     geneRareMaxValue=numeric(),
                                     geneRareMaxGene=numeric(),
                                     geneRareMinValue=numeric(),
                                     geneRareMinNum=numeric(),
                                     geneRareMeanValue=numeric(),
                                     
                                     #gene ratio
                                     geneRatioMaxValue=numeric(),
                                     geneRatioMaxNum=numeric(),
                                     geneRatioMinValue=numeric(),
                                     geneRatioMinNum=numeric(),
                                     geneRatioMeanValue=numeric()
                                     )
            return(MOJOV.SummaryObject)  
        }
        else
            stop("ROI should be a gene symbol from your file or \"scan\".")
    }
}

setClass(Class="MOJOV.Summary",
         representation=representation(
             ROI="character",
             MAF="numeric",
             bin="character",

             #total information
             totalGene="numeric",
             totalVariant="numeric",
             totalRareVar="numeric",
             
             #gene total
             geneVarNum="numeric",
             geneRareNum="numeric",
             geneList="character",
             geneDist="numeric",
             geneRatio="numeric",
             
             #gene variants
             geneVarMaxValue="numeric",
             geneVarMaxGene="character",
             geneVarMinValue="numeric",
             geneVarMinNum="numeric",
             geneVarMeanValue="numeric",
             
             #gene rare variants
             geneRareMaxValue="numeric",
             geneRareMaxGene="character",
             geneRareMinValue="numeric",
             geneRareMinNum="numeric",
             geneRareMeanValue="numeric",
             
             #gene ratio
             geneRatioMaxValue="numeric",
             geneRatioMaxNum="numeric",
             geneRatioMinValue="numeric",
             geneRatioMinNum="numeric",
             geneRatioMeanValue="numeric"
                           
             ),
         prototype=prototype()
         )

setMethod('show', signature="MOJOV.Summary", 
          definition=function(object)
              {
              if(object@ROI=="scan")
              {
                  symbols1=paste(c(rep("#",40),"\n"),collapse="")
                  cat(symbols1)
                  cat("Summary statistic infomations")
                  cat("\nGene :", object@ROI)
                  cat("\nThreshold value :",object@MAF)
                  cat("\nThe total number of variants :",object@totalVariant)
                  cat("\nThe total number of rare variants :",object@totalRareVar)
                  cat("\nThe total number of gene :",object@totalGene)
                  cat("\nThe average distribution for all genes :\n")
                  for(i in object@bin)
                  {
                      cat(i,":",object@geneDist[object@bin==i],"\n")
                  }
                  
                  cat("\nInformation for one gene.")
                  cat("\nAverage number of variants :",object@geneVarMeanValue)
                  cat("\nMax number of variants :",object@geneVarMaxGene,
                      object@geneVarMinValue)
                  cat("\nMin number of variants :",object@geneVarMinValue,
                      "     COUNTS :",object@geneVarMinNum)
                  
                  cat("\nAverage number of rare variants :",object@geneRareMeanValue)
                  cat("\nMax number of rare variants :",object@geneRareMaxGene,
                      object@geneRareMinValue)
                  cat("\nMin number of rare variants :",object@geneRareMinValue,
                      "     COUNTS :",object@geneRareMinNum)
                  
                  cat("\nAverage ratio for all genes :",object@geneRatioMeanValue)
                  cat("\nMax ratio for all genes :",object@geneRatioMaxValue,
                      "     COUNTS :",object@geneRatioMaxNum)
                  cat("\nMin ratio for all genes :",object@geneRatioMinValue,
                      "     COUNTS :",object@geneRatioMinNum)
                  
              }
              else
              {
                  symbols1=paste(c(rep("#",40),"\n"),collapse="")
                  cat(symbols1)
                  cat("Summary statistic infomations")
                  cat("\nGene :", object@ROI)
                  cat("\nThreshold value :",object@MAF) 
                  cat("\nThe number of variants :",object@geneVarNum)
                  cat("\nThe number of rare variants :",object@geneRareNum)
                  cat("\nThe ratio :",object@geneRatio)
                  cat("\nThe average distribution for all genes :\n")
                  for(i in object@bin)
                  {
                      cat(i,":",object@geneDist[object@bin==i],"\n")
                  }
              }
          })