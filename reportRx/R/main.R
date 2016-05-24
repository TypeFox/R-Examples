#' Plot KM curve
#'
#'This function will plot a KM curve with possible stratification. You can specifyif you want
#'a legend or confidence bands as well as the units of time used.
#'
#' @param data dataframe containing your data
#' @param response character vector with names of columns to use for response
#' @param group string specifiying the column name of stratification variable
#' @param pos what position you want the legend to be. Current option are bottomleft and topright
#' @param units string specifying what the unit of time is use lower case and plural
#' @param CI boolean to specify if you want confidence intervals
#' @param legend boolean to specify if you want a legend
#' @param title title of plot
#' @keywords plot
#' @export 
#' @examples
#' require(survival)
#' lung$sex<-factor(lung$sex)
#' plotkm(lung,c("time","status"))
#' plotkm(lung,c("time","status"),"sex")
plotkm<-function(data,response,group=1,pos="bottomleft",units="months",CI=F,legend=T, title=""){
  if(class(group)=="numeric"){  
    kfit<-survfit(as.formula(paste("Surv(",response[1],",",response[2],")~1",sep="")),data=data)
    sk<-summary(kfit)$table
    levelnames<-paste("N=",sk[1], ", Events=",sk[4]," (",round(sk[4]/sk[1],2)*100,"%)",sep="")
    if(title=="")  title<-paste("KM-Curve for ",nicename(response[2]),sep="")
    
  }else if(length(group)>1){
    return("Currently you can only stratify by 1 variable")
  }else{
    if(class(data[,group])!="factor")
      stop("group must be a vactor variable. (Or leave unspecified for no group)")
    lr<-survdiff(as.formula(paste("Surv(",response[1],",",response[2],")~", paste(group,collapse="+"),sep="")),data=data)
    lrpv<-1-pchisq(lr$chisq, length(lr$n)- 1)
    levelnames<-levels(data[,group])
    kfit<-survfit(as.formula(paste("Surv(",response[1],",",response[2],")~", paste(group,collapse="+"),sep="")),data=data)
    if(title=="") title<-paste("KM-Curve for ",nicename(response[2])," stratified by ", nicename(group),sep="")
    levelnames<-sapply(1:length(levelnames), function(x){paste(levelnames[x]," n=",lr$n[x],sep="")})
    
  }  
  
  
  plot(kfit,mark.time=T, lty=1:length(levelnames),xlab=paste("Time (",cap(units),")",sep=""),
       ylab="Suvival Probability ",cex=1.1, conf.int=CI,
       main=title)
  
  
  if(legend){
    if(class(group)=="numeric"){legend(pos,levelnames,lty=1:length(levelnames),bty="n")
    }else{ legend(pos,c(levelnames,paste("p-value=",pvalue(lrpv)," (Log Rank)",sep="")),
                  col=c(rep(1,length(levelnames)),"white"),lty=1:(length(levelnames)+1),bty="n")}
  }
}

#'Get event time summary dataframe
#'
#'This function will output a dataframe with usefull summary statistics from a coxph model
#'
#'@param data dataframe containing data
#'@param response character vector with names of columns to use for response
#'@param group string specifiying the column name of stratification variable
#'@param times numeric vector of times you want survival time provbabilities for.
#'@keywords dataframe
#'@export
#'@examples
#'require(survival)
#'lung$sex<-factor(lung$sex)
#'etsum(lung,c("time","status"),"sex")
#'etsum(lung,c("time","status"))
#'etsum(lung,c("time","status"),"sex",c(1,2,3))
etsum<- function(data,response,group=1,times=c(12,24)){
  if(class(group)=="numeric"){
  kfit<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",group,sep=""))  ,data=data))
  maxtime=max(kfit$time)
  times[times>maxtime]=maxtime
  kfit2<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",group,sep="")) ,data=data),times=times)
  tab<-as.data.frame(cbind(strata=as.character(kfit2$strata),times=kfit2$time,SR=paste(round(kfit2$surv*100,0)," (",round(kfit2$lower*100,0),"-",round(kfit2$upper*100,0),")",sep="")))
  tbl<-kfit2$table
  }else{
    if(class(data[,group])!="factor")
      stop("group variable must be factor or leave unspecified for no group")
    tab<-lapply(levels(data[,group]),function(level){
      subdata<-subset(data,data[,group]==level)
      kfit<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",1,sep=""))  ,data=subdata))
      maxtime=max(kfit$time)
      times[times>maxtime]=maxtime
      kfit2<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",1,sep="")) ,data=subdata),times=times)
      list(cbind(strata=paste0(group,"=",level),times=kfit2$time,SR=paste(round(kfit2$surv*100,0)," (",round(kfit2$lower*100,0),"-",round(kfit2$upper*100,0),")",sep="")),kfit2$table)})
    tbl=t(sapply(tab,"[[",2))
    rownames(tbl)=sapply(levels(data[,group]),function(level)paste0(group,"=",level))    
    tab=do.call(rbind.data.frame,lapply(tab,"[[",1))
}
  
  if(class(group)!="numeric"){
    kfit<-summary(survfit(as.formula(paste("Surv(",response[1],",",response[2],")~",group,sep=""))  ,data=data))
    med=by(data,data[,group],function(x) median(x[,response[1]],na.rm=T))
    min=by(data,data[,group],function(x) min(x[,response[1]],na.rm=T))
    max=by(data,data[,group],function(x) max(x[,response[1]],na.rm=T))
    survtimes<-data.frame(strata=as.character(kfit$strata),kfit$time)   
    minst<-round(as.numeric(by(survtimes,survtimes$strata,function(x) min (x[,2]))),1)
    maxst<-round(as.numeric(by(survtimes,survtimes$strata,function(x) max (x[,2]))),1)    
    tab<-cast(tab, strata ~ times)
    names<-names(tab)
    tab<-data.frame(tab)
    names(tab)<-names
    tab[,1]<-levels(data[,group])
    if(length(times)>1){
      indx<-c(0,sapply(sort(as.numeric(names(tab)[-1])),function(x){which(as.numeric(names(tab)[-1])==x)}))+1
      tab<-tab[,indx]
      tab<-tab[c(2:length(tab),1)]
    }else{
      tab<-tab[c(2:length(tab),1)]
    }
    noeventsindx<-ifelse(length(which(tbl[,4]==0))!=0,
                         which(tbl[,4]==0),NA)
    if(!is.na(noeventsindx)){
      for(i in noeventsindx){
        if(i==1){
          minst<-c(0,minst)
          maxst<-c(0,maxst)
        }else if(i>length(minst)){
          minst<-c(minst,0)
          maxst<-c(maxst,0)
        }else{
          minst<-c(minst[1:i-1],0,minst[i:length(minst)])
          maxst<-c(maxst[1:i-1],0,maxst[i:length(maxst)])
        }}}     
    
    
    tab<-cbind("n"=tbl[,1],"Events"=tbl[,4], "MedKM"=round(tbl[,5],1),
               "LCI"=round(tbl[,6],1), "UCI"=round(tbl[,7],1),
               "MedFU"=round(as.numeric(med),1),
               "MinFU"=round(as.numeric(min),1),"MaxFU"=round(as.numeric(max),1),
               "MinET"=minst,"MaxET"=maxst,tab)
    rownames(tab)<-NULL
  }else{
    med=median(data[,response[1]],na.rm=T)
    min=min(data[,response[1]],na.rm=T)
    max=max(data[,response[1]],na.rm=T)
    if(length(times)>1){
      tab<-data.frame(t(tab))
      rownames(tab)<-NULL      
      names(tab)<-as.numeric(as.matrix(tab[1,]))
      tab<-tab[-1,]      
    }else{
      rownames(tab)<-NULL
      names(tab)[2]<-times
      tab<-tab[-1]
    } 
    tab<-cbind("n"=tbl[1],"Events"=tbl[4],"MedKM"=round(tbl[5],1),"LCI"=round(tbl[6],1), "UCI"=round(tbl[7],1),
               "MedFU"=round(as.numeric(med),1),"MinFU"=round(as.numeric(min),1),"MaxFU"=round(as.numeric(max),1),
               "MinET"=round(min(kfit$time),1),"MaxET"=round(max(kfit$time),1),tab)
    rownames(tab)<-NULL
  }
  return(tab)
}

#'Print LaTeX event time summary
#'
#'Wrapper for the etsum function that prints paragraphs of text in LaTeX
#'
#'@param data dataframe containing data
#'@param response character vector with names of columns to use for response
#'@param group string specifiying the column name of stratification variable
#'@param times numeric vector of times you want survival time provbabilities for.
#'@param units string indicating the unit of time. Use lower case and plural.
#'@keywords print
#'@export 
#'@examples
#'require(survival)
#'lung$sex<-factor(lung$sex)
#'petsum(lung,c("time","status"),"sex")
#'petsum(lung,c("time","status"))
#'petsum(lung,c("time","status"),"sex",c(1,2,3),"months")
petsum<-function(data,response,group=1,times=c(12,14),units="months"){
  t<-etsum(data,response,group,times)
  
  #plotkm(nona,response,group)
  
  names<-names(t)
  if("strata"%in% names){
    strta<-sapply(t[,"strata"], function(x) paste(x,": ",sep=""))
    offset<-2
    ofst<-1
  }else{
    strta=matrix(c("",""))
    offset<-1
    ofst<-0
  }
  
  
  out<-sapply(seq_len(nrow(t)),function(i){
    
    if(is.na(t[i,3])) {km<-paste("The KM median event time has not been achieved due to lack of events.",sep="")
    }else if (!is.na(t[i,5])){km<-paste("The KM median event time is ",t[i,3]," with 95",sanitizestr("%")," confidence Interval (",t[i,4],",",t[i,5],").",sep="")
    }else{km<-paste("The KM median event time is ",t[i,3], " ",units, " with 95",sanitizestr("%")," confidence Interval (",t[i,4],",",t[i,10],").",sep="")}
    
    # if at least one event
    if(t[i,2]!=0){
      flet<-paste(" The first and last event times occurred at ",t[i,9],
                  " and ",t[i,10]," ",units," respectively. ",sep="")
      
      psindex=11:(ncol(t)-ofst)
      psindex=psindex[which(!is.na(t[i,psindex]))]
      if(length(psindex)>1){
        lastindex=psindex[length(psindex)]
        firstindex=psindex[-length(psindex)]
        ps<-paste("The ",paste(names[firstindex], collapse=", ")," and ",names[lastindex], " " , substring(units,1,nchar(units)-1),
                  " probabilities of 'survival' and their 95",sanitizestr("%")," confidence intervals are ",
                  paste(sapply(t[i,firstindex],function(x) paste(x)),collapse=", ")," and ", t[i,lastindex], " percent.",sep="")
        
      }else{
        ps<-paste("The ",names[psindex]," ", substring(units,1,nchar(units)-1),
                  " probability of 'survival' and 95",sanitizestr("%")," confidence interval is ",
                  t[i,psindex]," percent.",sep="")  
      }
      #if no events
    }else{
      km=""
      ps=""
      flet=""
    }
    
    
    out<-paste(lbld(sanitizestr(nicename(strta[i])))," There are ",t[i,1]," patients. There were ",t[i,2],
               " (",round(100*t[i,2]/t[i,1],0),sanitizestr("%"),") events. The median and range of the follow-up times is ",
               t[i,6]," (",t[i,7],"-",t[i,8],") ",units,". ", km, flet,ps,sep="")
    cat("\n",out,"\n")
  })
}

#'Get covariate summary dataframe
#'
#'Returns a dataframe corresponding to a descriptive table
#'
#'@param data dataframe containing data
#'@param covs character vector with the names of columns to include in table
#'@param maincov covariate to stratify table by
#'@param numobs named list overriding the number of people you expect to have the covariate
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@keywords dataframe
#'@export
covsum<-function(data,covs,maincov=NULL,numobs=NULL,markup=T,sanitize=T,nicenames=T){
  
  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity    
  }
  if(!sanitize) sanitizestr<-identity  
  if(!nicenames) nicename<-identity
  if(!is.null(maincov)){
    levels<-names(table(data[,maincov]))
    levels<-c(list(levels),as.list(levels))
  }else{
    levels<-"NOMAINCOVNULLNA"
  }
  N=nrow(data)
  if(!is.null(maincov)){
    nmaincov<-c(sum(table(data[,maincov])),table(data[,maincov]))
  }else{
    nmaincov<-N
    p<-NULL
  }
  out<-lapply(covs,function(cov){
    ismiss=F
    n<-sum(table(data[,cov]))
    
    #Set up the first coulmn
    factornames<-NULL
    if(is.null(numobs[[cov]]))  numobs[[cov]]<-nmaincov
    if(numobs[[cov]][1]-n>0) {ismiss=T
                              factornames<-c(factornames,"Missing")
    }
    #if the covariate is a factor
    if(is.factor(data[,cov])){
      factornames<-c(levels(data[,cov]),factornames)
      if(!is.null(maincov)){
        p<-try(lpvalue(fisher.test(data[,maincov],data[,cov])$p.value))
        if(class(p)=="try-error") p<-chisq.test(data[,maincov],data[,cov])$p.value
        p<-lpvalue(p)
      } 
      
      
      #set up the main columns
      onetbl<-mapply(function(sublevel,N){
        missing<-NULL
        if(sublevel[1]!="NOMAINCOVNULLNA"){
          subdata<-subset(data,subset=data[,maincov]%in%sublevel)
        }else{
          subdata<-data
        }
        table<-table(subdata[,cov])
        tbl<-table(subdata[,cov])
        n<-sum(tbl)
        prop<-round(tbl/n,2)*100
        prop<-sapply(prop,function(x){if(!is.nan(x)){x} else{0}})
        tbl<-mapply(function(num,prop){paste(num," (",prop,")",sep="")},tbl,prop)       
        if(ismiss) missing<-N-n 
        tbl<-c(tbl,lbld(missing))       
        return(tbl)
      },levels,numobs[[cov]])
      
      #if the covariate is not a factor
    }else{
      #setup the first column
      factornames<-c("Mean (sd)", "Median (Min,Max)",factornames)
      if(!is.null(maincov)){
        p<-try(anova(lm(data[,cov]~data[,maincov]))[5][[1]][1])
        if(class(p)=="try-error") p<-NA
        p<-lpvalue(p)}
      
      
      #set up the main columns
      onetbl<-mapply(function(sublevel,N){
        missing<-NULL
        if(sublevel[1]!="NOMAINCOVNULLNA"){
          subdata<-subset(data,subset=data[,maincov]%in%sublevel)
        }else{subdata<-data}
        summary<-round(summary(subdata[,cov]),1)
        meansd<-paste(summary[4]," (", round(sd(subdata[,cov],na.rm=T),1),")",sep="")
        mmm<-paste(summary[3]," (",summary[1],",",summary[6],")",sep="")
        
        
        #if there is a missing in the whole data
        if(ismiss){          
          n<-sum(table(subdata[,cov]))
          missing<-N-n 
        }
        tbl<-c(meansd,mmm,lbld(missing))
        
        return(tbl)}   
                     ,levels,numobs[[cov]])}
    
    #Add the first column to the main columns and get the matrix ready for later
    factornames<-addspace(sanitizestr(nicename(factornames)))    
    onetbl<-cbind(factornames,onetbl)
    
    if(!is.null(maincov)){
      onetbl<-rbind(c(lbld(sanitizestr(nicename(cov))),rep("",length(levels[[1]])+1)),onetbl)
      onetbl<-cbind(onetbl,c(p,rep("",nrow(onetbl)-1)))
    }else{
      onetbl<-rbind(c(lbld(sanitizestr(nicename(cov))),""),onetbl)
    }
    rownames(onetbl)<-NULL
    colnames(onetbl)<-NULL
    return(onetbl)})
  table<-do.call("rbind", lapply(out, data.frame, stringsAsFactors = FALSE))
  rownames(table)<-NULL
  if(!is.null(maincov)){
    colnames(table)<-c("Covariate",paste("Full Sample (n=",N,")",sep=""),
                       mapply(function(x,y){paste(x," (n=",y,")",sep="")},
                              names(table(data[,maincov])),table(data[,maincov])),"p-value")
  }else{
    colnames(table)<-c("Covariate",paste("n=",N,sep=""))
    
  }
  colnames(table)<-sanitizestr(colnames(table))
  return(table)
}

#'Print covariate summary Latex
#'
#'Returns a dataframe corresponding to a descriptive table
#'
#'@param data dataframe containing data
#'@param covs character vector with the names of columns to include in table
#'@param maincov covariate to stratify table by
#'@param numobs named list overriding the number of people you expect to have the covariate
#'@param TeX boolean indicating if you want to be able to view extra long tables in the LaTeX pdf. If TeX is T then the table will not convert properly to docx
#'@keywords print
#'@export
pcovsum<-function(data,covs,maincov=NULL,numobs=NULL,TeX=F){
  if(!TeX){
    print.xtable(xtable(covsum(data,covs,maincov,numobs)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
  }else{  
    
    print.xtable(xtable(covsum(data,covs,maincov,numobs)),include.rownames=F,sanitize.text.function=identity,table.placement="H",floating=FALSE,tabular.environment="longtable")
  }}

#'Get univariate summary dataframe
#'
#'Returns a dataframe corresponding to a univariate table
#'
#'@param response string vector with name of response 
#'@param covs character vector with the names of columns to fit univariate models to
#'@param data dataframe containing data
#'@param type string indicating he type of univariate model to fit. The function will try and guess what type you want based on your response. If you want to override this you can manually specify the type.
#'Options in clude "linear", "logistic", "coxph", "crr", "boxcox","logistic"
#'@param strata character vector of covariates to stratify by. Only used for coxph and crr
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@param testing boolean to indicate if you want to print out the covariates before the model fits.
#'This will allow you to see which model is not fitting if the function throws an error
#'@keywords dataframe
#'@export
uvsum<-function(response,covs,data,type=NULL,strata=1,markup=T,sanitize=T,nicenames=T,testing=F){
  
  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity    
  }
  if(!sanitize) sanitizestr<-identity  
  if(!nicenames) nicename<-identity
  
  if(class(strata)!="numeric") {strata<-sapply(strata,function(stra){paste("strata(",stra,")",sep="")})
  }else{strata<-""}
  
  if(!is.null(type)){
    if(type=="logistic"){beta<-"OR(95%CI)"
    }else if(type=="linear"|type=="boxcox"){beta<-"Estimate(95%CI)"
    }else if(type=="coxph"|type=="crr"){beta<-"HR(95%CI)"
    }else{stop("type must be either coxph, logisitc, linear, coxbox, crr (or NULL)")
    }}else
    {if(length(response)==2) {
      if(length(unique(data[,response[2]]))<3){type<-"coxph"
      }else{type<-"crr"}
      beta<-"HR(95%CI)"        
    }else if (length(unique(data[,response] ))==2) {type<-"logistic"
                                            beta<-"OR(95%CI)"                                            
    }else {type<-"linear"
           beta<-"Estimate(95%CI)"
    }
    }
  if(strata!="" & type!="coxph") stop("strata can only be used with coxph")
  
  
  
  out<-lapply(covs,function(cov){
    cov2<-cov
    if(testing) print(cov)
    if(is.factor(data[,cov])){
      levelnames<-sapply(sapply(sapply(levels(factor(data[,cov])),nicename),sanitizestr),addspace)
      cov<-lbld(sanitizestr(nicename(cov)))
      title<-NULL
      body<-NULL
      if(type=="coxph"){
        m2<-coxph(as.formula(paste(paste("Surv(",response[1],",",response[2],")",sep=""),"~",cov2,ifelse(strata=="","","+"),paste(strata,collapse="+"),sep="")),data=data)        
        hazardratio<-c("Reference",apply(matrix(summary(m2)$conf.int[,c(1,3,4)],ncol=3),1,psthr))    
        pvalue<-c("",sapply(summary(m2)$coef[,5],lpvalue))
        title<-c(cov,"","",lpvalue(summary(m2)$waldtest[3]))
      }else if(type=="crr"){        
        m2<-crrRx(as.formula(paste(paste(response,collapse="+"),"~",cov2,sep="")),data=data)        
        hazardratio<-c("Reference",apply(matrix(summary(m2)$conf.int[,c(1,3,4)],ncol=3),1,psthr))    
        pvalue<-c("",sapply(summary(m2)$coef[,5],lpvalue))        
        globalpvalue<-try(wald.test(b=m2$coef,Sigma=m2$var,Terms=seq_len(length(m2$coef)))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        title<-c(cov,"","",lpvalue(globalpvalue))
        
      }else if(type=="logistic"){
        m2<-glm(as.formula(paste(response,"~",cov2,sep="")),family="binomial",data=data)   
        #globalpvalue<-1-pchisq(2*(summary(m2)$null.deviance-summary(m2)$deviance),summary(m2)$df.null-summary(m2)$df.residual)
        globalpvalue<-try(wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        
        m<-summary(m2)$coefficients
        hazardratio<-c("Reference",apply(cbind(exp(m[-1,1]),exp(m[-1,1]-1.96*m[-1,2]),exp(m[-1,1]+1.96*m[-1,2])),1,psthr))        
        pvalue<-c("",sapply(m[-1,4],lpvalue))
        title<-c(cov,"","",lpvalue(globalpvalue))
        
      }else if(type=="linear"|type=="boxcox"){
        
        if(type=="linear"){m2<-lm(as.formula(paste(response,"~",cov2,sep="")),data=data)
        }else{m2<-boxcoxfitRx(as.formula(paste(response,"~",cov2,sep="")),data=data)}
        m<-summary(m2)$coefficients
        #globalpvalue<-anova(m2)[5][[1]][1])
        globalpvalue<-try(wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        
        
        hazardratio<-c("Reference",apply(cbind(m[-1,1],m[-1,1]-1.96*m[-1,2],m[-1,1]+1.96*m[-1,2]),1,psthr))        
        pvalue<-c("",sapply(m[-1,4],lpvalue))
        title<-c(cov,"","",lpvalue(globalpvalue))            
      }
      if(length(levelnames)==2){
        body<-cbind(levelnames,hazardratio,c("",""),c("",""))    
      }else{
        body<-cbind(levelnames,hazardratio,pvalue,rep("",length(levelnames)))      
      }
      out<-rbind(title,body)
      rownames(out)<-NULL
      colnames(out)<-NULL
      return(list(out,nrow(out)))
    }else
    {
      cov<-lbld(sanitizestr(nicename(cov)))
      if(type=="coxph"){
        m2<-coxph(as.formula(paste(paste("Surv(",response[1],",",response[2],")",sep=""),"~",cov2,ifelse(strata=="","","+"),paste(strata,collapse="+"),sep="")),data=data)        
        out<-matrix(c(cov,psthr(summary(m2)$conf.int[,c(1,3,4)]),"",lpvalue(summary(m2)$waldtest[3])),ncol=4)
        
      }else if(type=="crr"){
        
        m2<-crrRx(as.formula(paste(paste(response,collapse="+"),"~",cov2,sep="")),data=data)
        globalpvalue<-try(wald.test(b=m2$coef,Sigma=m2$var,Terms=seq_len(length(m2$coef)))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        out<-matrix(c(cov,psthr(summary(m2)$conf.int[,c(1,3,4)]),"",lpvalue(globalpvalue)),ncol=4)
        
      }else if(type=="logistic"){
        m2<-glm(as.formula(paste(response,"~",cov2,sep="")),family="binomial",data=data)    
        
        m<-summary(m2)$coefficients
        
        #globalpvalue<-1-pchisq(2*(summary(m2)$null.deviance-summary(m2)$deviance),summary(m2)$df.null-summary(m2)$df.residual)
        globalpvalue<-try(wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        
        out<-matrix(c(cov,psthr(c(exp(m[-1,1]),exp(m[-1,1]-1.96*m[-1,2]),exp(m[-1,1]+1.96*m[-1,2]))),"",lpvalue(globalpvalue)),ncol=4)
        
        
      }else if(type=="linear"|type=="boxcox"){
        if(type=="linear"){m2<-lm(as.formula(paste(response,"~",cov2,sep="")),data=data)
        }else{m2<-boxcoxfitRx(as.formula(paste(response,"~",cov2,sep="")),data=data)}
        #globalpvalue<-anova(m2)[5][[1]][1])
        globalpvalue<-try(wald.test(b=m2$coefficients[-1],Sigma=vcov(m2)[-1,-1],Terms=seq_len(length(m2$coefficients[-1])))$result$chi2[3]);
        if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
        
        m<-summary(m2)$coefficients
        
        out<-matrix(c(cov,psthr(c(m[-1,1],m[-1,1]-1.96*m[-1,2],m[-1,1]+1.96*m[-1,2])),"",lpvalue(globalpvalue)),ncol=4)
        
      }
      return(list(out,nrow(out)))}})
  table<-lapply(out,function(x){return(x[[1]])})
  table<-do.call("rbind", lapply(table, data.frame, stringsAsFactors = FALSE))
  
  colnames(table)<-sapply(c("Covariate",sanitizestr(beta),"p-value","Global p-value"),lbld)
  return(table)  
}

#'Print univariate summary LaTeX table
#'
#'Returns a LaTeX table of the univariate summary
#'
#'@param response string vector with name of response 
#'@param covs character vector with the names of columns to fit univariate models to
#'@param data dataframe containing data
#'@param type string indicating he type of univariate model to fit. The function will try and guess what type you want based on your response. If you want to override this you can manually specify the type. Options in clude "linear", "logistic", "coxph", "crr", "boxcox","logistic"
#'@param strata character vector of covariates to stratify by. Only used for coxph and crr
#'@param TeX boolean indicating if you want to be able to view extra long tables in the LaTeX pdf. If TeX is T then the table will not convert properly to docx
#'@keywords dataframe
#'@export
puvsum<-function(response,covs,data,type=NULL,strata=1,TeX=F){
  if(!TeX){
    print.xtable(xtable(uvsum(response,covs,data,type,strata)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
  }else{
    print.xtable(xtable(uvsum(response,covs,data,type,strata)),include.rownames=F,sanitize.text.function=identity,table.placement="H",floating=FALSE,tabular.environment="longtable")
  }
  
}

#'Get multivariate summary dataframe
#'
#'Returns a dataframe corresponding to a univariate table
#'
#'@param model fitted model object
#'@param data dataframe containing data
#'@param markup boolean indicating if you want latex markup
#'@param sanitize boolean indicating if you want to sanitize all strings to not break LaTeX
#'@param nicenames booling indicating if you want to replace . and _ in strings with a space
#'@keywords dataframe
#'@export
mvsum<-function(model,data,markup=T,sanitize=T,nicenames=T){
  if(!markup){
    lbld<-identity
    addspace<-identity
    lpvalue<-identity    
  }
  if(!sanitize) sanitizestr<-identity  
  if(!nicenames) nicename<-identity  
  
  call<-paste(deparse(summary(model)$call),collapse="")
  call<-unlist(strsplit(call,"~",fixed=T))[2]
  call<-unlist(strsplit(call,",",fixed=T))[1]
  if(substr(call,nchar(call),nchar(call))=="\"") call<-substr(call,1,nchar(call)-1)
  call<-unlist(strsplit(call,"\"",fixed=T))[1]  
  call<-unlist(strsplit(call,"+",fixed=T))
  call<-unlist(strsplit(call,"*",fixed=T))
  call<-unique(call)
  call<-call[which(is.na(sapply(call,function(cov){charmatch("strata(",cov)}))==T)]
  call<-gsub("\\s","", call)
  type<-class(model)[1]  
  
  if(type!="lme"){
  betanames<-attributes(summary(model)$coef)$dimnames[[1]]
  }else{
    betanames<-names(model$coef$fixed)
  }
  
  
 
    if(type=="glm"){ beta<-"OR(95%CI)"
                          betanames<-betanames[-1]
    }else if(type=="lm"|type=="lm"|type=="lme"){ beta<-"Estimate(95%CI)"
                              betanames<-betanames[-1]
    }else if(type=="coxph"|type=="crr"){beta<-"HR(95%CI)"
    }else{ stop("type must be either coxph, logisitc, lm, crr, lme (or NULL)")}
    
  ucall=unique(call)
 
#   indx<-as.vector(sapply(betanames,function(string){
#     
#     indx<-which(sapply(call,function(cov){charmatch(cov,string)})==1)
#     if(length(indx)==1) return(indx)
#     #If one  facorname is a subset of another
#     indx2<-which.max(sapply(call[indx],nchar))
#     if(length(indx2)==1) return(indx[indx2])
#     indx3<-which(sapply(call[indx2],function(c){substr(betaname,1,nchar(c))==c}))
#     if(length(indx3)==1)  return(call[indx[indx2[indx3]]])  
#     return(-1)
#   }))
  
  indx=matchcovariate(betanames,ucall)
  if(min(indx)==-1) stop("Factor name + level name is the same as another factor name. Please change. Will fix this issue later") 
  
 
  
  
  
  y<-betaindx(indx)
  if(type%in%c("lm","glm","lm","lme")){
    y<-lapply(y,function(x){ x+1})
    betanames<-c("intercept",betanames)
  }
  
  out<-lapply(y,function(covariateindex){
    
    #Get attribute names and split by ineractions
    betaname<-betanames[covariateindex]
    betaname<-strsplit(betaname,":",fixed=T)
    #get the covariate names
    oldcovname<-covnm(betaname[[1]],call)
    
    #get the levelnames
    levelnames<-unlist(lapply(betaname,function(level){
      paste(mapply(function(lvl,cn){result<-unlist(strsplit(lvl,cn,fixed=T))[2]
                                    out<-ifelse(is.na(result),cn,result)},level,oldcovname),collapse=":")}))
    levelnames<-addspace(sanitizestr(nicename(levelnames)))
    covariatename<-lbld(sanitizestr(nicename(paste(oldcovname,collapse=":"))))
    reference=NULL
    title=NULL
    body=NULL
    if(type=="lme"){
      globalpvalue<-try(wald.test(b=model$coef$fixed[covariateindex],Sigma=vcov(model)[covariateindex,covariateindex],Terms=seq_along(covariateindex))$result$chi2[3]);
    }
    else if(type!="crr"){
      globalpvalue<-try(wald.test(b=coef(model)[covariateindex],Sigma=vcov(model)[covariateindex,covariateindex],Terms=seq_along(covariateindex))$result$chi2[3]);
    }else{
      globalpvalue<-try(wald.test(b=model$coef[covariateindex],Sigma=model$var[covariateindex,covariateindex],Terms=seq_along(covariateindex))$result$chi2[3]);
      
    }
      if(class(globalpvalue) == "try-error") globalpvalue<-"NA"
    globalpvalue<-lpvalue(globalpvalue)
    
    if(type=="coxph"|type=="crr"){
      hazardratio<-c(apply(matrix(summary(model)$conf.int[covariateindex,c(1,3,4)],ncol=3),1,psthr))
      pvalues<-c(sapply(summary(model)$coef[covariateindex,5],lpvalue))

    }else if (type=="glm"){
      m<-summary(model)$coefficients      
      hazardratio<-apply(cbind(exp(m[covariateindex,1]),exp(m[covariateindex,1]-1.96*m[covariateindex,2]),
                               exp(m[covariateindex,1]+1.96*m[covariateindex,2])),1,psthr)        
      pvalues<-c(sapply(m[covariateindex,4],lpvalue))
    }else if (type=="lm"|type=="lm"){
      m<-summary(model)$coefficients      
      
      hazardratio<-apply(cbind(m[covariateindex,1],m[covariateindex,1]-1.96*m[covariateindex,2],
                               m[covariateindex,1]+1.96*m[covariateindex,2]),1,psthr)        
      pvalues<-sapply(m[covariateindex,4],lpvalue)
    }else if(type=="lme"){
      m<-summary(model)$tTable      
      hazardratio<-apply(cbind(m[covariateindex,1],m[covariateindex,1]-1.96*m[covariateindex,2],
                               m[covariateindex,1]+1.96*m[covariateindex,2]),1,psthr) 
      
      pvalues<-c(sapply(m[covariateindex,5],lpvalue))
    }
    
    #if not interaction
    
    if(length(betaname[[1]])==1){
      #if cts
      if(!is.factor(data[,oldcovname])){
        
        title<-c(nicename(covariatename),hazardratio,"",globalpvalue)
      }else if(length(levelnames)==1){
        title<-c(covariatename,"","",globalpvalue)        
        if(!is.null(data)) reference<-c(addspace(sanitizestr(names(table(data[,which(names(data)==oldcovname)]))[1])),"reference","","")
        body<-c(levelnames,hazardratio,"","")       
        
        
      }else{
        if(!is.null(data)) reference<-c(addspace(sanitizestr(names(table(data[,which(names(data)==oldcovname)]))[1])),"reference","","")
        title<-c(covariatename,"","",globalpvalue)
        body<-cbind(levelnames,hazardratio,pvalues,rep("",length(levelnames)))     
        
        #if interaction
      }}else{
        if(length(levelnames)!=1){
          title<-c(covariatename,"","",globalpvalue)
          body<-cbind(levelnames,hazardratio,pvalues,rep("",length(levelnames)))
        }else{
          title<-c(covariatename,hazardratio,"",globalpvalue)        
        }
      }    
    
    out<-rbind(title,reference,body)
    rownames(out)<-NULL
    colnames(out)<-NULL
    return(list(out,nrow(out)))
  })
    
  table<-lapply(out,function(x){return(x[[1]])})
  index<-unlist(lapply(out,function(x){return(x[[2]])}))  
  table<-do.call("rbind", lapply(table, data.frame, stringsAsFactors = FALSE))
  
  colnames(table)<-sapply(c("Covariate",sanitizestr(beta),"p-value","Global p-value"),lbld)
  return(table)   
}

#'Print multivariate summary LaTeX table
#'
#'Returns a LaTeX table of the multivariate summary.
#'
#'@param model fitted model object
#'@param data dataframe containing data
#'@keywords print
#'@export

pmvsum<-function(model,data){
  print.xtable(xtable(mvsum(model,data)),include.rownames=F,sanitize.text.function=identity,table.placement="H")
}

#' Convert .TeX to .docx
#' 
#' Covertes the knitr-compiled .TeX file to a .docx file
#' 
#' @param dir full path of .TeX file directory
#' @param fname .TeX file file name. Do not include extension
#' @param pdwd full path to pandoc
#' @param imwd full path to image magick. Only include if there is at least one graphic.
#' @keywords print
#' @export

makedocx<-function(dir,fname,pdwd,imwd=""){
  oldwd<-getwd()
  if(imwd!=""){
    setwd(imwd)
    command<-paste("mogrify -path ", dir,"figure\\ ", "-format png ", dir, "figure\\*.pdf",sep="" )
    shell(command)
  }
  setwd(pdwd)
  command<-paste("pandoc -o ",dir,fname,".docx ",dir,fname,".tex ",
                 "--default-image-extension=png",sep="")
  shell(command)
  setwd(oldwd)
}


#'Plot CI curve
#'
#'Plots a CI curve. Currently not very powerful. Only plots a single curve
#'
#'@param data dataframe containing data
#'@param response character vector or list of character vector. If a list it plot the '1' event for all outcomes on 
#'the same plot 
#'@param group string of the group want to stratify by
#'@param units units of time
#'@param main String corresponding to title
#'@param CI Bool If True will plot CI and only the '1' event. if F will plot all events except for the final one
#'@param legpos string indicating which position to put legend choies are "topright" etc
#'@param xlim numeric vector corresponding to xlimits. Default is NULL
#'@param outcomes character vector of the names of the different competing outcomes 
#'@keywords print
#'@export

plotci<-function (data, response, group=NULL, units = "months",main="Viral Infections",CI=F,legpos="topleft",xlim=NULL,outcomes=NULL){
  if(!is.null(group)){
    groups=levels(data[,group])
  }
  #If response is a list plot the '1' event for all outcomes on same plot
  if(class(response)!="list"){
    if(!is.null(group)){
      groups=levels(data[,group])
      fita <- cuminc(data[, response[1]], data[, response[2]],data[,group])
    }else{
      fita <- cuminc(data[, response[1]], data[, response[2]])
    }
    if(CI){
      plot(fita[[1]]$time, sapply(fita[[1]]$est + 1.96 * sqrt(fita[[1]]$var), 
                                  function(x) min(x, 1)), type = "l", lty = 2, main = paste("CI plot for ", 
                                 sanitizestr(nicename(response[2])), sep = ""), xlab = paste("Time (", 
                                 cap(units), ")", sep = ""), ylim = c(0, 1), ylab = paste("Incidence of ",                                                                                                                                                                                                                
                                 sanitizestr(nicename(response[2])), sep = ""),xlim=xlim)
      
      lines(fita[[1]]$time, fita[[1]]$est)
      lines(fita[[1]]$time, sapply(fita[[1]]$est - 1.96 * sqrt(fita[[1]]$var), 
                                   function(x) max(x, 0)), lty = 2)
    }else{
      plot(fita[[1]]$time, fita[[1]]$est, 
           type = "l",  main = paste("CI plot for ", 
                                     sanitizestr(nicename(response[2])), sep = ""), xlab = paste("Time (", 
                                     cap(units), ")", sep = ""), ylim = c(0, 1), ylab = paste("Incidence of ", 
                                     sanitizestr(nicename(response[2])), sep = ""),xlim=xlim)
      numoutcomes<-length(fita)-1
      if(numoutcomes>1){
        for (i in 2:numoutcomes){
          lines(fita[[i]]$time,fita[[i]]$est,lty=i,lwd=2)
        }
        legend(legpos, outcomes, lty = 1:numoutcomes, bty = "n",lwd=2)  
      }
    }
    
  }else{
    d<-lapply(response,function(respons){
      fita <- cuminc(data[, respons[1]], data[, respons[2]])
      list(fita[[1]]$time, fita[[1]]$est)})
    if(is.null(xlim)) xlim=c(0,ceiling(max(sapply(d,function(x) max(x[[1]])))))
    plot(1,type="n",xlim=xlim,ylim=c(0,1),
         ylab="Cumulative Incidence",xlab = paste("Time (",cap(units), ")", sep = ""),main=paste("Cumulative Incidence plot for",main))
    for(i in 1:length(d)){
      lines(d[[i]][[1]],d[[i]][[2]],lty=i,lwd=2)
    }
    legend(legpos, sapply(response,function(x) x[2]) , col =rep(1, length(response)), lty = 1:length(response), bty = "n",lwd=2)  
    
  }
}




#' Get CI cinfidence interval
#' 
#' Returns the confidence interval of a CI at a specified time. Currently not very powerful. Only works on single strata.
#' 
#' @param data dataframe containing data
#' @param response character vector of response
#' @param times numeric vector specifying single time to get CI for
#' @param units string specifying the unit of times
#' @param outcomes character vector specifying names of competing outcomes.
#' Leave NULL if there is only one outcome
#' @param decimals positive integer corresponding to the number of decimals 
#' @keywords print
#' @export
citime<-function (data, response, times, units="Years",outcomes=NULL,decimals=2) 
{
  out<-sapply(times,function(time){
    fita <- cuminc(data[, response[1]], data[, response[2]])
    numoutcomes<-length(fita)-1
    sapply(1:numoutcomes,function(i){
      index <- max(which(fita[[i]]$time <= time))
      est <- fita[[i]]$est[index]
      pm <- 1.96 * sqrt(fita[[i]]$var[index])
      psthr(c(est, max(est - pm, 0), min(est + pm, 1)),decimals)
    })
  })
  if (class(out)!="matrix")
    out<-t(out)
  out<-data.frame(out,stringsAsFactors=F)
  rownames(out)<-NULL
  if(!is.null(outcomes)){
    out<-cbind(outcomes,out)
    colnames(out)<-c("Outcome",paste(times,units))
  }else{
    colnames(out)<-paste(times,units)
  }
  return(out)
}


#' Create a forrest plot
#' 
#' Create a forrest plot. All entires with cutoff=T will be plotted with an NA
#' rather than their original value.
#' 
#' @param data dataframe containing data
#' @param xlab String corresponding to xlabel. By default is set to names(data)[2]
#' @param ylab String corresponding to ylabel. By default is set to names(data)[1]
#' @param main String corresponding to main title. By default is set to "Forest plot for subgroup analysis"
#' @param space numeric corresponding to offset of y label. Should be positive if y label is on top of the names of the y axis
#' @param bool A boolean vector. All entries with T will be invisible in the plot
#' @param xlim vector of length 2 corresponding to limits of x-axis. Default to NULL.
#' @keywords print
#' @export
forestplot<-function (data, xlab = NULL, ylab = NULL, main = NULL, space = 0, bool=F,xlim=NULL) 
{
  if (is.null(xlab)) 
    xlab <- names(data)[2]
  if (is.null(ylab)) 
    ylab <- names(data)[1]
  if (is.null(main)) 
    main <- "Forest plot for subgroup analysis"
  par(oma = c(0, space, 0, 0))
  l1 <- nrow(data)
  colors<- ifelse(bool,"white","black")
  if(is.null(xlim)) xlim<-c(0,max(data[!bool, 4]))
  plot(data[, 2], c(1:l1), col = colors, pch = "|", bg = colors, 
       yaxt = "n", xlim = xlim, ylab = "", xlab = "", 
       main = main)
  abline(v = 1, col = "red", lty = 2)
  segments(data[, 3], c(1:l1), data[, 4], c(1:l1),col=colors)
  axis(2, at = c(1:l1), labels = data[, 1], las = 1, cex.axis = 0.8)
  mtext(side = 1, xlab, line = 2)
  mtext(side = 2, ylab, line = space + 2.5)
}