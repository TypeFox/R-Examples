compare.clines <-
function(cline.data1=NULL, cline.data2=NULL, sig.test=FALSE,
                         n.reps=1000){
## verify that input data objects were supplied  
  if (is.null(cline.data1)==TRUE | is.null(cline.data2)==TRUE) stop("cline.data was not supplied")
  if (is.list(cline.data1)==FALSE | is.list(cline.data2)==FALSE) stop("cline.data is not a list object")
## make sure both data sets have the same number of loci
  if (dim(cline.data1$Loci.data)[1] != dim(cline.data2$Loci.data)[1])
    stop("the number of markers differs between data sets")
## determine if the range of values for hi.index1 > hi.index2
  if (min(cline.data1$hybrid.index) > min(cline.data2$hybrid.index))
    warning("regression model predictions are being extended beyond the range of the observed values used to fit the model")
  if (max(cline.data1$hybrid.index) < max(cline.data2$hybrid.index))
    warning("regression model predictions are being extended beyond the range of the observed values used to fit the model")
## compute likelihood ratios without significance testing
  if (sig.test==TRUE) cat("compare.clines permutation test working\n this may take a few minutes\n")
  if (sig.test==FALSE){
    ## set up data objects
    n.loci<-dim(cline.data2$Loci.data)[1]
    n.ind<-dim(cline.data2$Count.matrix)[2]
    ## compute fitted values for M1D2
    M1D2.fitted.AA<-array(dim=c(n.loci,n.ind))
    M1D2.fitted.Aa<-array(dim=c(n.loci,n.ind))
    M1D2.fitted.aa<-array(dim=c(n.loci,n.ind))
    for (i in 1:n.loci){
      reg.out<-multinom(cline.data1$Count.matrix[i,]~cline.data1$hybrid.index, trace=FALSE)
      ## for dominant data
      if (cline.data1$Loci.data[i,2]=="D" | cline.data1$Loci.data[i,2]=="d" | cline.data1$Loci.data[i,2]=="H" | cline.data1$Loci.data[i,2]=="h"){
        if (length(coef(reg.out))<2) warning("warning, invariant locus included")
        else{
          AA.slope<-coef(reg.out)[2]
          AA.int<-coef(reg.out)[1]
          Hx<-exp(AA.int+AA.slope*cline.data2$hybrid.index)
          M1D2.fitted.AA[i,]<-exp(AA.int+AA.slope*cline.data2$hybrid.index)/(1+Hx)
          M1D2.fitted.aa[i,]<-1-M1D2.fitted.AA[i,]
        }	
      }
      ## for co-dominant data
      ## this now uses the fit.c.clines function
      else if (cline.data1$Loci.data[i,2]=="C" | cline.data1$Loci.data[i,2]=="c"){
        if (length(coef(reg.out))>2){
          AA.slope<-coef(reg.out)[2,2]
          AA.int<-coef(reg.out)[2,1]
          Aa.slope<-coef(reg.out)[1,2]
          Aa.int<-coef(reg.out)[1,1]
          Hx<-exp(AA.slope*cline.data2$hybrid.index+AA.int) + exp(Aa.slope*cline.data2$hybrid.index+Aa.int)
          M1D2.fitted.AA[i,]<-exp(AA.slope*cline.data2$hybrid.index+AA.int)/(1+Hx)
          M1D2.fitted.Aa[i,]<-exp(Aa.slope*cline.data2$hybrid.index+Aa.int)/(1+Hx)
          M1D2.fitted.aa[i,]<-1-(M1D2.fitted.AA[i,] + M1D2.fitted.Aa[i,])
        }	
        else if (length(coef(reg.out))==2){
          if(sum(cline.data1$Count.matrix[i,]==2,na.rm=TRUE)>=1 & sum(cline.data1$Count.matrix[i,]==1,na.rm=TRUE)>=1){
            AA.slope<-coef(reg.out)[2]
            AA.int<-coef(reg.out)[1]
            Aa.slope<-NA
            Aa.int<-NA
            Hx<-exp(AA.slope*cline.data2$hybrid.index+AA.int)
            M1D2.fitted.AA[i,]<-exp(AA.slope*cline.data2$hybrid.index+AA.int)/(1+Hx)
            M1D2.fitted.Aa[i,]<-1-(exp(AA.slope*cline.data2$hybrid.index+AA.int)/(1+Hx))
            M1D2.fitted.aa[i,]<-rep(0,n.ind)
          }
          else if(sum(cline.data1$Count.matrix[i,]==2,na.rm=TRUE)>=1 & sum(cline.data1$Count.matrix[i,]==0,na.rm=TRUE)>=1){
            AA.slope<-coef(reg.out)[2]
            AA.int<-coef(reg.out)[1]
            Aa.slope<-NA
            Aa.int<-NA
            Hx<-exp(AA.slope*cline.data2$hybrid.index+AA.int)
            M1D2.fitted.AA[i,]<-exp(AA.slope*cline.data2$hybrid.index+AA.int)/(1+Hx)
            M1D2.fitted.Aa[i,]<-rep(0,n.ind)
            M1D2.fitted.aa[i,]<-1-(exp(AA.slope*cline.data2$hybrid.index+AA.int)/(1+Hx))
          }
          else if(sum(cline.data1$Count.matrix[i,]==1,na.rm=TRUE)>=1 & sum(cline.data1$Count.matrix[i,]==0,na.rm=TRUE)>=1){
            AA.slope<-NA
            AA.int<-NA
            Aa.slope<-coef(reg.out)[2]
            Aa.int<-coef(reg.out)[1]
            Hx<-exp(AA.slope*cline.data2$hybrid.index+AA.int)
            M1D2.fitted.AA[i,]<-rep(0,n.ind)
            M1D2.fitted.Aa[i,]<-exp(Aa.slope*cline.data2$hybrid.index+Aa.int)/(1+Hx)
            M1D2.fitted.aa[i,]<-1-(exp(Aa.slope*cline.data2$hybrid.index+Aa.int)/(1+Hx))
          }
        }
        else if (length(coef(reg.out))<2){
          if(sum(cline.data1$Count.matrix[i,]==2,na.rm=TRUE)>=1){
            M1D2.fitted.AA[i,]<-rep(1,n.ind)
            M1D2.fitted.Aa[i,]<-rep(0,n.ind)
            M1D2.fitted.aa[i,]<-rep(0,n.ind)
          }
          if(sum(cline.data1$Count.matrix[i,]==1,na.rm=TRUE)>=1){
            M1D2.fitted.AA[i,]<-rep(0,n.ind)
            M1D2.fitted.Aa[i,]<-rep(1,n.ind)
            M1D2.fitted.aa[i,]<-rep(0,n.ind)
          }
          if(sum(cline.data1$Count.matrix[i,]==0,na.rm=TRUE)>=1){
            M1D2.fitted.AA[i,]<-rep(0,n.ind)
            M1D2.fitted.Aa[i,]<-rep(0,n.ind)
            M1D2.fitted.aa[i,]<-rep(1,n.ind)
          }
        }
      }
    }
    ## compute L(M1|D2)
    ln.likelihood.M1<-numeric(n.loci)
    ln.likelihood.M2<-numeric(n.loci)
    ln.likelihood.ratio<-numeric(n.loci)
    for (i in 1:n.loci){
      probM1D2<-rep(NA,n.ind)
      probM2D2<-rep(NA,n.ind)
      ## for co-dominant markers    
      if (cline.data2$Loci.data[i,2]=="C" | cline.data2$Loci.data[i,2]=="c"){
        for (z in 1:n.ind){
          if (is.na(cline.data2$Count.matrix[i,z])==FALSE){
            if (cline.data2$Count.matrix[i,z]==2){
              probM1D2[z]<-M1D2.fitted.AA[i,z]
              probM2D2[z]<-cline.data2$Fitted.AA[i,z]            
            }
            else if (cline.data2$Count.matrix[i,z]==1){
              probM1D2[z]<-M1D2.fitted.Aa[i,z]
              probM2D2[z]<-cline.data2$Fitted.Aa[i,z]            
            }
            else if (cline.data2$Count.matrix[i,z]==0){
              probM1D2[z]<-M1D2.fitted.aa[i,z]
              probM2D2[z]<-cline.data2$Fitted.aa[i,z]            
            }
          }
          else if (is.na(cline.data2$Count.matrix[i,z])==TRUE){
            probM1D2[z]<-0.33
            probM2D2[z]<-0.33
          }
        }
      }
      else if (cline.data2$Loci.data[i,2]=="D" | cline.data2$Loci.data[i,2]=="d" | cline.data2$Loci.data[i,2]=="H" | cline.data2$Loci.data[i,2]=="h"){
        for (z in 1:n.ind){
          if (is.na(cline.data2$Count.matrix[i,z])==FALSE){
            if (cline.data2$Count.matrix[i,z]==2){
              probM1D2[z]<-M1D2.fitted.AA[i,z]
              probM2D2[z]<-cline.data2$Fitted.AA[i,z]            
            }
            else if (cline.data2$Count.matrix[i,z]==0){
              probM1D2[z]<-M1D2.fitted.aa[i,z]
              probM2D2[z]<-cline.data2$Fitted.aa[i,z]            
            }
          }
          else if (is.na(cline.data2$Count.matrix[i,z])==TRUE){
            probM1D2[z]<-0.5
            probM2D2[z]<-0.5
          }
        }      
      }
      ln.likelihood.M1[i]<-sum(log(probM1D2),na.rm=TRUE)
      ln.likelihood.M2[i]<-sum(log(probM2D2),na.rm=TRUE)
      ln.likelihood.ratio[i]<-ln.likelihood.M2[i]-ln.likelihood.M1[i]
    }
    likelihood.matrix<-array(ln.likelihood.ratio,dim=c(n.loci,1))
    rownames(likelihood.matrix)<-cline.data2$Loci.data[,1]
    colnames(likelihood.matrix)<-"log.likelihood.ratio"
    return (likelihood.matrix)
  }

## conduct permutation significance test, this is a full permutation, minimum and maximum values are not constrained
  else if (sig.test==TRUE){
    ## set up data objects
    n.loci<-dim(cline.data2$Loci.data)[1]
    n.ind<-dim(cline.data2$Count.matrix)[2]
    p.value<-numeric(n.loci)
    ## compute fitted values for M1D2
    M1D2.fitted.AA<-array(dim=c(n.loci,n.ind))
    M1D2.fitted.Aa<-array(dim=c(n.loci,n.ind))
    M1D2.fitted.aa<-array(dim=c(n.loci,n.ind))
    M2D2.fitted.AA<-array(dim=c(n.loci,n.ind))
    M2D2.fitted.Aa<-array(dim=c(n.loci,n.ind))
    M2D2.fitted.aa<-array(dim=c(n.loci,n.ind))
    
    for (i in 1:n.loci){
      reg.out<-multinom(cline.data1$Count.matrix[i,]~cline.data1$hybrid.index, trace=FALSE)
      ## for dominant data
      if (cline.data1$Loci.data[i,2]=="D" | cline.data1$Loci.data[i,2]=="d" | cline.data1$Loci.data[i,2]=="H" | cline.data1$Loci.data[i,2]=="h"){
        if (length(coef(reg.out))<2) warning("warning, invariant locus included")
        else{
          AA.slope<-coef(reg.out)[2]
          AA.int<-coef(reg.out)[1]
          Hx<-exp(AA.int[i]+AA.slope*cline.data2$hybrid.index)
          M1D2.fitted.AA[i,]<-exp(AA.int+AA.slope*cline.data2$hybrid.index)/(1+Hx)
          M1D2.fitted.aa[i,]<-1-M1D2.fitted.AA[i,]
        }	
      }
      ## for co-dominant data
      ## this now uses the fit.c.clines function
      else if (cline.data1$Loci.data[i,2]=="C" | cline.data1$Loci.data[i,2]=="c"){
        if (length(coef(reg.out))>2){
          AA.slope<-coef(reg.out)[2,2]
          AA.int<-coef(reg.out)[2,1]
          Aa.slope<-coef(reg.out)[1,2]
          Aa.int<-coef(reg.out)[1,1]
          Hx<-exp(AA.slope*cline.data2$hybrid.index+AA.int) + exp(Aa.slope*cline.data2$hybrid.index+Aa.int)
          M1D2.fitted.AA[i,]<-exp(AA.slope*cline.data2$hybrid.index+AA.int)/(1+Hx)
          M1D2.fitted.Aa[i,]<-exp(Aa.slope*cline.data2$hybrid.index+Aa.int)/(1+Hx)
          M1D2.fitted.aa[i,]<-1-(M1D2.fitted.AA[i,] + M1D2.fitted.Aa[i,])
        }	
        else if (length(coef(reg.out))==2){
          if(sum(cline.data1$Count.matrix[i,]==2,na.rm=TRUE)>=1 & sum(cline.data1$Count.matrix[i,]==1,na.rm=TRUE)>=1){
            AA.slope<-coef(reg.out)[2]
            AA.int<-coef(reg.out)[1]
            Aa.slope<-NA
            Aa.int<-NA
            Hx<-exp(AA.slope*cline.data2$hybrid.index+AA.int)
            M1D2.fitted.AA[i,]<-exp(AA.slope*cline.data2$hybrid.index+AA.int)/(1+Hx)
            M1D2.fitted.Aa[i,]<-1-(exp(AA.slope*cline.data2$hybrid.index+AA.int)/(1+Hx))
            M1D2.fitted.aa[i,]<-rep(0,n.ind)
          }
          else if(sum(cline.data1$Count.matrix[i,]==2,na.rm=TRUE)>=1 & sum(cline.data1$Count.matrix[i,]==0,na.rm=TRUE)>=1){
            AA.slope<-coef(reg.out)[2]
            AA.int<-coef(reg.out)[1]
            Aa.slope<-NA
            Aa.int<-NA
            Hx<-exp(AA.slope*cline.data2$hybrid.index+AA.int)
            M1D2.fitted.AA[i,]<-exp(AA.slope*cline.data2$hybrid.index+AA.int)/(1+Hx)
            M1D2.fitted.Aa[i,]<-rep(0,n.ind)
            M1D2.fitted.aa[i,]<-1-(exp(AA.slope*cline.data2$hybrid.index+AA.int)/(1+Hx))
          }
          else if(sum(cline.data1$Count.matrix[i,]==1,na.rm=TRUE)>=1 & sum(cline.data1$Count.matrix[i,]==0,na.rm=TRUE)>=1){
            AA.slope<-NA
            AA.int<-NA
            Aa.slope<-coef(reg.out)[2]
            Aa.int<-coef(reg.out)[1]
            Hx<-exp(AA.slope*cline.data2$hybrid.index+AA.int)
            M1D2.fitted.AA[i,]<-rep(0,n.ind)
            M1D2.fitted.Aa[i,]<-exp(Aa.slope*cline.data2$hybrid.index+Aa.int)/(1+Hx)
            M1D2.fitted.aa[i,]<-1-(exp(Aa.slope*cline.data2$hybrid.index+Aa.int)/(1+Hx))
          }
        }
        else if (length(coef(reg.out))<2){
          if(sum(cline.data1$Count.matrix[i,]==2,na.rm=TRUE)>=1){
            M1D2.fitted.AA[i,]<-rep(1,n.ind)
            M1D2.fitted.Aa[i,]<-rep(0,n.ind)
            M1D2.fitted.aa[i,]<-rep(0,n.ind)
          }
          if(sum(cline.data1$Count.matrix[i,]==1,na.rm=TRUE)>=1){
            M1D2.fitted.AA[i,]<-rep(0,n.ind)
            M1D2.fitted.Aa[i,]<-rep(1,n.ind)
            M1D2.fitted.aa[i,]<-rep(0,n.ind)
          }
          if(sum(cline.data1$Count.matrix[i,]==0,na.rm=TRUE)>=1){
            M1D2.fitted.AA[i,]<-rep(0,n.ind)
            M1D2.fitted.Aa[i,]<-rep(0,n.ind)
            M1D2.fitted.aa[i,]<-rep(1,n.ind)
          }
        }
      }
    }
    ## compute L(M1|D2)
    ln.likelihood.M1<-numeric(n.loci)
    ln.likelihood.M2<-numeric(n.loci)
    ln.likelihood.ratio<-numeric(n.loci)
    for (i in 1:n.loci){
      probM1D2<-rep(NA,n.ind)
      probM2D2<-rep(NA,n.ind)
      ## for co-dominant markers    
      if (cline.data2$Loci.data[i,2]=="C" | cline.data2$Loci.data[i,2]=="c"){
        for (z in 1:n.ind){
          if (is.na(cline.data2$Count.matrix[i,z])==FALSE){
            if (cline.data2$Count.matrix[i,z]==2){
              probM1D2[z]<-M1D2.fitted.AA[i,z]
              probM2D2[z]<-cline.data2$Fitted.AA[i,z]            
            }
            else if (cline.data2$Count.matrix[i,z]==1){
              probM1D2[z]<-M1D2.fitted.Aa[i,z]
              probM2D2[z]<-cline.data2$Fitted.Aa[i,z]            
            }
            else if (cline.data2$Count.matrix[i,z]==0){
              probM1D2[z]<-M1D2.fitted.aa[i,z]
              probM2D2[z]<-cline.data2$Fitted.aa[i,z]            
            }
          }
          else if (is.na(cline.data2$Count.matrix[i,z])==TRUE){
            probM1D2[z]<-0.33
            probM2D2[z]<-0.33
          }
        }
      }
      else if (cline.data2$Loci.data[i,2]=="D" | cline.data2$Loci.data[i,2]=="d" | cline.data2$Loci.data[i,2]=="H" | cline.data2$Loci.data[i,2]=="h"){
        for (z in 1:n.ind){
          if (is.na(cline.data2$Count.matrix[i,z])==FALSE){
            if (cline.data2$Count.matrix[i,z]==2){
              probM1D2[z]<-M1D2.fitted.AA[i,z]
              probM2D2[z]<-cline.data2$Fitted.AA[i,z]            
            }
            else if (cline.data2$Count.matrix[i,z]==0){
              probM1D2[z]<-M1D2.fitted.aa[i,z]
              probM2D2[z]<-cline.data2$Fitted.aa[i,z]            
            }
          }
          else if (is.na(cline.data2$Count.matrix[i,z])==TRUE){
            probM1D2[z]<-0.5
            probM2D2[z]<-0.5
          }
        }      
      }
      ln.likelihood.M1[i]<-sum(log(probM1D2),na.rm=TRUE)
      ln.likelihood.M2[i]<-sum(log(probM2D2),na.rm=TRUE)
      ln.likelihood.ratio[i]<-ln.likelihood.M2[i]-ln.likelihood.M1[i]

      ## Permutations for neutral data set
      n.total<-dim(cline.data1$Count.matrix)[2] + dim(cline.data2$Count.matrix)[2]
      allCount<-cbind(cline.data1$Count.matrix,cline.data2$Count.matrix)
      allHI<-c(cline.data1$hybrid.index,cline.data2$hybrid.index)
      perm.ln.likelihood.ratio<-numeric(n.reps)
      for(k in 1:n.reps){
        ## sample individual indexes for permutation
        permData<-sample(1:n.total,n.total,replace=FALSE)
        sam1<-permData[1:dim(cline.data1$Count.matrix)[2]]
        sam2<-permData[(1 + dim(cline.data1$Count.matrix)[2]):n.total]

        ## Regression for first sample
        reg.out<-multinom(allCount[i,sam1]~allHI[sam1], trace=FALSE)
        ## for dominant data
        if (cline.data1$Loci.data[i,2]=="D" | cline.data1$Loci.data[i,2]=="d" | cline.data1$Loci.data[i,2]=="H" | cline.data1$Loci.data[i,2]=="h"){
          if (length(coef(reg.out))<2) warning("warning, invariant locus included")
          else{
            AA.slope<-coef(reg.out)[2]
            AA.int<-coef(reg.out)[1]
            Hx<-exp(AA.int[i]+AA.slope*allHI[sam2])
            M1D2.fitted.AA[i,]<-exp(AA.int+AA.slope*allHI[sam2])/(1+Hx)
            M1D2.fitted.aa[i,]<-1-M1D2.fitted.AA[i,]
          }	
        }
        ## for co-dominant data
        else if (cline.data1$Loci.data[i,2]=="C" | cline.data1$Loci.data[i,2]=="c"){
          if (length(coef(reg.out))>2){
            AA.slope<-coef(reg.out)[2,2]
            AA.int<-coef(reg.out)[2,1]
            Aa.slope<-coef(reg.out)[1,2]
            Aa.int<-coef(reg.out)[1,1]
            Hx<-exp(AA.slope*allHI[sam2]+AA.int) + exp(Aa.slope*allHI[sam2]+Aa.int)
            M1D2.fitted.AA[i,]<-exp(AA.slope*allHI[sam2]+AA.int)/(1+Hx)
            M1D2.fitted.Aa[i,]<-exp(Aa.slope*allHI[sam2]+Aa.int)/(1+Hx)
            M1D2.fitted.aa[i,]<-1-(M1D2.fitted.AA[i,] + M1D2.fitted.Aa[i,])
          }	
          else if (length(coef(reg.out))==2){
            if(sum(allCount[i,sam1]==2,na.rm=TRUE)>=1 & sum(allCount[i,sam1]==1,na.rm=TRUE)>=1){
              AA.slope<-coef(reg.out)[2]
              AA.int<-coef(reg.out)[1]
              Aa.slope<-NA
              Aa.int<-NA
              Hx<-exp(AA.slope*allHI[sam2]+AA.int)
              M1D2.fitted.AA[i,]<-exp(AA.slope*allHI[sam2]+AA.int)/(1+Hx)
              M1D2.fitted.Aa[i,]<-1-(exp(AA.slope*allHI[sam2]+AA.int)/(1+Hx))
              M1D2.fitted.aa[i,]<-rep(0,n.ind)
            }
            else if(sum(allCount[i,sam1]==2,na.rm=TRUE)>=1 & sum(allCount[i,sam1]==0,na.rm=TRUE)>=1){
              AA.slope<-coef(reg.out)[2]
              AA.int<-coef(reg.out)[1]
              Aa.slope<-NA
              Aa.int<-NA
              Hx<-exp(AA.slope*allHI[sam2]+AA.int)
              M1D2.fitted.AA[i,]<-exp(AA.slope*allHI[sam2]+AA.int)/(1+Hx)
              M1D2.fitted.Aa[i,]<-rep(0,n.ind)
              M1D2.fitted.aa[i,]<-1-(exp(AA.slope*allHI[sam2]+AA.int)/(1+Hx))
            }
            else if(sum(allCount[i,sam1]==1,na.rm=TRUE)>=1 & sum(allCount[i,sam1]==0,na.rm=TRUE)>=1){
              AA.slope<-NA
              AA.int<-NA
              Aa.slope<-coef(reg.out)[2]
              Aa.int<-coef(reg.out)[1]
              Hx<-exp(AA.slope*allHI[sam2]+AA.int)
              M1D2.fitted.AA[i,]<-rep(0,n.ind)
              M1D2.fitted.Aa[i,]<-exp(Aa.slope*allHI[sam2]+Aa.int)/(1+Hx)
              M1D2.fitted.aa[i,]<-1-(exp(Aa.slope*allHI[sam2]+Aa.int)/(1+Hx))
            }
          }
          else if (length(coef(reg.out))<2){
            if(sum(allCount[i,sam1]==2,na.rm=TRUE)>=1){
              M1D2.fitted.AA[i,]<-rep(1,n.ind)
              M1D2.fitted.Aa[i,]<-rep(0,n.ind)
              M1D2.fitted.aa[i,]<-rep(0,n.ind)
            }
            if(sum(allCount[i,sam1]==1,na.rm=TRUE)>=1){
              M1D2.fitted.AA[i,]<-rep(0,n.ind)
              M1D2.fitted.Aa[i,]<-rep(1,n.ind)
              M1D2.fitted.aa[i,]<-rep(0,n.ind)
            }
            if(sum(allCount[i,sam1]==0,na.rm=TRUE)>=1){
              M1D2.fitted.AA[i,]<-rep(0,n.ind)
              M1D2.fitted.Aa[i,]<-rep(0,n.ind)
              M1D2.fitted.aa[i,]<-rep(1,n.ind)
            }
          }
          
          ## compute L(M1|D2)
          probM1D2<-rep(NA,n.ind)
          ## for dominant markers
          if (cline.data1$Loci.data[i,2]=="D" | cline.data1$Loci.data[i,2]=="d" | cline.data1$Loci.data[i,2]=="H" | cline.data1$Loci.data[i,2]=="h"){
            for (z in 1:n.ind){
              zSam2<-sam2[z]
              if (is.na(allCount[i,zSam2]==FALSE)){
                if (allCount[i,zSam2]==2) probM1D2[z]<-M1D2.fitted.AA[i,z]
                else if (allCount[i,zSam2]==0) probM1D2[z]<-M1D2.fitted.aa[i,z]
              }
              else if (is.na(allCount[i,zSam2])==TRUE) probM1D2[z]<-0.5
            }
          }
          ## for co-dominant markers
          else if (cline.data1$Loci.data[i,2]=="C" | cline.data1$Loci.data[i,2]=="c"){
            for (z in 1:n.ind){
              zSam2<-sam2[z]
              if (is.na(allCount[i,zSam2])==FALSE){
                if (allCount[i,zSam2]==2) probM1D2[z]<-M1D2.fitted.AA[i,z]
                else if (allCount[i,zSam2]==1) probM1D2[z]<-M1D2.fitted.Aa[i,z]
                else if (allCount[i,zSam2]==0) probM1D2[z]<-M1D2.fitted.aa[i,z]
              }
              else if (is.na(allCount[i,zSam2])==TRUE) probM1D2[z]<-0.33
            }
          }
          
          replnLikM1D2<-sum(log(probM1D2))

          ## Now calculate M2D2

          ## Regression for first sample
          reg.out<-multinom(allCount[i,sam2]~allHI[sam2], trace=FALSE)
          if (cline.data1$Loci.data[i,2]=="D" | cline.data1$Loci.data[i,2]=="d" | cline.data1$Loci.data[i,2]=="H" | cline.data1$Loci.data[i,2]=="h"){
            if (length(coef(reg.out))<2) warning("warning, invariant locus included")
            else{
              AA.slope<-coef(reg.out)[2]
              AA.int<-coef(reg.out)[1]
              Hx<-exp(AA.int[i]+AA.slope*allHI[sam2])
              M2D2.fitted.AA[i,]<-exp(AA.int+AA.slope*allHI[sam2])/(1+Hx)
              M2D2.fitted.aa[i,]<-1-M2D2.fitted.AA[i,]
            }	
          }
          else if (cline.data1$Loci.data[i,2]=="C" | cline.data1$Loci.data[i,2]=="c"){
            if (length(coef(reg.out))>2){
              AA.slope<-coef(reg.out)[2,2]
              AA.int<-coef(reg.out)[2,1]
              Aa.slope<-coef(reg.out)[1,2]
              Aa.int<-coef(reg.out)[1,1]
              Hx<-exp(AA.slope*allHI[sam2]+AA.int) + exp(Aa.slope*allHI[sam2]+Aa.int)
              M2D2.fitted.AA[i,]<-exp(AA.slope*allHI[sam2]+AA.int)/(1+Hx)
              M2D2.fitted.Aa[i,]<-exp(Aa.slope*allHI[sam2]+Aa.int)/(1+Hx)
              M2D2.fitted.aa[i,]<-1-(M2D2.fitted.AA[i,] + M2D2.fitted.Aa[i,])
            }	
            else if (length(coef(reg.out))==2){
              if(sum(allCount[i,sam2]==2,na.rm=TRUE)>=1 & sum(allCount[i,sam2]==1,na.rm=TRUE)>=1){
                AA.slope<-coef(reg.out)[2]
                AA.int<-coef(reg.out)[1]
                Aa.slope<-NA
                Aa.int<-NA
                Hx<-exp(AA.slope*allHI[sam2]+AA.int)
                M2D2.fitted.AA[i,]<-exp(AA.slope*allHI[sam2]+AA.int)/(1+Hx)
                M2D2.fitted.Aa[i,]<-1-(exp(AA.slope*allHI[sam2]+AA.int)/(1+Hx))
                M2D2.fitted.aa[i,]<-rep(0,n.ind)
              }
              else if(sum(allCount[i,sam2]==2,na.rm=TRUE)>=1 & sum(allCount[i,sam2]==0,na.rm=TRUE)>=1){
                AA.slope<-coef(reg.out)[2]
                AA.int<-coef(reg.out)[1]
                Aa.slope<-NA
                Aa.int<-NA
                Hx<-exp(AA.slope*allHI[sam2]+AA.int)
                M2D2.fitted.AA[i,]<-exp(AA.slope*allHI[sam2]+AA.int)/(1+Hx)
                M2D2.fitted.Aa[i,]<-rep(0,n.ind)
                M2D2.fitted.aa[i,]<-1-(exp(AA.slope*allHI[sam2]+AA.int)/(1+Hx))
              }
              else if(sum(allCount[i,sam2]==1,na.rm=TRUE)>=1 & sum(allCount[i,sam2]==0,na.rm=TRUE)>=1){
                AA.slope<-NA
                AA.int<-NA
                Aa.slope<-coef(reg.out)[2]
                Aa.int<-coef(reg.out)[1]
                Hx<-exp(AA.slope*allHI[sam2]+AA.int)
                M2D2.fitted.AA[i,]<-rep(0,n.ind)
                M2D2.fitted.Aa[i,]<-exp(Aa.slope*allHI[sam2]+Aa.int)/(1+Hx)
                M2D2.fitted.aa[i,]<-1-(exp(Aa.slope*allHI[sam2]+Aa.int)/(1+Hx))
              }
            }
            else if (length(coef(reg.out))<2){
              if(sum(allCount[i,sam2]==2,na.rm=TRUE)>=1){
                M2D2.fitted.AA[i,]<-rep(1,n.ind)
                M2D2.fitted.Aa[i,]<-rep(0,n.ind)
                M2D2.fitted.aa[i,]<-rep(0,n.ind)
              }
              if(sum(allCount[i,sam2]==1,na.rm=TRUE)>=1){
                M2D2.fitted.AA[i,]<-rep(0,n.ind)
                M2D2.fitted.Aa[i,]<-rep(1,n.ind)
                M2D2.fitted.aa[i,]<-rep(0,n.ind)
              }
              if(sum(allCount[i,sam2]==0,na.rm=TRUE)>=1){
                M2D2.fitted.AA[i,]<-rep(0,n.ind)
                M2D2.fitted.Aa[i,]<-rep(0,n.ind)
                M2D2.fitted.aa[i,]<-rep(1,n.ind)
              }
            }
          }
    
          ## compute L(M2|D2)
          probM2D2<-rep(NA,n.ind)
          ## for dominant markers
          if (cline.data1$Loci.data[i,2]=="D" | cline.data1$Loci.data[i,2]=="d" | cline.data1$Loci.data[i,2]=="H" | cline.data1$Loci.data[i,2]=="h"){
            for (z in 1:n.ind){
              zSam2<-sam2[z]
              if (is.na(allCount[i,zSam2]==FALSE)){
                if (allCount[i,zSam2]==2) probM2D2[z]<-M2D2.fitted.AA[i,z]
                else if (allCount[i,zSam2]==0) probM2D2[z]<-M2D2.fitted.aa[i,z]
              }
              else if (is.na(allCount[i,zSam2])==TRUE) probM2D2[z]<-0.5
            }
          }
          ## for co-dominant markers
          else if (cline.data1$Loci.data[i,2]=="C" | cline.data1$Loci.data[i,2]=="c"){
            for (z in 1:n.ind){
              zSam2<-sam2[z]
              if (is.na(allCount[i,zSam2])==FALSE){
                if (allCount[i,zSam2]==2) probM2D2[z]<-M2D2.fitted.AA[i,z]
                else if (allCount[i,zSam2]==1) probM2D2[z]<-M2D2.fitted.Aa[i,z]
                else if (allCount[i,zSam2]==0) probM2D2[z]<-M2D2.fitted.aa[i,z]
              }
              else if (is.na(allCount[i,zSam2])==TRUE) probM2D2[z]<-0.33
            }
          }
          replnLikM2D2<-sum(log(probM2D2))
          perm.ln.likelihood.ratio[k]<-replnLikM2D2 - replnLikM1D2
        }
        ## This ends the rep loop
      }
      p.value[i]<-sum(perm.ln.likelihood.ratio > ln.likelihood.ratio[i])/n.reps
      cat("permutation test for locus: ", cline.data1$Loci.data[i,1],"\n")
      ## This ends the locus loop
    }
    ## This ends the locus loop
    likelihood.matrix<-cbind(ln.likelihood.ratio,p.value)
    rownames(likelihood.matrix)<-cline.data2$Loci.data[,1]
    colnames(likelihood.matrix)<-c("log.likelihood.ratio","PValues")
    return (likelihood.matrix)
  }
}

