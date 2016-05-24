
mat.ROCcalc  <-  function (truth, data, evalseq,method="trap") 
  {
    numThreasholds = length(evalseq)
    TPF = rep(NA,  numThreasholds)
    TNF = rep(NA,  numThreasholds)
    
    #  1. Get ROC True Positive Fraction (TPF) and True Negative Fraction (TNF) frequencies across range of threasholds
    
    
    #    for (i in 1:numThreasholds) {
    #        pred = ifelse(data <= evalseq[i], 1, 0)
    #        TPF[i] = mean(pred[truth == 1])
    #        TNF[i] = mean(1 - pred[truth == 0])
    #    }
    
    for (i in 1:numThreasholds) {
      pred = ifelse(data <= evalseq[i], 1, 0)
      intx = pred[truth == 1]
      TPF[i] = sum(intx)/length(intx)
      intx = 1 - pred[truth == 0]
      TNF[i] = sum(intx)/length(intx)
    }
    
    
    #  2. Calculate the optimal SCD value that minimizes TPF and FPE
    
    FPE=1-TNF
    slope=(TPF-FPE)
    maxslopeloc=which(slope==max(slope))
    optDissimilarity = evalseq[maxslopeloc]
    
    
    #  3. Calculate AUC of ROC curve using tapezoidal rule
    
    if(method=="trap")
    {
      AUC=-1*polyarea(cbind(c(0,FPE,1,1,0),c(0,TPF,1,0,0)))
      wilcoxon=NULL
    }
    else if (method=="wilcox")
    {
      nfalse=length(which(truth==0))
      ntrue=length(which(truth==1))
      wilcoxon=wilcox.test(data[truth==1],data[truth==0],conf.int=T)
      AUC=1-(wilcoxon$statistic/(nfalse*ntrue))
    }
    #  4. Calculate standard error of AUC using Hanley J.A., McNeil B.J. 1982. Radiology 143:29-36.
    
    Q1=AUC/(2-AUC)
    Q2=(2*(AUC*AUC))/(1+AUC)
    nfalse=length(which(truth==0))
    ntrue=length(which(truth==1))
    denomSE=ntrue*nfalse
    SEarea=sqrt( ((AUC*(1-AUC)) + ((ntrue-1)*(Q1-AUC*AUC))+((nfalse-1)*(Q2-AUC*AUC)) ) /denomSE  )
    
    #  5. Output list object with all ROC information
    
    list( TNF = TNF, TPF = TPF, FPE=FPE, slope=slope, threasholds = evalseq, optDissVal=optDissimilarity, maxPos = maxslopeloc, AUC=AUC, SEAUC=SEarea,method=method,wilcoxon=wilcoxon)
  }
