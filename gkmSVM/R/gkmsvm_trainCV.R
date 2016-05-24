

# repeated CV
# nCV: number of CV folds
# nrepeat: number of repeated CVs
# cv: CV group label. An array of length (npos+nneg), containing CV group number (between 1 an nCV) for each sequence. (optional)
# C: a vector of all values of C (SVM parameter) to be tested. (optinal)
# showPlots: generate plots (default==TRUE)
# outputPDFfn: filename for output PDF, default=NA (no PDF output)
# outputCVpredfn: filename for output cross validation predictions, default=NA (no output)
# outputROCfn: filename for output auROC (Area Under an ROC Curve) and auPRC (Area Under the Precision Recall Curve) values, default=NA (no output)


gkmsvm_trainCV = function (kernelfn, posfn, negfn, svmfnprfx=NA, nCV=5, nrepeat=1, cv=NA, Type="C-svc", C=1, shrinking=FALSE, showPlots=TRUE, outputPDFfn=NA,  outputCVpredfn=NA, outputROCfn=NA,  ...){
  
  #uses some codes by Dongwon Lee 
  auPRC <- function (perf) {
    rec <- perf@x.values
    prec <- perf@y.values
    result <- list()
    for (i in 1:length(perf@x.values)) {
      result[i] <- list(sum((rec[[i]][2:length(rec[[i]])] - rec[[i]][2:length(rec[[i]])-1])*prec[[i]][-1]))
    }
    return(result)
  }
  
  rocprc <- function(preds,labs, output=NA, showPlots=TRUE) {
    
    linewd <- 1
    wd <- 4
    ht <- 4
    fig.nrows <- 1 
    fig.ncols <- 2
    pt <- 10
    cex.general <- 1 
    cex.lab <- 0.9
    cex.axis <- 0.9
    cex.main <- 1.2
    cex.legend <- 0.8
    
    pred <- ROCR::prediction(preds, labs)
    perf_roc <- ROCR::performance(pred, 'tpr', 'fpr')
    perf_prc <- ROCR::performance(pred, 'prec', 'rec')
    
    perf_auc <- ROCR::performance(pred, 'auc')
    prcs <- auPRC(perf_prc)
    avgauc <- 0
    avgprc <- 0
    
    nCV=length(preds)
    
    for(j in 1:nCV) {
      avgauc <- avgauc + perf_auc@y.values[[j]]
      avgprc <- avgprc + prcs[[j]]
    }        
    
    avgauc <- avgauc/nCV
    avgprc <- avgprc/nCV
    
    if (showPlots){
      if(!is.na(output)){
        pdf(output, width=wd*fig.ncols, height=ht*fig.nrows)
      }
      
      par(xaxs="i", yaxs="i", mar=c(3.5,3.5,2,2)+0.1, mgp=c(2,0.8,0), mfrow=c(fig.nrows, fig.ncols))
      
      ROCR::plot(perf_roc, colorize=FALSE, main="ROC curve", 
                 xlab="1-Specificity", ylab="Sensitivity", cex.lab=1.2,avg='vertical',spread.estimate='stddev')
      text(0.2, 0.1, paste("AUC=", format(avgauc, digits=3, nsmall=3)))
      
      ROCR::plot(perf_prc, colorize=FALSE, main="P-R curve", 
                 xlab="Recall", ylab="Precision", cex.lab=1.2, xlim=c(0,1), ylim=c(0,1),avg='vertical',spread.estimate='stddev')
      text(0.2, 0.1, paste("AUC=", format(avgprc, digits=3, nsmall=3)))
      
      if(!is.na(output)){
        dev.off()
      }
    }
    #cat(paste("auROC=", format(avgauc, digits=3, nsmall=3), '\n'))
    #cat(paste("auPRC=", format(avgprc, digits=3, nsmall=3), '\n'))
    res = c(avgauc,avgprc); 
    names(res)= c('avgauc','avgprc'); 
    return(res);
  }
  
  #posfn="/Users/mghandi/gkmsvm/test/pos500.fa"
  #negfn="/Users/mghandi/gkmsvm/test/neg500.fa"
  #kernelfn= '/Users/mghandi/gkmsvm/test/kernel500.txt'
  #gkmsvm_trainCV(kernelfn, posfn, negfn, svmfnprfx=NA, nCV=5, cv=NA, Type="C-svc", C=10^((-5:5)/2), showPlots=TRUE, outputPDFfn=NA, ...){
  #gkmsvm_trainCV(kernelfn, posfn, negfn )
  #gkmsvm_trainCV(kernelfn, posfn, negfn, C=1, outputPDFfn = '~/plots/ROC.pdf' )
  #gkmsvm_trainCV(kernelfn, posfn, negfn, nCV=5,nrepeat=7, outputPDFfn = '~/plots/ROC.pdf' )
  #gkmsvm_trainCV(kernelfn, posfn, negfn, nCV=5,nrepeat=7, outputPDFfn = '~/plots/ROC2.pdf' )
  
  
  #  library(seqinr)
  #  library(kernlab)
  #  library(utils)
  outres = list(); 
  
  if (requireNamespace("seqinr", quietly = TRUE)&
      requireNamespace("utils", quietly = TRUE)&
      requireNamespace("ROCR", quietly = TRUE)&
      requireNamespace("kernlab", quietly = TRUE)){
        
        
        
        #  negfn= '/Users/mghandi/gkmsvm/test/testneg9.fa'
        #  posfn= '/Users/mghandi/gkmsvm/test/testpos9.fa'
        #  kernelfn= '/Users/mghandi/gkmsvm/test/test9kernel.txt'
        
        pos = seqinr::read.fasta(posfn)
        npos = length(unique(names(pos))) #length(pos)
        neg = seqinr::read.fasta(negfn)
        nneg = length(unique(names(neg))) #length(neg)
        nseq = npos+nneg; 
        
        if(length(which(duplicated(names(pos))))>0){
          print(paste("Error: duplicated sequence ID in", posfn))
          print(names(pos)[which(duplicated(names(pos)))])
          return;
        }
        if(length(which(duplicated(names(neg))))>0){
          print(paste("Error: duplicated sequence ID in", negfn))
          print(names(neg)[which(duplicated(names(neg)))])
          return;
        }
        
        mat <- data.matrix( utils::read.table(file=kernelfn, fill=TRUE, col.names=paste("V", 1:nseq)))
        mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
        rownames(mat)=colnames(mat)
        
        seqnames = c(unique(names(pos)), unique(names(neg))) #c(names(pos), names(neg))
        
        y = c(rep(1, npos), rep(0, nneg)); names(y)=rownames(mat)
        #       cv= sample(nCV, nseq, replace =TRUE)
        if(is.na(cv)){
          cv= 1+sample(nseq, nseq)%%nCV
        }else{
          cv = cv-min(cv)+1
          nCV=max(cv)
        }
        
        if(nCV<2){
          print("Error:  nCV should be >= 2!")
          return;
        }
        if(nrepeat>=nseq){
          nrepeat = nseq-1; print(paste('Warning: nrepeat was reduced to',nrepeat));
        }
        
        aucss=c();
        Cs= C;#10^((-5:5)/2);
        mxaucs=-1;
        for(ic in 1:length(Cs)){
          
          preds = rep(NA,nseq);
          Lpreds <- list()
          Llabs <- list()
          Lseqnams <- list()
          
          idx = 1:nseq; 
          
          for(j in 1:nrepeat){
            jcv = c(cv[-(1:j)],cv[1:j]);
            for(i in 1:nCV){
              ii = which(jcv==i)
              iidx = idx[-ii];
              iK <- kernlab::as.kernelMatrix(mat[iidx, iidx]);
              #isvp <- kernlab::ksvm(iK, y[iidx], type=Type, C=Cs[ic], ...)
              isvp <- kernlab::ksvm(iK, y[iidx], type="C-svc", C=Cs[ic], shrinking=shrinking,...)
              
              alpha = unlist(isvp@alpha )
              kk = iidx[unlist(isvp@SVindex)]
              jj = which(kk>npos); 
              alpha[jj]= -alpha[jj];
              pp= mat[ii,kk]%*%as.matrix(alpha,ncol=1)- isvp@b
              #ll =y[ii] 
              preds[ii]=pp;
              
              Lpreds[[i+(j-1)*nCV]]=pp;
              Llabs[[i+(j-1)*nCV]]=y[ii]
              Lseqnams[[i+(j-1)*nCV]]=seqnames[ii]
            }
          }      
          aucs= rocprc(Lpreds, Llabs,showPlots = FALSE)#, output="~/plots/rocprc.pdf") 
          ri = c(ic,Cs[ic]); names(ri)=c('','C');
          print(c(ri,aucs)); 
          aucss= rbind(aucss,c(ri,aucs) )
          if(aucs[1]>mxaucs){
            mxaucs=aucs[1];
            bestLpreds=Lpreds
            bestLlabs=Llabs
            bestLseqnams=Lseqnams
          }
          
        }
        Copt = aucss[order(-aucss[,3])[1],2]; 
        res = aucss[order(-aucss[,3])[1],2:4];
        
        outres$aucss = res; 
        
        if(!is.na(outputCVpredfn)){
          
          
          cvpred =unlist(bestLpreds)
          names(cvpred) = unlist(bestLseqnams)
          mnpred = rep(0, length(seqnames));
          names(mnpred)=seqnames;
          mi = match(names(cvpred), seqnames)
          labels =unlist(bestLlabs)[match(seqnames,names(cvpred))]
          sdpred = mnpred
          for(i in 1:length(seqnames)){
            ii = which(mi==i)
            mnpred[i]=mean(cvpred[ii])
            sdpred[i]=sd(cvpred[ii])
          }
          #          res = cbind(seqnames,labels, format(round(cbind(mnpred, sdpred),5),nsmall = 5)); 
          res = cbind(seqnames,format(round(mnpred,5),nsmall = 5),2*labels-1,cv-1 ); 
          colnames(res)= c('seqID', 'cvpred_mean', 'label', 'cv_set')
          write.table(res, file=paste(outputCVpredfn), quote = FALSE,row.names = FALSE, col.names=FALSE,sep='\t');
          outres$cvpred = data.frame(seqID=seqnames, cvpred_mean=mnpred, cvpred_sd=sdpred, label =2*labels-1, cv_set=cv-1); 
          #boxplot(as.numeric(res[,3])~as.numeric(res[,2]))
        }
        if(!is.na(outputROCfn)){
          write.table(aucss, file = paste(outputROCfn), col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
        }
        
        if(showPlots){
          if(!is.na(outputPDFfn)){
            pdf(outputPDFfn, width=8, height=4)
          }            
          if(length(Cs)>1){
            par(mfrow=c(1,2))
            plot(aucss[,2],aucss[,3], xlab ='C', ylab='ROC AUC', log='x' ); abline(v=Copt,lty=2)
            plot(aucss[,2],aucss[,4], xlab ='C', ylab='PRC AUC', log='x' ); abline(v=Copt,lty=2)
          }
          rocprc(bestLpreds, bestLlabs,showPlots = TRUE)#, output="~/plots/rocprc.pdf") 
          
          if(!is.na(outputPDFfn)){
            dev.off();
          }
        }
        if (!is.na(svmfnprfx)){
          gkmsvm_train(kernelfn, posfn, negfn, svmfnprfx, Type=Type, C=Copt);
        }
      }
  
    return(outres); 
  
}