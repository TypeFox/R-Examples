#####################################################################
#This function tests two of K markers in a pair-wise manner
# and do multiple testing adjustment
######################################################################

Pairwise.DiagTest3Grp.Test <- function(dat,paired=FALSE,type=c("VUS","Youden"),p=0,q=0,mu=0,conf.level=0.95,alternative=c("two.sided","less","greater"),p.adjust.method=c("none","holm","hochberg","hommel","bonferroni","BH","BY","fdr"),digits=3,...){
            
    type <- match.arg(type)##allows incomplete input in the argument, e.g., "V" for "VUS" and "Y" for "Youden"
    
    if(!type%in%c("VUS","Youden")) stop("type must be VUS or Youden!!")
    
    dat.str <- class(dat)
    if(dat.str=="list" & paired) stop("Paired data (where markers measured on the same setof subjects) should be organized as a data frame with samples at row and group& marker measurements at column, not a list!")
    if(dat.str=="data.frame" & !paired) stop("unpaired data (where markers measured on different sets subjects) should be organized as a list (with each component for a marker), not a data frame!")

    ####check data structure
    if(dat.str=="list")
      {
        K <- length(dat)
        markerID <- names(dat)
      }
    if(dat.str=="data.frame")
      {
        K <- ncol(dat)-1##the first column is group membership
        markerID <- names(dat)[-1]
      }
    if(K<2) stop("must have >=2 markers!")
    
    p.adjust.method <- match.arg(p.adjust.method)

    ###do pairwise testing
    print.mat <- matrix(NA,3*(K-1),K-1)
    pval.mat <- NULL
    for (i1 in 1:(K-1))
      {
        for(i2 in (i1+1):K)
          {
            if(dat.str=="list") dat0 <- dat[c(i1,i2)]
            if (dat.str=="data.frame") dat0 <- dat[,c(1,i1+1,i2+1)]            
            test.2marker <- DiagTest3Grp.Test(dat=dat0,paired=paired,type=type,p=p,q=q,mu=mu,conf.level=0.95,alternative=alternative)
            #print(test.2marker)
            stat <- test.2marker$statistic
            print.mat[(i1-1)*3+1,i2-1] <- round(stat,digits=3)
            pval <- test.2marker$p.value
            pval.mat <- rbind(pval.mat,data.frame(i1=i1,i2=i2,marker1=markerID[i1],summary.measure1=test.2marker$estimate[1],marker2=markerID[i2],summary.measure2=test.2marker$estimate[2],CI.lower=as.numeric(test.2marker$conf.int[1]),CI.upper=as.numeric(test.2marker$conf.int[2]),statistic=stat,p.value=pval))
            print.mat[(i1-1)*3+2,i2-1] <- format.pval(pval,digits=digits)

          }
      }
    
    ##multiple testing adjustment
    adjust.p.value <- p.adjust(p=pval.mat$p.value,method=p.adjust.method)
    pval.mat <- data.frame(pval.mat,adjust.p.value=adjust.p.value)
    rownames(pval.mat) <- NULL

    ###add multiple testing adjusted p value at the 3rd row for each marker
    for(idx in 1:nrow(pval.mat))
      {

        print.mat[(pval.mat$i1[idx]-1)*3+3,pval.mat$i2[idx]-1] <- format(pval.mat$adjust.p.value[idx],digits=digits)
      }
    print.mat <- data.frame(print.mat,stringsAsFactors=F)
    names(print.mat) <- markerID[-1]
    print.mat <- data.frame(MarkerID=c(sapply(markerID[-length(markerID)],function(x) c(x,NA,NA))),Row=rep(c("Z-statistic","raw P","adjusted P"),length(markerID)-1),print.mat,stringsAsFactors=F)

    ###make a heatmap to exhibit the p-values after adjustment
    temp.mat2 <- print.mat[seq(3,nrow(print.mat),by=3),-c(1,2)]
    for (ii in 1:ncol(temp.mat2)) storage.mode(temp.mat2[,ii]) <- "double"
    markerID <- print.mat$MarkerID[seq(1,nrow(print.mat),by=3)]
    rownames(temp.mat2) <- markerID
    heatmap.2(as.matrix(temp.mat2),symm=T,Rowv=FALSE, dendrogram="none",cellnote=temp.mat2,trace="none",main="Adjusted P-values",...)
    
    return(list(print.matrix=print.mat,pval.matrix=pval.mat))
    
  }
