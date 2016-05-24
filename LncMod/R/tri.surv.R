tri.surv <-
function(tri,exp.sur,train,test,index=1){
  #a dataframe named "wr" was constructed to contain result
  wr<-as.data.frame(matrix(ncol=15,nrow=length(index)),stringsAsFactors=F)
  colnames(wr)<-c("modulator","effector","target","coef_modulator","p_modulator","coef_effector","p_effector","coef_target","p_target","N_train1","N_train2","dif_train","N_test1","N_test2","dif_test")
  
  #calculation and plot
  for(i in index){
    exp_train<-matrix(nrow=length(train),ncol=3)
    exp_train[,1]<-exp.sur[train,tri[i,1]]
    exp_train[,2]<-exp.sur[train,tri[i,2]]
    exp_train[,3]<-exp.sur[train,tri[i,3]]
    colnames(exp_train)<-tri[i,]
    rownames(exp_train)<-train
    exp_train<-t(exp_train)
    exp_test<-matrix(nrow=length(test),ncol=3)
    exp_test[,1]<-exp.sur[test,tri[i,1]]
    exp_test[,2]<-exp.sur[test,tri[i,2]]
    exp_test[,3]<-exp.sur[test,tri[i,3]]
    colnames(exp_test)<-tri[i,]
    rownames(exp_test)<-test
    exp_test<-t(exp_test)
    #
    q<-summary(coxph(Surv(OS,living)~exp.sur[train,tri[i,1]],exp.sur[train,]))
    wr[i,1]=tri[i,1]
    wr[i,4]=q$coefficients[,1] 
    wr[i,5]=q$coefficients[,5] 
    q<-summary(coxph(Surv(OS,living)~exp.sur[train,tri[i,2]],exp.sur[train,]))
    wr[i,2]=tri[i,2]
    wr[i,6]=q$coefficients[,1] 
    wr[i,7]=q$coefficients[,5] 
    q<-summary(coxph(Surv(OS,living)~exp.sur[train,tri[i,3]],exp.sur[train,]))
    wr[i,3]=tri[i,3]
    wr[i,8]=q$coefficients[,1] 
    wr[i,9]=q$coefficients[,5] 
    
    train_score<-as.data.frame(matrix(ncol=5,nrow=length(train)))
    train_score[,1]<-train
    train_score[,2]<-exp.sur[train,"living"]
    train_score[,3]<-exp.sur[train,"OS"]
    train_score[,4]<-wr[i,5]*exp_train[1,]+wr[i,7]*exp_train[2,]+wr[i,9]*exp_train[3,]
    cuff_train<-median(train_score[,4])
    index1<-train_score[,4]<=cuff_train
    index2<-!(index1)
    train_score[index1,5]<-1
    train_score[index2,5]<-2
    
    index1_sum<-sum(index1)
    index2_sum<-sum(index2)
    wr[i,10]<-index1_sum
    wr[i,11]<-index2_sum
    sort_score<-sort(train_score[,4])

    ord_train<-order(train_score[,4])
    ord_exptrain<-exp_train[,train[ord_train]]
name<-paste("survival",i,tri[i,1],tri[i,2],tri[i,3],sep="_")
    pdf(file=paste(name,"_exp_train",".pdf",sep=""))
    pheatmap(ord_exptrain,cluster_rows=F,cluster_cols=F,cellheight=35,scale="none",show_colnames=F,color = colorRampPalette(c("forestgreen","white","orangered"))(500))
    dev.off()
    
    pdf(file=paste(name,"_train_score",".pdf",sep=""))
    barplot(sort_score,space=0)
    axis(3,at=c(1,index1_sum,length(train))-0.5,labels=c(1,index1_sum,length(train)))
    dev.off()
    
    dif<-survdiff(Surv(train_score[,3],train_score[,2])~train_score[,5],train_score)
    p <- 1-pchisq(dif[[5]],1)
    wr[i,12]<-p
    pdf(file=paste(name,"_train_curve",".pdf",sep=""))
    plot(survfit(Surv(train_score[,3],train_score[,2])~train_score[,5],train_score),col=c('forestgreen','orangered'),lwd=2,xlab="Time:days")
    p<-paste("p-value=",round(p,3))
    text(max(train_score[,3])*0.7,0.9,p,pos=4)
    legend(max(train_score[,3])*0.7,0.85,c("lowrisk","highrisk"),lty=1:1,col=c('forestgreen','orangered'))
    dev.off()

    test_score<-as.data.frame(matrix(ncol=5,nrow=length(test)))
    test_score[,1]<-test
    test_score[,2]<-exp.sur[test,"living"]
    test_score[,3]<-exp.sur[test,"OS"]
    test_score[,4]<-wr[i,5]*exp_test[1,]+wr[i,7]*exp_test[2,]+wr[i,9]*exp_test[3,]
    index1<-test_score[,4]<=cuff_train
    index2<-!(index1)
    test_score[index1,5]<-1
    test_score[index2,5]<-2
    
    index1_sum<-sum(index1)
    index2_sum<-sum(index2)
    wr[i,13]<-index1_sum
    wr[i,14]<-index2_sum
    
    sort_score<-sort(test_score[,4])
    ord_test<-order(test_score[,4])
    ord_exptest<-exp_test[,test[ord_test]]
    pdf(file=paste(name,"_exp_test",".pdf",sep=""))
    pheatmap(ord_exptest,cluster_rows=F,cluster_cols=F,cellheight=35,scale="none",show_colnames=F,color = colorRampPalette(c("forestgreen","white","orangered"))(500))
    dev.off()
    
    pdf(file=paste(name,"_test_score",".pdf",sep=""))
    barplot(sort_score,space=0)
    axis(3,at=c(1,index1_sum,length(test))-0.5,labels=c(1,index1_sum,length(test)))
    dev.off()
    
    dif<-survdiff(Surv(test_score[,3],test_score[,2])~test_score[,5],test_score)
    p <- 1-pchisq(dif[[5]],1)
    wr[i,15]<-p
    pdf(file=paste(name,"_test_curve",".pdf",sep=""))
    plot(survfit(Surv(test_score[,3],test_score[,2])~test_score[,5],test_score),col=c('forestgreen','orangered'),lwd=2,xlab="Time:days")
    p<-paste("p-value=",round(p,3))
    text(max(test_score[,3])*0.7,0.9,p,pos=4)
    legend(max(test_score[,3])*0.7,0.85,c("lowrisk","highrisk"),lty=1:1,col=c('forestgreen','orangered'))
    dev.off()
    }
cat("The plot(s) is on the working directory!\n\n")
    return(wr)
}
