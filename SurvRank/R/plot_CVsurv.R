#' Main function of SurvRank
#'
#' This function creates a pdf file containing some summary plots of the survival ranking analysis
#' @param cv.ob output of the \code{\link{CVrankSurv_fct}}
#' @param data same list used as input in \code{\link{CVrankSurv_fct}}
#' @param file name of the pdf file to be created
#' @param ... other arguments, not used know
#' @keywords SurvRank
#' @export

plot_CVsurv=function(cv.ob,data,file="test.surv.pdf",...){
  pdf(file=file)
  par(cex.lab=1.7,cex.axis=1.7,cex.main=1.7)
  # Plot feature selection method with model selection criterion
  plot(1, type="n", axes=F, xlab="", ylab="")
  text(1,labels=paste(cv.ob$method[1]),cex=3)
  Var1=value=len=auc=m=lowCI=upCI=variable=name=NULL
  seAUC_fct = function(cv.ob){
    # calculates SE for # features selected
    auc_df <- data.frame(len=factor(unlist(lapply(cv.ob$accuracy$used.rank,length))),auc=c(cv.ob$accuracy$auc.out))
    len = levels(auc_df[,1])
    se_a=m_a=NULL
    for(i in 1:length(len)){
      m_a[i] = mean(auc_df[which(len[i]==auc_df[,1]),2],na.rm=T)
      se_a[i] = sd(auc_df[which(len[i]==auc_df[,1]),2],na.rm=T)/sqrt(sum(!is.na(auc_df[which(len[i]==auc_df[,1]),2])))
      if(is.na(se_a[i])){se_a[i]=0}
    }
    ms_a = data.frame(lowCI=m_a-1.96*se_a,m=m_a,upCI=m_a+1.96*se_a,len)
    return(ms_a)
  }
  # Number of features per fold and their selected frequency
    table <- table(unlist(lapply(cv.ob$accuracy$used.rank,length)))
    df1<-data.frame(table)
    df1=reshape::melt(df1)
    p1<-ggplot2::ggplot(data=df1,aes(x=factor(Var1),y=value))+geom_bar(stat="identity",col="darkblue",width=0.1)+xlab("# features")+ylab("freq")+
    ggtitle(paste(cv.ob$method))+theme(text=element_text(size=25),axis.title.y=element_text(size=20))
    print(p1)
    df2 <- data.frame(len=unlist(lapply(cv.ob$accuracy$used.rank,length)),auc=c(cv.ob$accuracy$auc.out))
    #df2=reshape::melt(df2)
  # Number of features per fold and their outer AUCs as scatterplot
    p2.1<-ggplot2::ggplot()+geom_point(data=df2,aes(x=len,y=auc),col="darkblue")+ylab("Survival AUC")+xlab("# features")+ggtitle(paste(cv.ob$method))+
    theme(text=element_text(size=25),axis.title.y=element_text(size=20)) +geom_abline(intercept=0.5,slope=0)
    print(p2.1)
    lim2 = seAUC_fct(cv.ob)
    df2 = merge(x = df2,lim2,by.x = "len",by.y = "len")
  # Number of features per fold and their outer AUCs as boxplot
    p2.2<-ggplot2::ggplot(data=df2,aes(x=factor(len),y=auc))+geom_boxplot(col="darkblue")+ylab("Survival AUC")+xlab("# features")+ggtitle(paste(cv.ob$method))+
    geom_abline(intercept=0.5,slope=0)+theme(text=element_text(size=25),axis.title.y=element_text(size=20))
    p2.2 <- p2.2 + geom_point(aes(factor(len),m),size=3,col="red") + geom_errorbar(aes(factor(len),ymin=lowCI,ymax=upCI),linetype=2,width=0.25,col="red")
    print(p2.2)
  # Boxplot of the AUCs
    df3 = data.frame(auc=apply(cv.ob$accuracy$auc.out,2,mean,na.rm=T))
    df3=reshape::melt(df3)
    p3<-ggplot2::ggplot(data=df3,aes(x=factor(variable),y=value))+geom_boxplot(col="darkblue")+ylab("Mean tAUC")+xlab("")+ggtitle(paste(cv.ob$method))+
    theme(text=element_text(size=25),axis.title.y=element_text(size=20))+stat_summary(fun.y="mean",geom="point",shape=23,size=3,col="darkblue",fill="darkblue")
    print(p3)
  # Heatmap of ranks per fold
    gplots::heatmap.2(cv.ob$rank.mat,Colv = F,trace = "n",dendrogram = "row",breaks=c(1:12))
  # Heatmap of selected features and AUC
    v<-as.vector(cv.ob$accuracy$auc.out)
    mat.auc<-cv.ob$out.mat
    for (i in 1:length(v)){
      mat.auc[,i]<-v[i]*mat.auc[,i]
    }
    colors<-rev(heat.colors(256))
    gplots::heatmap.2(mat.auc,Colv=F,trace="n",dendrogram="row",col=colors)
  # Survival curves
  risk<-riskscore_fct(cv.ob,data)
  # Order scatterplot of selected features and selection frequency
    v<-which(cv.ob$weighted$rel.freq.w>=0.25)
    df4<-data.frame(name=rownames(cv.ob$weighted[v,]),value=cv.ob$weighted[v,1])
    #df4=reshape::melt(df4)
    p4<-ggplot2::ggplot(data=df4,aes(x=factor(c(1:length(value))),y=value,label=name))+geom_point(col="darkblue")+xlab("selected features")+
    ylab("selection frequency")+geom_abline(intercept=0.5,slope=0)+ggtitle(paste(cv.ob$method))+scale_x_discrete(breaks = 1:nrow(df4), labels=df4[,1])+
    theme(axis.text.x = element_text( angle = 90, size = 8),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))
    print(p4)

  dev.off()
}
