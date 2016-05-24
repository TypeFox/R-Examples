## This demo creates figure 1
do.fig1<-function(IRDATA){
    age<-as.character(IRDATA$Age.in.yrs.at.first.vaccination)
    antigen<-as.character(IRDATA$Antigen)
    uage<-sort(unique(age))
    pick.fig1<- age<"1" & antigen!="A" & antigen!="C" & antigen!="FIM" &
        antigen!="W135" & antigen!="Y" & antigen!="PRP*"

    alevels<-c("1","3","4","5","6B","7F","9V","14","18C","19F","23F","MenC",
        "PRP","DT","FHA","PRN","PT","TT","HBs","Polio-1","Polio-2","Polio-3")
    fig1.data<-IRDATA[pick.fig1,c("Antigen","n",
        "GMT.95.pct.interval.low.limit","GMT.95.pct.interval.high.limit")]
    cf<-calc.foldrange(fig1.data$n,fig1.data$GMT.95.pct.interval.low.limit,
        fig1.data$GMT.95.pct.interval.high.limit)

    nfig1<-dim(fig1.data)[[1]]
    x<-COL<-rep(NA,nfig1)
    bit<-.4
    bit2<-1
    space<-1
    xlevels<-c((0:10)*bit,10*bit+space+bit2*(0:7),10*bit+space+bit2*7+space+bit*(0:2))
    collevels<-c(rep("blue",11),"slateblue4","aquamarine","green","red","red1","red2","orange","black","pink3","pink","pink1")
    for (i in 1:length(alevels)){
        x[fig1.data$Antigen==alevels[i]]<-xlevels[i]
        COL[fig1.data$Antigen==alevels[i]]<-collevels[i]
    }

    fig1.data<-data.frame(antigen=ordered(fig1.data$Antigen,levels=alevels),
        foldrange=cf[,"foldrange.byt"]) 
    par(oma=c(4,0,0,0))
    plot(range(x),log10(range(fig1.data$foldrange)),type="n",ylab="95% fold range",xlab="",axes=FALSE)
    box()
    axis(2,at=c(1,2,3,4),labels=c("10","100","1000","10000"))
    axis(1,at=xlevels[c(6,12:19,21)],
        labels=c("Pneumonococcal","Meningococcal","HiB",
        "Diptheria Toxin","Pertussis FHA","Pertussis PRN",
        "Pertussis Toxin","Tetanus Toxin","HBSAg","Polio"),las=2,cex=.5)
    points(x,log10(fig1.data$foldrange),col=COL,cex=1.5,pch=18)
}
data(irdata)
do.fig1(irdata)
## Each of the following lines are ways to output the graph in different formats
## change the filename as is appropriate, and remove the # to run 
dev.print(postscript,"H:\\main\\malaria\\saul\\combination\\figures\\figure1.eps",family="Bookman")
#dev.print(pdf,"H:\\main\\malaria\\saul\\combination\\figures\\figure1.pdf")
#dev.print(win.metafile,"H:\\main\\malaria\\saul\\combination\\figures\\figure1.wmf")
 