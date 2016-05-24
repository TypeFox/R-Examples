Plot.cashflow<-function(cashflow)
{
  #return(list(array.rbns.boot=array.rbns.boot,array.ibnr.boot=array.ibnr.boot,
  #            Mat_rbns=Mat_rbns,Mat_ibnr=Mat_ibnr,Mat_total=Mat_total) ) 

  Mat_rbns<-cashflow$Mat.rbns
  Mat_ibnr<-cashflow$Mat.ibnr
  Mat_total<-cashflow$Mat.total
  dd<-ncol(Mat_rbns)
  m<-(dd+1)/2
  par(mfrow=c(3,2))  
  par(mar=c(3,2.5,1.5,1.5),oma=c(2,0.5,0.5,0.2),mgp=c(1.5,0.5,0))
    # c(bottom, left, top, right)
    
  tit<-'RBNS cashflow'
  labelx<-'Future calendar'
  boxplot(Mat_rbns[,-(m:dd)],outline=T,col='gray',names=as.character(1:(m-1)),
            main=tit,xlab=labelx,ylab='')
  hist(Mat_rbns[,dd],breaks=30,main='RBNS reserve',xlab='')
    
  labelx<-'Future calendar'
  tit<-'IBNR cashflow'
  boxplot(Mat_ibnr[,-(m+4:dd)],outline=T,col='gray',names=as.character(1:(m+3)),
          main=tit,xlab=labelx,ylab='')
  hist(Mat_ibnr[,dd],breaks=30,main='IBNR reserve',xlab='')
    
  tit<-'RBNS+IBNR cashflow'
  labelx<-'Future calendar'
  boxplot(Mat_total[,-(m+4:dd)],outline=T,col='gray',names=as.character(1:(m+3)),
            main=tit,xlab=labelx,ylab='')
  hist(Mat_total[,dd],breaks=30,main='Total reserve',xlab='')

  par(mfrow=c(1,1)) 

}  