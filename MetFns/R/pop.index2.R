pop.index2<-function(data,year,month.beg,month.end=month.beg,day.beg,
day.end=day.beg,time.beg=0,time.end=2359,shw,k1=0.01,k2=1,num,
add.plot=FALSE,xlim1=NULL,xlim2=NULL,xinc=NULL,ylim1=NULL,ylim2=NULL,yinc=NULL)
{ 
   if(!is.data.frame(data) || !is.numeric(c(k1,k2,num))|| !is.logical(add.plot))
      stop("invalid input parameter(s) specification: check data/k1/k2/num/add.plot")  
     
   sol1<-solar.long(year,month.beg,day.beg,dec.time(time.beg))
   sol2<-solar.long(year,month.end,day.end,dec.time(time.end))

   data(popind,envir=environment())
   popind<-get("popind",envir=environment())

   data(popind.err,envir=environment())
   popind.err<-get("popind.err",envir=environment())


   data.shw<-filter(data,shw=shw,sol.low=sol1, sol.up=sol2)
   results<-as.data.frame(replicate(7,numeric(0)))
   names(results)<-c("start","stop","sollong","nINT","nSHW",
                                      "pop.index","r.error")

   
   mag.val<--6:7

   
   p1<-sol1
   Solar.long<-Vectorize(solar.long)
   Dec.time<-Vectorize(dec.time)
   
   while(p1<sol2){
         for(i in seq(k1,k2,0.01)){
             p2<-p1+i
             data.sel<-filter.sol(data.shw,p1,p2)
             if(nrow(data.sel)>0){
             datasel<-data.sel[2*(data.sel$sollong-Solar.long(data.sel$year,
                               data.sel$month,data.sel$day,Dec.time(data.sel$start),3))<=i,]
             if(((sum(datasel$N)>=num)||(p2==p1+k2))&& sum(datasel$N)>=10){
                 mag.distrib<-matrix(rep(0,14*nrow(datasel)),ncol=14) 
                 colnames(mag.distrib)<-names(datasel)[match("m6",names(datasel)):match("p7",names(datasel))]
                 coefdm<-rep(0,nrow(datasel))
                 mean.deltam<-0
   
                 for(j in 1:nrow(datasel)){
                     select<-datasel[j, which(names(datasel)=="m6"):which(names(datasel)=="p7")]
                     coefdm[j]<-datasel$N[j]*datasel$lmg[j]-sum(mag.val*select)}

                 mean.deltam<-sum(coefdm)/sum(datasel$N)
                 pop.index<-spline(popind$avdeltam,popind$r,method="natural",xout=mean.deltam)$y
                 r.error<-round(krigeInterp(popind.err$r,log(popind.err$n),popind.err$r.err,xo=pop.index,
                                             yo=log(sum(datasel$N)),extrap=TRUE)$z,6)
         
                 start<-as.character(sollong_date(year,p1,month.beg,month.end,                                                                                            day.beg,day.end,time.beg,time.end))
   
                 if(p2<=sol2){
                    stop<-as.character(sollong_date(year,p2,month.beg,month.end,                                                                                             day.beg,day.end,time.beg,time.end))
                    sollong<-round(mean(c(p1,p2)),3)
                  }
                 else {
                     stop<-as.character(sollong_date(year,sol2,month.beg,month.end,                                                                                           day.beg,day.end,time.beg,time.end))
                     sollong<-round(mean(c(p1,sol2)),3)}
         
   
   
                 nINT<-nrow(datasel)
                 nSHW<-sum(datasel$N)
  
                 results<-rbind(results,data.frame(start,stop,sollong,nINT,nSHW,pop.index,r.error))
                 break
   
        }}}
       p1<-p2
                 
         
       
   }
  names(results)[5]<-paste("n",shw,sep="") 
  results$pop.index<-round(results$pop.index,2)
  results$r.error<-round(results$r.error,2)
  
  if(add.plot && !is.null(xlim1) && !is.null(xlim2) && !is.null(xinc) && 
     !is.null(ylim1) && !is.null(ylim2) && !is.null(yinc) ) {
     graph.data(results$sollong,results$pop.index,results$r.error,"Population index",xlim1,xlim2,xinc,ylim1,ylim2,yinc)
  }


  results
 
}


