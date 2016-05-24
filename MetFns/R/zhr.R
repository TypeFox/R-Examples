zhr<-function(data,year,month.beg,month.end=month.beg,day.beg,day.end=day.beg,time.beg=0,
time.end=2359,shw,r=NULL,Ralpha=NULL,Delta=NULL,k1=0.01,k2=1,num,C=1,data2=NULL,
add.plot=FALSE,xlim1=NULL,xlim2=NULL,xinc=NULL,ylim1=NULL,ylim2=NULL,yinc=NULL)
{ 
   if(!is.data.frame(data) || (is.null(r)&&is.null(data2)) || !is.numeric(c(k1,k2,C,num)) || !is.logical(add.plot))
      stop("invalid input parameter(s) specification: check data/r/data2/k1/k2/C/num/add.plot")
   
   data(shw_list,envir=environment())
   shw_list<-get("shw_list",envir=environment())   
   V<-shw_list$V[shw_list$Shw==shw]

   if(!is.null(data2)){
      rdata=pop.index2(data2,year, month.beg, month.end, day.beg, day.end, time.beg, time.end,shw,k1,k2,num)}

   sol1<-solar.long(year,month.beg,day.beg,dec.time(time.beg))
   sol2<-solar.long(year,month.end,day.end,dec.time(time.end))

   data.shw<-filter(data,shw=shw,sol.low=sol1, sol.up=sol2)
   results<-as.data.frame(replicate(9,numeric(0)))
   names(results)<-c("start","stop","sollong","nINT","nSHW",
                     "ZHR","st.err","density","dens.err")

   data.h<-sinh(data.shw,shw,Ralpha,Delta)
   Solar.long<-Vectorize(solar.long)
   Dec.time<-Vectorize(dec.time)
   

   p1=sol1
   
   while(p1<sol2){
         for(i in seq(k1,k2,0.01)){
             p2=p1+i
             data.sel<-filter.sol(data.h,p1,p2)
             if(nrow(data.sel)>0){
             datasel<-data.sel[2*(data.sel$sollong-Solar.long(data.sel$year,
                               data.sel$month,data.sel$day,Dec.time(data.sel$start),3))<=i,]
             nINT<-nrow(datasel)
             if(((sum(datasel$N)>=num)||(p2==p1+k2)) && nINT>0){
                  start<-as.character(sollong_date(year,p1, month.beg,month.end,                                                                                                         day.beg,day.end,time.beg,time.end))
   
                  if(p2<=sol2){
                     stop<-as.character(sollong_date(year,p2,month.beg,month.end,                                                                                                          day.beg,day.end,time.beg,time.end))
                     sollong<-round(mean(c(p1,p2)),3)
                  }
                  else {
                     stop<-as.character(sollong_date(year,sol2,month.beg,month.end,                                                                                                        day.beg,day.end,time.beg,time.end))
                     sollong<-round(mean(c(p1,sol2)),3)}

                 if(is.null(r)){
                    ind<-rdata$sollong==sollong
                    r<-ifelse(any(ind),rdata$pop.index[ind],
                             spline(rdata$sollong,rdata$pop.index,method="natural",xout=sollong)$y)}
          
                 nSHW<-sum(datasel$N)
                 T<-sum((datasel$sine.h*datasel$Teff)/(datasel$F*r^(6.50-datasel$lmg)))
                 ZHR<-(nSHW+C)/T
                 st.err<-ZHR/sqrt(nSHW+C)
                 density<-(10.65*r-12.15)*ZHR/(3600*178700*r^(-1.82)*V)*10^9
                 dens.err<-density*st.err/ZHR
  
                 results=rbind(results,data.frame(start,stop,sollong,nINT,nSHW,ZHR,st.err,density,
                               dens.err))
                 break
   
         }}}
        p1=p2      
       
   }
                                                       
  names(results)[5]<-paste("n",shw,sep="")  
  results[,6:9]=round(results[,6:9],1)
  
  if(add.plot && !is.null(xlim1) && !is.null(xlim2) && !is.null(xinc) && 
     !is.null(ylim1) && !is.null(ylim2) && !is.null(yinc) ) {
     graph.data(results$sollong,results$ZHR,results$st.err,"ZHR(Corrected hourly meteor rate)",
                 xlim1,xlim2,xinc,ylim1,ylim2,yinc)
  }
  
  
  results
}
