pop.index<-function(data,year,month.beg,month.end=month.beg,day.beg,day.end=day.beg,
time.beg=0, time.end=2359,shw, mag.range=-6:7,k,
add.plot=FALSE,xlim1=NULL,xlim2=NULL,xinc=NULL,ylim1=NULL,ylim2=NULL,yinc=NULL)
{ 
   if(!is.data.frame(data) || !is.numeric(k) || !is.logical(add.plot))
      stop("invalid input parameter(s) specification: check data/k/add.plot")  
     
   sol1<-solar.long(year,month.beg,day.beg,dec.time(time.beg),3)
   sol2<-solar.long(year,month.end,day.end,dec.time(time.end),3)

   

   data.shw<-filter(data,shw=shw,sol.low=sol1, sol.up=sol2)
   points=seq(sol1,sol2,by=k)
   Solar.long<-Vectorize(solar.long) 
   Dec.time<-Vectorize(dec.time)

   results<-as.data.frame(matrix(rep(NA,8*(length(points)-1)),ncol=8))
   names(results)=c("start","stop","sollong", "mag","nINT","nSHW",
                                      "pop.index","sigma.r")

   deltam<-(-4:74)/10
   p<-c(0.00046,0.00074,0.0011,0.0016,0.0023,0.0033,0.0046,0.0063,0.0081,0.0100,0.0122,
   0.015,0.018,0.020,0.023,0.026,0.030,0.034,0.039,0.044,0.049,0.056,0.063,0.071,0.079,
   0.088,0.10,0.11,0.13,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.29,0.32,0.35,0.38,0.40,0.43,
   0.46,0.49,0.52,0.53,0.55,0.59,0.64,0.66,0.67,0.69,0.71,0.73,0.74,0.76,0.77,0.78,0.79,
   0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.935,0.94,0.95,
   0.96,0.97,0.98)
   mag.val<--6:7
   

   for(i in 1:(length(points)-1)){
 
       data.sel<-filter.sol(data.shw,points[i],points[i+1])
       results$start[i]<-as.character(sollong_date(year,points[i],month.beg,month.end,                                                                                      day.beg,day.end,time.beg,time.end))
   
       results$stop[i]<-as.character(sollong_date(year,points[i+1],month.beg,month.end,                                                                                    day.beg,day.end,time.beg,time.end))
       
       results$sollong[i]<-round(mean(c(points[i],points[i+1])),3)
       
       if(nrow(data.sel)>0){
  
          datasel<-data.sel[2*(data.sel$sollong-Solar.long(data.sel$year,
                            data.sel$month,data.sel$day,Dec.time(data.sel$start),3))<=k,]
   
          results$nINT[i]<-nrow(datasel)
          results$nSHW[i]<-sum(datasel[,which(names(data.sel)=="N")])
   

          if(nrow(datasel)>0){
  
             mag.distrib<-matrix(rep(0,14*nrow(datasel)),ncol=14)
   
             colnames(mag.distrib)<-names(data.sel)[match("m6",names(datasel)):match("p7",names(datasel))]
   
             for(j in 1:nrow(data.sel)){
                 select<-datasel[j, which(names(datasel)=="m6"):which(names(datasel)=="p7")]
                 deltam.obs<-datasel$lmg[j]-mag.val
   
                 indc<-match(round(deltam.obs,1),deltam)
                 coef<-p[indc]
                 coef[deltam.obs>7.4]<-1
                 mag.distrib[j,]<-as.matrix(select/coef)
                 mag.distrib[j,is.na(mag.distrib[j,])]<-0 }

   
                 counts<-t(apply(mag.distrib,2,sum))
                 cum.freq<-cumsum(counts)

                 ind<-mag.val %in% mag.range & (!counts %in% seq(0,2,by=0.5)) & mag.val<=5

                 m<-mag.val[ind]
                 results$mag[i]<-ifelse(length(m)>0,paste(m[1],":",m[length(m)],sep=""),"-6:7")
                 if(length(m)>4){
                    cum.freq<-cum.freq[ind]
                    mag0.ind<-match(0,m)
                    y<-log(cum.freq/cum.freq[mag0.ind])

                    a_est<-sum(m*y)/sum(m^2)
                    results$pop.index[i]<-exp(a_est)

                    var.est<-1/(length(m)-2)*sum((y-a_est*m)^2)
                    var.a<-var.est/sum(m^2)
                    results$sigma.r[i]<-results$pop.index[i]*sqrt(var.a) }
 
          }}
   }
   
  names(results)[6]<-paste("n",shw,sep="") 
  results$pop.index<-round(results$pop.index,2)
  results$sigma.r<-round(results$sigma.r,2)
  
  if(add.plot && !is.null(xlim1) && !is.null(xlim2) && !is.null(xinc) && 
     !is.null(ylim1) && !is.null(ylim2) && !is.null(yinc) ) {
    graph.data(results$sollong,results$pop.index,results$sigma.r,"Population index",xlim1,xlim2,xinc,ylim1,ylim2,yinc)
  }
  
  results
 
}


