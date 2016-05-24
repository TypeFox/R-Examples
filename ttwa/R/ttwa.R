ttwa <-
function(df,origin,destination,flow,conti=NULL,minimum_link=0.00001,size_center,list_center=NULL)
{
  if (!is.null(conti))
    {
      data<-df[,c(origin,destination,flow,conti)]
      names(data)<-c("origin","destination","flow","conti")
    }
  else
    {
      data<-df[,c(origin,destination,flow)]
      names(data)<-c("origin","destination","flow")   
     }
  
  pb_0<-aggregate(flow~origin,data,sum)
  pb_0<-pb_0$origin[pb_0$flow==0]
  data<-data[!(data$origin %in% pb_0 )|!(data$destination %in% pb_0 ),]
  
  
  temp<-data.table(data,key="origin")
  temp<-.link_rate(temp)
  log<-data.frame()
  zoning<-data.frame(id=unique(c(data$origin,data$destination)),zone=unique(c(data$origin,data$destination)),stringsAsFactors = F)
  
  cat("Greedy processus :   \n")
  
  pb <- txtProgressBar(min = 0, max = nrow(zoning), style = 3)
  
  for (i in 1:nrow(zoning))
  {
    link_maximum<-.link_max(temp,size_center,conti,list_center)
    if (nrow(link_maximum)==0){break}
    if (link_maximum$link>=minimum_link)
    {
      temp<-.aggregate_data(temp,link_maximum,conti)
      log<-rbind(log,link_maximum)
      zoning$zone[zoning$zone==link_maximum$origin[1]]=link_maximum$destination[1]
      setTxtProgressBar(pb, i)
    }
    else {break}
  }
  close(pb)
  
  statistic<-aggregate(data.frame(t_outflow=temp$flow,rate_steady=(temp$flow*(temp$origin==temp$destination))),list(zone=temp$origin),sum)
  statistic$rate_steady<-statistic$rate_steady/statistic$t_outflow
  statistic<-merge(statistic,aggregate(data.frame(size=rep(1,nrow(zoning))) ,by=list(zone=zoning$zone),sum))
  
  return(structure(list(zoning=zoning,log=log,data=temp,statistic=statistic),class="zoning"))
}
