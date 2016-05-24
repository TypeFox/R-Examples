.aggregate_data <-
function(data,link_maximum,conti,flow="flow",origin="origin",destination="destination")
{
  data[data==link_maximum$origin[1]]=link_maximum$destination[1]
  if (!is.null(conti))
    {
      data=(data[,list(sum(flow),sum(conti)),by=list(origin,destination)])
      setnames(data,c("origin","destination","V1","V2"),c("origin","destination","flow","conti"))
      data$conti=(data$conti>=1)
    }
  else 
    {
    data=(data[,list(sum(flow)),by=list(origin,destination)])
    setnames(data,c("origin","destination","V1"),c("origin","destination","flow"))
    }
  data<-.link_rate(data)
  return(data)
}
.link_max <-
function(data,size_center,conti,list_center)
{
  if (!is.null(conti))
    {
      if (is.null(list_center))
        {data<-data[(data$destination!=data$origin)&(data$sdestinationck_origin<size_center)&(data$conti),]}
      else 
        {data<-data[(data$destination!=data$origin)&(data$sdestinationck_origin<size_center)&(data$conti)&(!(data$origin %in% list_center)),]}
    }
  else
    {
      if (is.null(list_center)) 
        {data<-data[(data$destination!=data$origin)&(data$sdestinationck_origin<size_center),]}
      else
        {data<-data[(data$destination!=data$origin)&(data$sdestinationck_origin<size_center)&(!(data$origin %in% list_center)),]}
    }
  return(data[which.max(data$link),])      
}
.link_rate <-
function(data,flow="flow",origin="origin")
{
  sdestinationck_origin<-(data[,sum(flow),by=list(origin)])
  setnames(sdestinationck_origin,c("origin","V1"),c("origin","sdestinationck_origin"))
  setkey(sdestinationck_origin,"origin")
  setkey(data,"origin")
  data=merge(data,sdestinationck_origin) 
  data$link<-data$flow/data$sdestinationck_origin
  return(data)
}
