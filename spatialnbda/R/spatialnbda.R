#' Formats the data for NBDA
#' @param a diffusion times
#' @param b social network
#' @param c order of acq
#' @param d group
#' @param e diffusion
#' @param f event
#' @param g interaction covariate based on user specified social network
#' @param h spatial covariate


#' @export

FormatData =   function(a,b,c,d,e,f,g,h,networkdata){
  orderacq = as.matrix(c)
  foragetimes = as.matrix(a)
  groups = as.matrix(d)
  diffusions = as.matrix(e)
  n = length(foragetimes)
  events = c(1:n)
  Status = rep(0,n)
  nid = altnid = length(unique(orderacq))
  if(missing(g)){g=matrix(1,nrow=altnid,ncol=altnid)}
  if(missing(g)){nid = altnid}
  if(missing(b)){b=matrix(1,nrow=altnid,ncol=altnid)}
  if(missing(h)){h = rep(0,n)}
  if(missing(networkdata)){networkdata= g}
  b = as.matrix(b)
  g= as.matrix(g)
  h = as.matrix(h)
  nid = nrow(g)
  
  
  
  create_passport = function(index) {
    # in which group am I
    my_group = dataset[index,2]
    # which unique id do I have
    my_id = dataset[index,4]
    # create my passport
    my_passport = paste(my_group,my_id,sep="")
    #my_passport =paste(4,5,sep=".")
    return(my_passport)
  }
  
  iCovariateHomogenous = function(index) {
    mytime = dataset[index,5]
    if(length(mytime)==0){mytime=0}
    groupfilter = which(dataset[,2]==dataset[index,2])
    diffusionfilter= which(dataset[,3]==dataset[index,3]) 
    
    data1 = dataset[diffusionfilter,]
    data2 = data1[which(data1[,2]==dataset[index,2]),] 
    
    if(length(data2[,1]>0)){
      subdataset = data2
      teachers = subdataset[subdataset[,5] < mytime,4] 
      numteachers = length(teachers)
    }else{numteachers=0}
    
    return(numteachers)
  }
  
  
  iCovariate = function(index) {
    mygroup = dataset[index,2]
    
    if(is.list(networkdata)==TRUE){     
      SocialNetwork = networkdata[[mygroup]]     
    }
    
    
    mytime = dataset[index,5]
    if(length(mytime)==0){mytime=0}
    myid = dataset[index,4]
    SocialNetwork=g
    groupfilter = which(dataset[,2]==dataset[index,2])
    diffusionfilter= which(dataset[,3]==dataset[index,3]) 
    data1 = dataset[diffusionfilter,]
    data2 = data1[which(data1[,2]==dataset[index,2]),] 
    
    if(length(data2[,1]>0)){
      subdataset = data2
      teachers = subdataset[subdataset[,5] < mytime,4]   
      numteachers = length(teachers)
    }else{numteachers=0}
    
    
    if(numteachers>0){        
      interactions =  sum(SocialNetwork[myid,c((teachers))])         
    }else{interactions=0}
    
    return(interactions)
  }
  
  
  Followers = function(){
    cents = array(0,c(n,nid))
    for(index in 1:n){
      mytime = dataset[index,5]
      myid = dataset[index,4]
      
      groupfilter = which(dataset[,2]==dataset[index,2])
      diffusionfilter = which(dataset[,3]==dataset[index,3])
      
      data1 = dataset[diffusionfilter,]
      data2 = data1[which(data1[,2]==dataset[index,2]),]
      
      
      if(length(data2[,1]>0)){
        followers = data2
        myfollowers = followers[which(followers[,5]>mytime),4]
      }else{myfollowers=0}
      cents[index,c(myfollowers)] = rep(1,length(myfollowers))
    }
    return(cents)
  }
  
  
  
  dataset = cbind(events,groups,diffusions,orderacq,foragetimes,Status)
  
  
  
  
  
  Censoreds = function(){
    G  = unique(c(dataset[,2]));lg = length(G)
    D = unique(c(dataset[,3])); ld = length(D)
    cargo = data.frame()
    
    for(i in 1: lg){
      
      for(j in 1:ld){
        
        data1 = dataset[which(dataset[,2]==i),]
        data2 = data1[which(data1[,3]==j),]
        myids = unique(c(data2[,4]))
        mytimes = data2[,5]
        candidates = 1:nid
        cs = candidates[-myids]; lcs = length(cs)
        
        if(lcs>0){
          maxtime = max(mytimes)
          csstatus = rep(1,lcs)
          csevents = rep(0,lcs)
          csgroups = rep(i,lcs)
          csdiffusions = rep(j,lcs)
          csID = cs
          cstimes = rep(maxtime+1,lcs)
          csnaives = rep(0,lcs)
          csiC = rep(0,lcs)
          csiCh = rep(0,lcs)
          csBigT = rep(1,lcs)
          entry = cbind(csevents,csgroups,csdiffusions,csID,cstimes,csstatus,csBigT)
          
          cargo = rbind(cargo,entry)
        }
        
        
      }
      
    }
    return(cargo)
  }
  
  
  BigT = function(){   
    bigt = 0
    T = array(0,c(n,1))
    for(index in 1:n){
      mytime = dataset[index,5]
      myid = dataset[index,4]
      
      groupfilter = which(dataset[,2]==dataset[index,2])
      diffusionfilter = which(dataset[,3]==dataset[index,3])
      
      data1 = dataset[diffusionfilter,]
      data2 = data1[which(data1[,2]==dataset[index,2]),]
      
      
      if(length(data2[,1]>0)){
        marker = which(data2[,5]==mytime)
        if(marker[1]==1){bigt = 0}
        if(marker[1]>1){bigt = data2[(marker[1]),5] - data2[(marker[1]-1),5] }
      }
      
      
      T[index] = bigt
    }
    return(T)
  }
  
  T = BigT()
  dataset = cbind(dataset,T)
  censoreds = Censoreds()
  if(dim(censoreds)[1]>0){
    colnames(censoreds)=c("events","groups","diffusions","Ids","times","status","T")
    colnames(dataset)=c("events","groups","diffusions","Ids","times","status","T")
    dataset = rbind(dataset,censoreds)
    statussymbol = dataset[,6]
    dataset[statussymbol==1,1] = length(dataset[statussymbol==0,1]) + c(1:length(dataset[statussymbol==1,1]))
    len = length(dataset[,1])
    
    iCh = iC = passports = array(0,c(len,1))
    for(i in 1:len){
      iCh[i] = iCovariateHomogenous(i)  
      iC[i] = iCovariate(i)
      passports[i] = create_passport(i)
    }
    naives = Followers()
    
    passports = as.numeric(factor(passports))
    #colnames(censoreds)=c("events","groups","diffusions","Ids","times","status","T","iC","iCh")
    dataset = cbind(dataset,iC,iCh,passports)
    colnames(dataset)=c("events","groups","diffusions","Ids","times","status","T","iC","iCh","Passports")
    #dataset = rbind(dataset,censoreds)
  }
  
  
  
  if(dim(censoreds)[1]==0){
    
    colnames(dataset)=c("events","groups","diffusions","Ids","times","status","T")
    statussymbol = dataset[,6]
    dataset[statussymbol==1,1] = length(dataset[statussymbol==0,1]) + c(1:length(dataset[statussymbol==1,1]))
    len = length(dataset[,1])
    
    iCh = iC = passports = array(0,c(len,1))
    for(i in 1:len){
      iCh[i] = iCovariateHomogenous(i)  
      iC[i] = iCovariate(i)
      passports[i] = create_passport(i)
    }
    naives = Followers()
    
    
    
    #naives_list = list()
    #for(i in 1:n){
    # naives_list = list(naives_list,naives[1,])  
    #}
    
    passports = as.numeric(factor(passports))
    #colnames(censoreds)=c("events","groups","diffusions","Ids","times","status","T","iC","iCh")
    dataset = cbind(dataset,iC,iCh,h,passports)
    colnames(dataset)=c("events","groups","diffusions","Ids","times","status","T","iC","iCh","spatialC","Passports")
    colnames(naives)=c(as.character(1:nid))
    #dataset = rbind(dataset,censoreds)
  }
  
  #dataset = cbind(dataset,naives)
  
  list(dataset,naives)
  
}