path.rpart.new=function(rpart.object,nodes){
   path.ori=path.rpart(rpart.object,nodes)
   path.ori=path.ori[[1]]
   info=data.frame(rep(NA,length(path.ori)-1),rep(NA,length(path.ori)-1),rep(NA,length(path.ori)-1))
   
   for(i in 1:(length(path.ori)-1)){
     path = path.ori[i+1]
     location1 = gregexpr('<',path)[[1]][1]
     location2 = gregexpr('>=',path)[[1]][1]
     location3 = ifelse(location2<0,gregexpr('=',path)[[1]][1],-1)
     status= max(location1,location2,location3)
     category= (location1>0)*0+(location2>0)*1+(location3>0)*2
     info[i,1]=  substr(path,1,status-1)
     info[i,2]=  category
     info[i,3]=  ifelse(location1>0 | location2>0 , as.numeric(substr(path,status+2,nchar(path)))  ,substr(path,status+1,nchar(path)))
        
   }
    return(info)  
}