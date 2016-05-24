vif.cp <-function(data,l,siglev=0.05){
      if(!is.numeric(data))stop("The data must be of numeric type")
      if(!is.numeric(l))stop("The length of the partition shoule be numeric")
      if(length(data)>100000)stop("The length of the data should not exceed 100000")
      n<-length(data)
      if(l<10)stop("The length of partition should be at least 10")
      if(l>n/4)stop("There are at least 4 partition.")
        
      storage.mode(data)<-"double"
      storage.mode(l)<-"integer"
      storage.mode(siglev)<-"double"
      m<-0
      storage.mode(m)<-"integer"
      cpset<-numeric(100)
      storage.mode(cpset)<-"integer"

      storage.mode(n)<-"integer"      
      zzz <- .Fortran("vif",y=data,n=n,l=l,siglev=siglev,cpset=cpset,m=m,PACKAGE="VIFCP")
      cp<-zzz$cpset
      m<-zzz$m
      return(cp[1:m])
}
