CalcDist<-function(data,index,attrs.nominal)
{
  dd=dim(data)
  len=0
  if(length(attrs.nominal)==0)
  {
    data.num=data
  }
  else
  {
    index.nominal=which(colnames(data)%in%colnames(data[,attrs.nominal,drop=FALSE]))
    data.num=data[,-index.nominal]
  }

  #handle the numeric features
  if(length(attrs.nominal)<(dd[2]-1))
  {
    data.num=data.num[,-ncol(data.num)]
    #data.num=normalize(data.num)
    dn=dim(data.num)
    na.feature=which(is.na(data.num[index,]))
  }

  if(length(attrs.nominal)>0)
  {
    no.feature=which(is.na(data[index,attrs.nominal]))
  }

    vrem.dist=sapply(1:dd[1], function(z) {

      rast=0
      len=0
    if(length(attrs.nominal)<(dd[2]-1))
    {
    na.feature1=which(is.na(data.num[z,]))
    na.end=union(na.feature,na.feature1)
    len=dn[2]-length(na.end)
    #check na.end==dn[2]
    if(len!=0)
    {
      if(length(na.end)==0)
      {
        rast=data.num[z,]-data.num[index,]
      }
      else
      {
        rast=data.num[z,-na.end]-data.num[index,-na.end]
      }
    rast=rast**2
    rast=sum(rast)
    }
    }

    sum=0
    len.no=0
    if(length(attrs.nominal)>0)
    {
    no.feature1=which(is.na(data[z,attrs.nominal]))
    na.end=union(no.feature,no.feature1)
    len.no=length(attrs.nominal)-length(na.end)
    #check na.end==length(attrs.nominal)
    if(len.no!=0)
    {
      if(length(na.end)==0)
      {
        sum=data[index,attrs.nominal]!=data[z,attrs.nominal]
      }
      else
      {
        sum=data[index,attrs.nominal[-na.end]]!=data[z,attrs.nominal[-na.end]]
      }
      sum=length(which(sum))
    }
    }

    if((len==0)&&(len.no==0))
    {
       dist=NA
    }
    else
    {
      len=len+len.no
      dist=(rast+sum)/len
      dist=sqrt(dist)
    }
  })

  vrem.dist
}

generate.data.miss<-function(data,percent=5,filename=NULL)
{
  #add factors - nominal attributes
  n.nom=sample(1:2,nrow(data),replace=TRUE)
  data=cbind(data[,-ncol(data)],n.nom,data[,ncol(data)])
  colnames(data)[c(ncol(data)-1,ncol(data))]=c("Nominal","Class")
  data[,ncol(data)-1]=as.factor(data[,ncol(data)-1])

  dd=dim(data)

  percent=ceiling(percent*dd[1]*(dd[2]-1)/100)
  set.seed(123)
  index.row=sample(1:dd[1],percent,replace=TRUE)

  vrem=table(index.row)
  for(i in 1:length(vrem))
  {
    index.col=sample(1:(dd[2]-1),vrem[i])
    data[as.numeric(names(vrem)[i]),index.col]=NA
  }

  len=which(is.na(data))
  #write.table(data,file="leukemia_miss.txt",sep='\t',row.names = FALSE)
  if(length(filename)>0)
  {
    write.table(data,file=filename,sep='\t',row.names = FALSE)
  }
  return(data)
}

input_miss<-function(matrix,method.subst="near.value",attrs.nominal=numeric(),delThre=0.2)
{
  dd=dim(matrix)

  if(length(attrs.nominal)>0)
  {
    for(i in 1:length(attrs.nominal))
    {
      matrix[,attrs.nominal[i]]=as.factor(matrix[,attrs.nominal[i]])
    }
  }

  # data=matrix[,-c(attrs.nominal,dd[2])]
  data=matrix
  dd=dim(data)


  flag.miss=TRUE
  #delete the genes
  na_values=sapply(1:(dd[2]-1), function(z) length(which(is.na(data[,z]))) )
  index=which(na_values>delThre*dd[1])
  #global matrix
  if(length(index)<(dd[2]-1))
  {
  if(length(index)>0)
  {
    data=data[,-index]
  }

  #substitute missing values
  dd=dim(data)
  switch(method.subst,
         #there is the error when the number of features=1
         del = {


         },
         mean.value={

            index.nominal=which(colnames(data)%in%colnames(matrix[,attrs.nominal,drop=FALSE]))
           if(length(index.nominal)>0)
           {
             index.number=setdiff(1:(dd[2]-1),index.nominal)
             for(ij in 1:length(index.nominal))
             {
               na.feature=which(is.na(data[,index.nominal[ij]]))
               if(length(na.feature)>0)
               {
               tt=table(data[-na.feature,index.nominal[ij]])
               tt1=names(tt)[which.max(tt)]
               data[na.feature,index.nominal[ij]]=tt1
               }
             }
           }
           else
           {
             index.number=1:(dd[2]-1)
           }

            if(length(index.number)>0)
            {

              mean.value=apply(data[,index.number,drop=FALSE],2,mean,na.rm=TRUE)

            for(ij in 1:length(index.number))
            {
              na.feature=which(is.na(data[,index.number[ij]]))
              if(length(na.feature)>0)
              {
                data[na.feature,index.number[ij]]=mean.value[ij]
              }

            }
            }
         },
         median.value={

         },
         near.value={

           index.nominal=which(colnames(data)%in%colnames(matrix[,attrs.nominal,drop=FALSE]))
           d.vrem=data[,-c(index.nominal,ncol(data))]

           #normalize the numeric data
           normalize <- function(x) {
             x <- as.matrix(x)
             minAttr=apply(x, 2, min,na.rm=TRUE)
             maxAttr=apply(x, 2, max,na.rm=TRUE)
             #x<-x-rep(minAttr,each=nrow(x))
             #x<-x/rep(maxAttr-minAttr,each=nrow(x))
             x <- sweep(x, 2, minAttr, FUN="-")
             x=sweep(x, 2,  maxAttr-minAttr, "/")
             attr(x, 'normalized:min') = minAttr
             attr(x, 'normalized:max') = maxAttr
             return (x)
           }

           d.vrem=normalize(d.vrem)
           data[,-c(index.nominal,ncol(data))]=d.vrem



           for(i in 1:dd[1])
           {
             if(any(is.na(data[i,-ncol(data)])))
             {
               na.feature=which(is.na(data[i,-ncol(data)]))
               dist.value=CalcDist(data,i,attrs.nominal)
               index.no.na=which(!is.na(dist.value))
               index.no.na=index.no.na[-i]
               sort.dist=sort(dist.value[index.no.na],index.return=TRUE)



               vrem=0
                 len=length(index.no.na)
                 #check for len==0
                 for(k in 1:length(na.feature))
                 {
                   vrem=0
                   vrem1=0
                   iter=0
                   nom.value=NULL

                   for(jk in 1:len)
                   {
                     #10 nearest neighbors
                     if(iter==10) break
                     dat.vrem=data[index.no.na[sort.dist$ix[jk]],na.feature[k]]
                      if(!is.na(dat.vrem))
                      {
                        if(names(data)[na.feature[k]]%in%colnames(data[,attrs.nominal,drop=FALSE]))
                        {
                          nom.value=c(nom.value,dat.vrem)
                        }
                        else
                        {
                          vrem=vrem+(1/(sort.dist$x[jk]+0.01))*dat.vrem
                          vrem1=vrem1+1/(sort.dist$x[jk]+0.01)
                        }
                      iter=iter+1
                      }
                   }

                   if(iter==0) {
                     #vrem=NA
                     flag.miss=FALSE
                   }
                   else
                   {
                   if(names(data)[na.feature[k]]%in%colnames(data[,attrs.nominal,drop=FALSE]))
                   {
                     n.value=table(nom.value)
                     n.index=which.max(n.value)
                     data[i,na.feature[k]]=names(n.value)[n.index]
                   }
                   else
                   {
                     data[i,na.feature[k]]=vrem/vrem1
                   }
                   }
                 }
             }
             #Iteration
             #cat(paste("Iteration",i,"\n"))
           }


           #reverse normalization
           #normalize the numeric data
           unnormalize <- function(x,minx,maxx) {
             x <- as.matrix(x)

             #x<-x*rep(maxAttr-minAttr,each=nrow(x))
             #x<-x+rep(minAttr,each=nrow(x))
             x=sweep(x, 2,  maxx-minx, "*")
             x <- sweep(x, 2, minx, FUN="+")
             return (x)
           }

           minAttr=attr(d.vrem,'normalized:min')
           maxAttr=attr(d.vrem,'normalized:max')

           data[,-c(index.nominal,ncol(data))]=unnormalize(data[,-c(index.nominal,ncol(data))],minAttr,maxAttr)
         }
  )
  }
  else
  {
    #delete all features
    flag.miss=FALSE
  }


  return(list(data=data,flag.miss=flag.miss))
}
