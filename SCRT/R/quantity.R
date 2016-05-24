quantity<-function(design,MT,limit,starts=file.choose(new=FALSE)){

  if(design=="AB"|design=="ABA"|design=="ABAB"){
    phase<-length(strsplit(design,c(""))[[1]])
    quantity<-choose(MT-phase*limit+phase-1,phase-1)
  }
  
  if(design=="CRD"){
    quantity<-choose(MT,MT/2)
  }
  
  if(design=="RBD"){
    quantity<-2^(MT/2)
  }
  
  if(design=="ATD"){
    quantityCRD<-choose(MT,MT/2)
    if(MT<=20){
      index<-1:MT 
      index.a<-matrix(combn(index,(MT/2)),ncol=quantityCRD)
      index.b<-matrix(index.a[,ncol(index.a):1],ncol=quantityCRD)
      if(MT/2<=limit){
        stop<-1:ncol(index.a)
      }
      if(MT/2>limit){
        dist.a<-numeric()
        for(it in 2:nrow(index.a)){
          dist.a<-rbind(dist.a,index.a[it,]-index.a[it-1,])
        }
        dist.b<-numeric()
        for(it in 2:nrow(index.b)){
          dist.b<-rbind(dist.b,index.b[it,]-index.b[it-1,])
        }
        dist.check.a<-dist.a==1
        dist.check.b<-dist.b==1
        sum.a<-numeric()
        for(itr in limit:nrow(dist.check.a)){
          sum.a2<-0
          for(itr2 in 1:limit){
            sum.a2<-sum.a2+dist.check.a[itr-itr2+1,]
          }
          sum.a<-rbind(sum.a,sum.a2)
        }
        sum.b<-numeric()
        for(itr in limit:nrow(dist.check.b)){
          sum.b2<-0
          for(itr2 in 1:limit){
            sum.b2<-sum.b2+dist.check.b[itr-itr2+1,]
          }
          sum.b<-rbind(sum.b,sum.b2)
        }
        sum.a.check<-sum.a==limit
        sum.b.check<-sum.b==limit
        sum.rows.a<-numeric()
        for(it in 1:ncol(sum.a.check)){
          sum.rows.a[it]<-sum(sum.a.check[,it])
        }
        sum.rows.b<-numeric()
        for(it in 1:ncol(sum.b.check)){
          sum.rows.b[it]<-sum(sum.b.check[,it])
        }
        check.stop<-sum.rows.a+sum.rows.b!=0
        stop<-order(check.stop)[1:sum(check.stop==F)]
      }
      quantity<-length(stop)
    }
    if(MT>20){
      fileATD<-tempfile(pattern="ATDassignments",tmpdir=tempdir())
      fileCRD<-tempfile(pattern="CRDassignments",tmpdir=tempdir())
      N<-c(rep(0,MT/2),rep(1,MT/2))
      assignment<-matrix(0,ncol=MT)
      assignment<-rbind(rep(c(0,1),MT/2))
      write.table(assignment,file=fileCRD,append=TRUE,col.names=FALSE,row.names=FALSE)
      CRD<-read.table(fileCRD)
      write.table(assignment,file=fileATD,append=TRUE,col.names=FALSE,row.names=FALSE)
      assignments<-read.table(fileATD)
      repeat{
        assignment<-matrix(0,ncol=MT)
        assignment<-rbind(sample(N,MT,replace=FALSE))
        copy<-numeric()
        for(itr in 1:nrow(CRD)){
          copy2<-numeric(MT)
          for(it in 1:MT){
            copy2[it]<-assignment[1,it]==CRD[itr,it]
          }
          copy<-c(copy,prod(copy2))
        }
        if(sum(copy)==0){
          write.table(assignment,file=fileCRD,append=TRUE,col.names=FALSE,row.names=FALSE)
          CRD<-read.table(fileCRD)
          check<-numeric()
          for(itr in 1:(MT-limit)){
            check2<-0
            for(it in itr:(itr+limit)){
              check2<-check2+assignment[,it]
            }
            check<-cbind(check,check2)
          }
          if(sum(check==(limit+1)|check==0)==0){
            write.table(assignment,file=fileATD,append=TRUE,col.names=FALSE,row.names=FALSE)
            assignments<-read.table(fileATD)
          }
        }
        if(nrow(CRD)==quantityCRD)break
      }
      quantity<-nrow(assignments)
      unlink(fileCRD,recursive=FALSE)
      unlink(fileATD,recursive=FALSE)
    }
  }
  
  if(design=="MBD"){
    points<-read.table(starts)
    N<-nrow(points)
    readLines(con=starts,n=N)->startpoints
    limits<-list()
    for(it in 1:N){
      limits[[it]]<-startpoints[it]
    }
    for(it in 1:N){
      limits[[it]]<-strsplit(limits[[it]],"\t")
    }
    number<-numeric(N)
    for(it in 1:N){
      number[it]<-length(limits[[it]][[1]])
    }
    quantity<-factorial(N)*prod(number)
  }
  
  return(quantity)

}