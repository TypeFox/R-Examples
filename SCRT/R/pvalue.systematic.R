pvalue.systematic<-function(design,statistic,save="no",limit,data=read.table(file.choose(new=FALSE)),starts=file.choose(new=FALSE)){

  if(design=="CRD"){
    observed.a<-data[,2][data[,1]=="A"]
    observed.b<-data[,2][data[,1]=="B"]
    MT<-nrow(data)
    quantity<-choose(MT,MT/2)
    if(MT<=22){
      if(save=="yes"){
        file<-file.choose(new=FALSE)
      }
      observed<-data[,2]
      scores.a<-combn(observed,MT/2)
      mean.a<-numeric(ncol(scores.a))
      for(it in 1:ncol(scores.a)){
        mean.a[it]<-mean(scores.a[,it])
      }
      mean.b<-rev(mean.a)
      if(statistic=="A-B"){
        distribution<-numeric(length(mean.a))
        for(it in 1:length(mean.a)){
          distribution[it]<-mean.a[it]-mean.b[it]
        }
        distribution<-sort(distribution)
        observed.statistic<-mean(observed.a)-mean(observed.b)
      }
      if(statistic=="B-A"){
        distribution<-numeric(length(mean.a))
        for(it in 1:length(mean.a)){
          distribution[it]<-mean.b[it]-mean.a[it]
        }
        distribution<-sort(distribution)
        observed.statistic<-mean(observed.b)-mean(observed.a)
      }
      if(statistic=="|A-B|"){
        distribution<-numeric(length(mean.a))
        for(it in 1:length(mean.a)){
          distribution[it]<-abs(mean.a[it]-mean.b[it])
        }
        distribution<-sort(distribution)
        observed.statistic<-abs(mean(observed.a)-mean(observed.b))
      }
      test<-distribution>=observed.statistic
      p.value<-sum(test)/quantity
      if(save=="yes"|save=="check"){
        write.table(distribution,file=file,col.names=FALSE,row.names=FALSE)
      }
      return(p.value)
    }
    if(MT>22){
      fileCRD<-tempfile(pattern="CRDassignments",tmpdir=tempdir())
      file.a<-tempfile(pattern="ascores",tmpdir=tempdir())
      file.b<-tempfile(pattern="bscores",tmpdir=tempdir())
      if(save=="yes"){
        file<-file.choose(new=FALSE)
      }
      N<-c(rep("A",MT/2),rep("B",MT/2))
      assignment<-matrix(0,ncol=MT)
      assignment<-rbind(sample(N,MT,replace=FALSE))
      write.table(assignment,file=fileCRD,append=TRUE,col.names=FALSE,row.names=FALSE)
      assignments<-read.table(fileCRD)
      assignment<-as.vector(assignment)
      score.a<-data[,2][assignment=="A"]
      score.a<-t(as.matrix(score.a))
      write.table(score.a,file=file.a,append=TRUE,col.names=FALSE,row.names=FALSE)
      scores.a<-read.table(file.a)
      score.b<-data[,2][assignment=="B"]
      score.b<-t(as.matrix(score.b))
      write.table(score.b,file=file.b,append=TRUE,col.names=FALSE,row.names=FALSE)
      scores.b<-read.table(file.b)
      repeat{
        assignment<-matrix(0,ncol=MT)
        assignment<-rbind(sample(N,MT,replace=FALSE))
        copy<-numeric()
        for(itr in 1:nrow(assignments)){
          copy2<-numeric(MT)
          for(it in 1:MT){
            copy2[it]<-assignment[1,it]==assignments[itr,it]
          }
          copy<-c(copy,prod(copy2))
        }
        if(sum(copy)==0){
          write.table(assignment,file=fileCRD,append=TRUE,col.names=FALSE,row.names=FALSE)
          assignments<-read.table(fileCRD)
          assignment<-as.vector(assignment)
          score.a<-data[,2][assignment=="A"]
          score.a<-t(as.matrix(score.a))
          write.table(score.a,file=file.a,append=TRUE,col.names=FALSE,row.names=FALSE)
          scores.a<-read.table(file.a)
          score.b<-data[,2][assignment=="B"]
          score.b<-t(as.matrix(score.b))
          write.table(score.b,file=file.b,append=TRUE,col.names=FALSE,row.names=FALSE)
          scores.b<-read.table(file.b)
        }
        if(nrow(assignments)==quantity)break
      }
      scores.a<-as.matrix(scores.a)
      scores.b<-as.matrix(scores.b)
      if(statistic=="A-B"){
        distribution<-numeric(quantity)
        for(it in 1:quantity){
          distribution[it]<-mean(scores.a[it,])-mean(scores.b[it,])
        }
        distribution<-sort(distribution)
        observed.statistic<-mean(observed.a)-mean(observed.b)
      }
      if(statistic=="B-A"){
        distribution<-numeric(quantity)
        for(it in 1:quantity){
          distribution[it]<-mean(scores.b[it,])-mean(scores.a[it,])
        }
        distribution<-sort(distribution)
        observed.statistic<-mean(observed.b)-mean(observed.a)
      }
      if(statistic=="|A-B|"){
        distribution<-numeric(quantity)
        for(it in 1:quantity){
          distribution[it]<-abs(mean(scores.a[it,])-mean(scores.b[it,]))
        }
        distribution<-sort(distribution)
        observed.statistic<-abs(mean(observed.a)-mean(observed.b))
      }
      test<-distribution>=observed.statistic
      p.value<-sum(test)/quantity
      unlink(fileCRD,recursive=FALSE)
      unlink(file.a,recursive=FALSE)
      unlink(file.b,recursive=FALSE)
      if(save=="yes"|save=="check"){
        write.table(distribution,file=file,col.names=FALSE,row.names=FALSE)
      }
      return(p.value)
    }
  }
  
  if(design=="RBD"){
    observed.a<-data[,2][data[,1]=="A"]
    observed.b<-data[,2][data[,1]=="B"]
    MT<-nrow(data)
    quantity<-2^(MT/2)
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    observed1<-rbind(observed.a,observed.b)
    observed2<-rbind(observed.b,observed.a)
    scores.a<-numeric()
    for(it in 1:(MT/2)){
      scores.a<-cbind(scores.a,cbind(rep(cbind(rep(observed1[,it],rep(2^it/2,2))),2^(MT/2)/2^it)))
    }
    scores.b<-numeric()
    for(it in 1:(MT/2)){
      scores.b<-cbind(scores.b,cbind(rep(cbind(rep(observed2[,it],rep(2^it/2,2))),2^(MT/2)/2^it)))
    }
    if(statistic=="A-B"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-mean(scores.a[it,])-mean(scores.b[it,])
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.a)-mean(observed.b)
    }
    if(statistic=="B-A"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-mean(scores.b[it,])-mean(scores.a[it,])
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b)-mean(observed.a)
    }
    if(statistic=="|A-B|"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-abs(mean(scores.a[it,])-mean(scores.b[it,]))
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(mean(observed.a)-mean(observed.b))
    }
    test<-distribution>=observed.statistic
    p.value<-sum(test)/quantity
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE)
    }
    return(p.value)
  }
  
  if(design=="ATD"){
    observed.a<-data[,2][data[,1]=="A"]
    observed.b<-data[,2][data[,1]=="B"]
    MT<-nrow(data)
    quantityCRD<-choose(MT,MT/2)
    if(MT<=20){
      if(save=="yes"){
        file<-file.choose(new=FALSE)
      }
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
            sum.b2<-sum.b2 + dist.check.b[itr-itr2+1,]
          }
          sum.b<-rbind(sum.b,sum.b2)
        }
        sum.a.check<-sum.a==limit
        sum.b.check<-sum.b==limit
        sum.rows.a<-numeric(ncol(sum.a.check))
        for(it in 1:ncol(sum.a.check)){
          sum.rows.a[it]<-sum(sum.a.check[,it])
        }
        sum.rows.b<-numeric(ncol(sum.b.check))
        for(it in 1:ncol(sum.b.check)){
          sum.rows.b[it]<-sum(sum.b.check[,it])
        }
        check.stop<-sum.rows.a+sum.rows.b!=0
        stop<-order(check.stop)[1:sum(check.stop==F)]
      }
      indexes.a<-numeric()
      for(it in 1:length(stop)){
        indexes.a<-rbind(indexes.a,index.a[,stop[it]])
      }
      indexes.b<-numeric()
      for(it in 1:length(stop)){
        indexes.b<-rbind(indexes.b,index.b[,stop[it]])
      }
      scores.a<-numeric()
      for(it in 1:ncol(indexes.a)){
        scores.a<-cbind(scores.a,data[,2][indexes.a[,it]])
      }
      scores.b<-numeric()
      for(it in 1:ncol(indexes.b)){
        scores.b<-cbind(scores.b,data[,2][indexes.b[,it]])
      }
      quantity<-nrow(scores.a)
      if(statistic=="A-B"){
        distribution<-numeric(nrow(scores.a))
        for(it in 1:nrow(scores.a)){
          distribution[it]<-mean(scores.a[it,])-mean(scores.b[it,])
        }
        distribution<-sort(distribution)
        observed.statistic<-mean(observed.a)-mean(observed.b)
      }
      if(statistic=="B-A"){
        distribution<-numeric(nrow(scores.a))
        for(it in 1:nrow(scores.a)){
          distribution[it]<-mean(scores.b[it,])-mean(scores.a[it,])
        }
        distribution<-sort(distribution)
        observed.statistic<-mean(observed.b)-mean(observed.a)
      }
      if(statistic=="|A-B|"){
        distribution<-numeric(nrow(scores.a))
        for(it in 1:nrow(scores.a)){
          distribution[it]<-abs(mean(scores.a[it,])-mean(scores.b[it,]))
        }
        distribution<-sort(distribution)
        observed.statistic<-abs(mean(observed.a)-mean(observed.b))
      }
      test<-distribution>=observed.statistic
      p.value<-sum(test)/quantity
      if(save=="yes"|save=="check"){
        write.table(distribution,file=file,row.names=FALSE,col.names=FALSE)
      }
      return(p.value)
    }
    if(MT>20){
      fileCRD<-tempfile(pattern="CRDassignments",tmpdir=tempdir())
      file.a<-tempfile(pattern="ascores",tmpdir=tempdir())
      file.b<-tempfile(pattern="bscores",tmpdir=tempdir())
      if(save=="yes"){
        file<-file.choose(new=FALSE)
      }
      N<-c(rep(0,MT/2),rep(1,MT/2))
      assignment<-matrix(0,ncol=MT)
      assignment<-rbind(rep(c(0,1),MT/2))
      write.table(assignment,file=fileCRD,append=TRUE,col.names=FALSE,row.names=FALSE)
      CRD<-read.table(fileCRD)
      for(it in 1:(length(assignment))){
        if(assignment[,it]==0){
          assignment[,it]<-"A"
        }
        else{
          assignment[,it]<-"B"
        }
      }
      assignment<-as.vector(assignment)
      score.a<-data[,2][assignment=="A"]
      score.a<-t(as.matrix(score.a))
      write.table(score.a,file=file.a,append=TRUE,col.names=FALSE,row.names=FALSE)
      scores.a<-read.table(file.a)
      score.b<-data[,2][assignment=="B"]
      score.b<-t(as.matrix(score.b))
      write.table(score.b,file=file.b,append=TRUE,col.names=FALSE,row.names=FALSE)
      scores.b<-read.table(file.b)
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
            for(it in 1:(length(assignment))){
              if(assignment[,it]==0){
                assignment[,it]<-"A"
              }
              else{
                assignment[,it]<-"B"
              }
            }
            assignment<-as.vector(assignment)
            score.a<-data[,2][assignment=="A"]
            score.a<-t(as.matrix(score.a))
            write.table(score.a,file=file.a,append=TRUE,col.names=FALSE,row.names=FALSE)
            scores.a<-read.table(file.a)
            score.b<-data[,2][assignment=="B"]
            score.b<-t(as.matrix(score.b))
            write.table(score.b,file=file.b,append=TRUE,col.names=FALSE,row.names=FALSE)
            scores.b<-read.table(file.b)
          }
        }
        if(nrow(CRD)==quantityCRD)break
      }
      scores.a<-as.matrix(scores.a)
      scores.b<-as.matrix(scores.b)
      quantity<-nrow(scores.a)
      if(statistic=="A-B"){
        distribution<-numeric(nrow(scores.a))
        for(it in 1:nrow(scores.a)){
          distribution[it]<-mean(scores.a[it,])-mean(scores.b[it,])
        }
        distribution<-sort(distribution)
        observed.statistic<-mean(observed.a)-mean(observed.b)
      }
      if(statistic=="B-A"){
        distribution<-numeric(nrow(scores.a))
        for(it in 1:nrow(scores.a)){
          distribution[it]<-mean(scores.b[it,])-mean(scores.a[it,])
        }
        distribution<-sort(distribution)
        observed.statistic<-mean(observed.b)-mean(observed.a)
      }
      if(statistic=="|A-B|"){
        distribution<-numeric(nrow(scores.a))
        for(it in 1:nrow(scores.a)){
          distribution[it]<-abs(mean(scores.a[it,])-mean(scores.b[it,]))
        }
        distribution<-sort(distribution)
        observed.statistic<-abs(mean(observed.a)-mean(observed.b))
      }
      test<-distribution>=observed.statistic
      p.value<-sum(test)/quantity
      unlink(fileCRD,recursive=FALSE)
      unlink(file.a,recursive=FALSE)
      unlink(file.b,recursive=FALSE)
      if(save=="yes"|save=="check"){
        write.table(distribution,file=file,col.names=FALSE,row.names=FALSE)
      }
      return(p.value)
    }
  }
  
  if(design=="AB"){
    observed.a<-data[,2][data[,1]=="A"]
    observed.b<-data[,2][data[,1]=="B"]
    observed<-data[,2]
    MT<-nrow(data)
    quantity<-choose(MT-2*limit+1,1)    
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    index.a<-limit:(MT-limit)
    scores.a<-list()
    for(it in 1:quantity){
      scores.a[[it]]<-c(observed[1:index.a[it]])
    }
    scores.b<-list()
    for(it in 1:quantity){
      scores.b[[it]]<-c(observed[(index.a[it]+1):MT])
    }
    if(statistic=="A-B"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-mean(scores.a[[it]])-mean(scores.b[[it]])
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.a)-mean(observed.b)
    }
    if(statistic=="B-A"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-mean(scores.b[[it]])-mean(scores.a[[it]])
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b)-mean(observed.a)
    }
    if(statistic=="|A-B|"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-abs(mean(scores.a[[it]])-mean(scores.b[[it]]))
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(mean(observed.a)-mean(observed.b))
    }
    test<-distribution>=observed.statistic
    p.value<-sum(test)/quantity
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
      return(p.value)
    }
    if(save=="no"){
      return(p.value)
    }
  }
  
  if(design=="ABA"){
    observed.a1<-data[,2][data[,1]=="A1"]
    observed.b1<-data[,2][data[,1]=="B1"]
    observed.a2<-data[,2][data[,1]=="A2"]
    observed.a<-c(observed.a1,observed.a2)
    observed<-data[,2]
    MT<-nrow(data)
    quantity<-choose(MT-3*limit+2,2)    
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    index1<-1:(MT-3*limit+1)
    index2<-rev(index1)
    index.a<-numeric()
    for(it in 1:length(index1)){
      index.a<-c(index.a,(rep((limit-1+index1[it]),index2[it])))
    }
    scores.a1<-list()
    for(it in 1:quantity){
      scores.a1[[it]]<-c(observed[1:(index.a[it])])
    }
    index.b<-numeric() 
    for(itr in index1){
      for(it in itr:(MT-3*limit+1)){
        index.b<-c(index.b,2*limit-1+it)
      }
    }
    scores.b1<-list()
    for(it in 1:quantity){
      scores.b1[[it]]<-c(observed[(index.a[it]+1):(index.b[it])])
    }
    scores.a2<-list()
    for(it in 1:quantity){
      scores.a2[[it]]<-c(observed[(index.b[it]+1):(MT)])
    }
    scores.a<-list()
    for(it in 1:quantity){
      scores.a[[it]]<-c(scores.a1[[it]],scores.a2[[it]])
    }
    mean.a<-numeric(quantity)
    for(it in 1:quantity){
      mean.a[it]<-mean(scores.a[[it]])
    }
    mean.b<-numeric(quantity)
    for(it in 1:quantity){
      mean.b[it]<-mean(scores.b1[[it]])
    }
    pmean.a<-numeric(quantity)
    for(it in 1:quantity){
      pmean.a[it]<-(mean(scores.a1[[it]])+mean(scores.a2[[it]]))/2
    }
    mean.a1<-numeric(quantity)
    for(it in 1:quantity){
      mean.a1[it]<-mean(scores.a1[[it]])
    }
    mean.a2<-numeric(quantity)
    for(it in 1:quantity){
      mean.a2[it]<-mean(scores.a2[[it]])
    }
    if(statistic=="A-B"){	
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-mean.a[it]-mean.b[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.a)-mean(observed.b1)
    }
    if(statistic=="B-A"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-mean.b[it]-mean.a[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b1)-mean(observed.a)
    }
    if(statistic=="|A-B|"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-abs(mean.a[it]-mean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(mean(observed.a)-mean(observed.b1))
    }
    if(statistic=="PA-PB"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-pmean.a[it]-mean.b[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-((mean(observed.a1)+mean(observed.a2))/2)-mean(observed.b1)
    }
    if(statistic=="PB-PA"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-mean.b[it]-pmean.a[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b1)-((mean(observed.a1)+mean(observed.a2))/2)
    }
    if(statistic=="|PA-PB|"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-abs(pmean.a[it]-mean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(((mean(observed.a1)+mean(observed.a2))/2)-mean(observed.b1))
    }
    if(statistic=="AA-BB"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-(mean.a1[it]+mean.a2[it])-(mean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-(mean(observed.a1)+mean(observed.a2))-(mean(observed.b1))
    }
    if(statistic=="BB-AA"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-mean.b[it]-(mean.a1[it]+mean.a2[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-(mean(observed.b1))-(mean(observed.a1)+mean(observed.a2))
    }
    if(statistic=="|AA-BB|"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-abs((mean.a1[it]+mean.a2[it])-(mean.b[it]))
      }
      distribution<-sort(distribution)
      observed.statistic<-abs((mean(observed.a1)+mean(observed.a2))-(mean(observed.b1)))
    }
    test<-distribution>=observed.statistic
    p.value<-sum(test)/quantity
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
      return(p.value)
    }
    if(save=="no"){
      return(p.value)
    }
  }
  
  if(design=="ABAB"){
    observed.a1<-data[,2][data[,1]=="A1"]
    observed.b1<-data[,2][data[,1]=="B1"]
    observed.a2<-data[,2][data[,1]=="A2"]
    observed.b2<-data[,2][data[,1]=="B2"]
    observed.a<-c(observed.a1,observed.a2)
    observed.b<-c(observed.b1,observed.b2)
    observed<-data[,2]
    MT<-nrow(data)
    quantity<-choose(MT-4*limit+3,3)    
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    index1<-1:(MT-4*limit+1)
    index2<-rev(cumsum(index1))
    index.a1<-numeric()
    for(it in 1:length(index1)){
      index.a1<-c(index.a1,(rep((limit+index1[it]-1),index2[it])))
    }
    scores.a1<-list()
    for(it in 1:quantity){
      scores.a1[[it]]<-c(observed[1:(index.a1[it])])
    }
    index.b1<-numeric() 
    for(itr in index1){
      for(it in (itr-1):(MT-4*limit)){
        index.b1<-c(index.b1,rep((2*limit+it),(MT-4*limit+1-it)))
      }
    }
    scores.b1<-list()
    for(it in 1:quantity){
      scores.b1[[it]]<-c(observed[(1+index.a1[it]):index.b1[it]])
    }
    indexa2<-numeric()
    for(it in 1:length(index1)){
      indexa2<-c(indexa2,(index1[it]:length(index1)))
    }
    index.a2<-numeric()
    for(it in 1:length(indexa2)){
      index.a2<-c(index.a2,(4*limit-limit-1+(indexa2[it]:length(index1))))
    }
    scores.a2<-list()
    for(it in 1:quantity){
      scores.a2[[it]]<-c(observed[(1+index.b1[it]):index.a2[it]])
    }
    scores.b2<-list()
    for(it in 1:quantity){
      scores.b2[[it]]<-c(observed[(1+index.a2[it]):MT])
    }
    scores.a<-list()
    for(it in 1:quantity){
      scores.a[[it]]<-c(scores.a1[[it]],scores.a2[[it]])
    }
    scores.b<-list()
    for(it in 1:quantity){
      scores.b[[it]]<-c(scores.b1[[it]],scores.b2[[it]])
    }
    mean.a<-numeric(quantity)
    for(it in 1:quantity){
      mean.a[it]<-mean(scores.a[[it]])
    }
    mean.b<-numeric(quantity)
    for(it in 1:quantity){
      mean.b[it]<-mean(scores.b[[it]])
    }
    pmean.a<-numeric(quantity)
    for(it in 1:quantity){
      pmean.a[it]<-(mean(scores.a1[[it]])+mean(scores.a2[[it]]))/2
    }
    pmean.b<-numeric(quantity)
    for(it in 1:quantity){
      pmean.b[it]<-(mean(scores.b1[[it]])+mean(scores.b2[[it]]))/2
    }
    mean.a1<-numeric(quantity)
    for(it in 1:quantity){
      mean.a1[it]<-mean(scores.a1[[it]])
    }
    mean.a2<-numeric(quantity)
    for(it in 1:quantity){
      mean.a2[it]<-mean(scores.a2[[it]])
    }
    mean.b1<-numeric(quantity)
    for(it in 1:quantity){
      mean.b1[it]<-mean(scores.b1[[it]])
    }
    mean.b2<-numeric(quantity)
    for(it in 1:quantity){
      mean.b2[it]<-mean(scores.b2[[it]])
    }
    if(statistic=="A-B"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-mean.a[it]-mean.b[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.a)-mean(observed.b)
    }
    if(statistic=="B-A"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-mean.b[it]-mean.a[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b)-mean(observed.a)
    }
    if(statistic=="|A-B|"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-abs(mean.a[it]-mean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(mean(observed.a)-mean(observed.b))
    }
    if(statistic=="PA-PB"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-pmean.a[it]-pmean.b[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-((mean(observed.a1)+mean(observed.a2))/2)-((mean(observed.b1)+mean(observed.b2))/2)
    }
    if(statistic=="PB-PA"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-pmean.b[it]-pmean.a[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-((mean(observed.b1)+mean(observed.b2))/2)-((mean(observed.a1)+mean(observed.a2))/2)
    }
    if(statistic=="|PA-PB|"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-abs(pmean.a[it]-pmean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(((mean(observed.a1)+mean(observed.a2))/2)-((mean(observed.b1)+mean(observed.b2))/2))
    }
    if(statistic=="AA-BB"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-(mean.a1[it]+mean.a2[it])-(mean.b1[it]+mean.b2[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-(mean(observed.a1)+mean(observed.a2))-(mean(observed.b1)+mean(observed.b2))
    }
    if(statistic=="BB-AA"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-(mean.b1[it]+mean.b2[it])-(mean.a1[it]+mean.a2[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-(mean(observed.b1)+mean(observed.b2))-(mean(observed.a1)+mean(observed.a2))
    }
    if(statistic=="|AA-BB|"){
      distribution<-numeric(quantity)
      for(it in 1:quantity){
        distribution[it]<-abs((mean.a1[it]+mean.a2[it])-(mean.b1[it]+mean.b2[it]))
      }
      distribution<-sort(distribution)
      observed.statistic<-abs((mean(observed.a1)+mean(observed.a2))-(mean(observed.b1)+mean(observed.b2)))
    }
    test<-distribution>=observed.statistic
    p.value<-sum(test)/quantity
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
      return(p.value)
    }
    if(save=="no"){
      return(p.value)
    }
  }
  
  if(design=="MBD"){
    N<-ncol(data)/2
    MT<-nrow(data)
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
    coord<-list()
    for(itr in 1:length(number)){
      cor<-numeric()
      for(it in 1:number[itr]){
        cor<-c(cor,paste(itr,it,sep=""))
      }
      coord[[itr]]<-cor
    }
    startpt<-numeric(N)
    for(it in 1:N){
      if(number[it]!=1){
        startpt[it]<-sample(coord[[it]],1)
      }
      else{
        startpt[it]<-coord[[it]]
      }
    }
    fileSTARTPTS<-tempfile(pattern="startpoints",tmpdir=tempdir())
    startpt1<-rbind(startpt)
    write.table(startpt1,file=fileSTARTPTS,append=TRUE,col.names=FALSE,row.names=FALSE)
    startpts<-read.table(fileSTARTPTS)
    repeat{
      startpt<-numeric(N)
      for(it in 1:N){
        if(number[it]!=1){
          startpt[it]<-sample(coord[[it]],1)
        }
        else{
          startpt[it]<-coord[[it]]
        }
      }
      copy<-numeric()
      for(itr in 1:nrow(startpts)){
        copy2<-numeric(N)
        for(it in 1:N){
          copy2[it]<-startpt[it]==startpts[itr,it]
        }
        copy<-c(copy,prod(copy2))
      }
      if(sum(copy)==0){
        startpt1<-rbind(startpt)
        write.table(startpt1,file=fileSTARTPTS,append=TRUE,col.names=FALSE,row.names=FALSE)
        startpts<-read.table(fileSTARTPTS)
      }
      if(nrow(startpts)==prod(number))break
    }
    fileCOMBSTARTPOINTS<-tempfile(pattern="combstartpoints",tmpdir=tempdir())
    combstartpts1<-sample(startpts[1,],replace=FALSE)
    write.table(combstartpts1,file=fileCOMBSTARTPOINTS,append=TRUE,col.names=FALSE,row.names=FALSE)
    combstartpts<-read.table(fileCOMBSTARTPOINTS)
    for(iter in 1:nrow(startpts)){
      repeat{
        combstartpts1<-sample(startpts[iter,],replace=FALSE)
        copy<-numeric()
        for(itr in 1:nrow(combstartpts)){
          copy2<-numeric(N)
          for(it in 1:N){
            copy2[it]<-combstartpts1[it]==combstartpts[itr,it]
          }
          copy<-c(copy,prod(copy2))
        }
        if(sum(copy)==0){
          write.table(combstartpts1,file=fileCOMBSTARTPOINTS,append=TRUE,col.names=FALSE,row.names=FALSE)
          combstartpts<-read.table(fileCOMBSTARTPOINTS)
        }
        if(nrow(combstartpts)==iter*factorial(N))break
      }
    }
    for(itrow in 1:nrow(combstartpts)){
      for(itcol in 1:ncol(combstartpts)){
        for(it in 1:N){
          for(itr in 1:number[it]){
            if(combstartpts[itrow,itcol]==coord[[it]][itr]){combstartpts[itrow,itcol]<-limits[[it]][[1]][itr]}
          }
        }
      }
    }
    fileASSIGNMENTS<-tempfile(pattern="assignments",tmpdir=tempdir())
    assignment<-combstartpts[1,]
    write.table(assignment,file=fileASSIGNMENTS,append=TRUE,col.names=FALSE,row.names=FALSE)
    assignments<-read.table(fileASSIGNMENTS)
    for(iter in 2:nrow(combstartpts)){
      assignment<-numeric(N)
      for(it in 1:N){
        assignment[it]<-combstartpts[iter,it]
      }
      copy<-numeric()
      for(itr in 1:nrow(assignments)){
        copy2<-numeric(N)
        for(it in 1:N){
          copy2[it]<-assignment[it]==assignments[itr,it]
        }
        copy<-c(copy,prod(copy2))
      }
      if(sum(copy)==0){
        assignment1<-rbind(assignment)
        write.table(assignment1,file=fileASSIGNMENTS,append=TRUE,col.names=FALSE,row.names=FALSE)
        assignments<-read.table(fileASSIGNMENTS)
      }
    }
    observed.a<-list()
    for(it in 1:N){
      observed.a[[it]]<-data[,it*2][data[,(it*2)-1]=="A"]
    }
    observed.b<-list()
    for(it in 1:N){
      observed.b[[it]]<-data[,it*2][data[,(it*2)-1]=="B"]
    }
    if(statistic=="A-B"){
      differences<-numeric(N)
      for(it in 1:N){
        differences[it]<-mean(observed.a[[it]])-mean(observed.b[[it]])
      }
    }
    if(statistic=="B-A"){
      differences<-numeric(N)
      for(it in 1:N){
        differences[it]<-mean(observed.b[[it]])-mean(observed.a[[it]])
      }
    }
    if(statistic=="|A-B|"){
      differences<-numeric(N)
      for(it in 1:N){
        differences[it]<-abs(mean(observed.b[[it]])-mean(observed.a[[it]]))
      }
    }
    observed.statistic<-mean(differences)
    scores.a<-list()
    for(iter in 1:nrow(assignments)){
      ascores<-list()
      for(it in 1:N){
        ascores[[it]]<-data[1:(assignments[iter,it]-1),it*2]
      }
      scores.a[[iter]]<-ascores
    }
    scores.b<-list()
    for(iter in 1:nrow(assignments)){
      bscores<-list()
      for(it in 1:N){
        bscores[[it]]<-data[assignments[iter,it]:MT,it*2]
      }
      scores.b[[iter]]<-bscores
    }
    if(statistic=="A-B"){
      differs<-list()
      for(iter in 1:nrow(assignments)){
        differ<-numeric(N)
        for(it in 1:N){
          differ[it]<-mean(scores.a[[iter]][[it]])-mean(scores.b[[iter]][[it]])
        }
        differs[[iter]]<-differ
      }
    }
    if(statistic=="B-A"){
      differs<-list()
      for(iter in 1:nrow(assignments)){
        differ<-numeric(N)
        for(it in 1:N){
          differ[it]<-mean(scores.b[[iter]][[it]])-mean(scores.a[[iter]][[it]])
        }
        differs[[iter]]<-differ
      }
    }
    if(statistic=="|A-B|"){
      differs<-list()
      for(iter in 1:nrow(assignments)){
        differ<-numeric(N)
        for(it in 1:N){
          differ[it]<-abs(mean(scores.a[[iter]][[it]])-mean(scores.b[[iter]][[it]]))
        }
        differs[[iter]]<-differ
      }
    }
    distribution<-numeric(nrow(assignments))
    for(it in 1:nrow(assignments)){
      distribution[it]<-mean(differs[[it]])
    }
    distribution<-sort(distribution)
    test<-distribution>=observed.statistic
    p.value<-sum(test)/nrow(assignments)
    if(save=="yes"){
      fileSAVE<-file.choose(new=FALSE)
    }
    if(save=="yes"|save=="check"){
      write.table(distribution,file=fileSAVE,col.names=FALSE,row.names=FALSE,append=FALSE)
      return(p.value)
      unlink(fileSTARTPTS,recursive=FALSE)
      unlink(fileCOMBSTARTPOINTS,recursive=FALSE)
      unlink(fileASSIGNMENTS,recursive=FALSE)
    }
    if(save=="no"){
      return(p.value)
      unlink(fileSTARTPTS,recursive=FALSE)
      unlink(fileCOMBSTARTPOINTS,recursive=FALSE)
      unlink(fileASSIGNMENTS,recursive=FALSE)
    }
  }

}
