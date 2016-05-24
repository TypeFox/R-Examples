assignments<-function(design,save="no",MT,limit,starts=file.choose(new=FALSE)){

  if(design=="CRD"){
    quantity<-choose(MT,MT/2)
    if(MT<=22){
      options(max.print=999999999)
      if(save=="yes"){
        file<-file.choose(new=FALSE)
      }
      index<-1:MT
      index.a<-combn(index,MT/2)
      assignments<-matrix("B",quantity,MT)
      for(itr in 1:quantity){
        for(it in 1:(MT/2)){
          assignments[itr,index.a[it,itr]]<-"A"
        }
      }
      if(save=="no"){
        return(assignments)
      }
      if(save=="yes"|save=="check"){
        write.table(assignments,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
        return(assignments)
      }
    }
    if(MT>22){
      if(save=="yes"){
        file<-file.choose(new=FALSE)
      }
      if(save=="no"){
        file<-tempfile(pattern="CRDassignments",tmpdir=tempdir())
      }
      N<-c(rep("A",MT/2),rep("B",MT/2))
      assignment<-matrix(0,ncol=MT)
      assignment<-rbind(sample(N,MT,replace=FALSE))
      write.table(assignment,file=file,append=TRUE,col.names=FALSE,row.names=FALSE)
      assignments<-read.table(file)
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
          write.table(assignment,file=file,append=TRUE,col.names=FALSE,row.names=FALSE)
          assignments<-read.table(file)
          if(nrow(assignments)==quantity)break
        }
      }
      return(assignments)
      if(save=="no"){
        unlink(file,recursive=FALSE)
      }
    }
  }
  
  if(design=="RBD"){
    options(max.print=999999999)
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    Na<-rep("A",MT/2)
    Nb<-rep("B",MT/2)
    Nab<-rbind(Na,Nb)
    Nba<-rbind(Nb,Na)
    assignment.odd<-numeric()
    for(it in 1:(MT/2)){
      assignment.odd<-cbind(assignment.odd,cbind(rep(cbind(rep(Nab[,it],rep(2^it/2,2))),2^(MT/2)/2^it)))
    }
    assignment.even<-numeric()
    for(it in 1:(MT/2)){
      assignment.even<-cbind(assignment.even,cbind(rep(cbind(rep(Nba[,it],rep(2^it/2,2))),2^(MT/2)/2^it)))
    }
    assignments<-numeric()
    for(it in 1:(MT/2)){
      assignments<-cbind(assignments,assignment.odd[,it],assignment.even[,it])
    }
    if(save=="no"){
      return(assignments)
    }
    if(save=="yes"|save=="check"){
      write.table(assignments,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
      return(assignments)
    }
  }
  
  if(design=="ATD"){
    quantityCRD<-choose(MT,MT/2)
    if(MT<=20){
      options(max.print=999999999)
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
      indexes.a<-numeric()
      for(it in 1:length(stop)){
        indexes.a<-rbind(indexes.a,index.a[,stop[it]])
      }
      indexes.b<-numeric()
      for(it in 1:length(stop)){
        indexes.b<-rbind(indexes.b,index.b[,stop[it]])
      }
      assignments<-matrix(0,nrow(indexes.a),MT)
      for(itr in 1:nrow(indexes.a)){
        for(it in 1:ncol(indexes.a)){
          assignments[itr,indexes.a[itr,it]]<-"A"
          assignments[itr,indexes.b[itr,it]]<-"B"
        }
      }
      if(save=="no"){
        return(assignments)
      }
      if(save=="yes"|save=="check"){
        write.table(assignments,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
        return(assignments)
      }
    }
    if(MT>20){
      if(save=="yes"){
        fileATD<-file.choose(new=FALSE)
      }
      if(save=="no"){
        fileATD<-tempfile(pattern="ATDassignments",tmpdir=tempdir())
      }
      fileCRD<-tempfile(pattern="CRDassignments",tmpdir=tempdir())
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
            for(it in 1:(length(assignment))){
              if(assignment[,it]==0){
                assignment[,it]<-"A"
              }
              else{
                assignment[,it]<-"B"
              }
            }
            write.table(assignment,file=fileATD,append=TRUE,col.names=FALSE,row.names=FALSE)
            assignments<-read.table(fileATD)
          }
        }
        if(nrow(CRD)==quantityCRD)break
      }
      return(assignments)
      unlink(fileCRD,recursive=FALSE)
      if(save=="no"){
        unlink(fileATD,recursive=FALSE)
      }
    }
  }
  
  if(design=="AB"){
    options(max.print=999999999)
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    quantity<-choose(MT-2*limit+1,1)
    assignments<-matrix("A",quantity,MT)
    index.b<-(limit+1):(MT-(limit-1))
    for(it in 1:quantity){
      assignments[it,index.b[it]:MT]<-"B"
    }
    if(save=="no"){
      return(assignments)
    }
    if(save=="yes"|save=="check"){
      write.table(assignments,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
      return(assignments)
    }
  }
  
  if(design=="ABA"){
    options(max.print=999999999)
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    quantity<-choose(MT-3*limit+2,2)
    assignments<-matrix("A",quantity,MT)
    index1<-1:(MT-3*limit+1)
    index2<-rev(index1)
    index.b.1<-numeric()
    for(it in 1:length(index1)){
      index.b.1<-c(index.b.1,rep(index1[it],index2[it]))
    }
    index.b.2<-numeric()
    for(itr in index1){
      for(it in itr:(MT-3*limit+1)){
        index.b.2<-c(index.b.2,2*limit-1+it)
      }
    }
    for(it in 1:quantity){
      assignments[it,(limit+index.b.1[it]):(index.b.2[it])]<-"B"
    }
    if(save=="no"){
      return(assignments)
    }
    if(save=="yes"|save=="check"){
      write.table(assignments,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
      return(assignments)
    }
  }
  
  if(design=="ABAB"){
    options(max.print=999999999)
    memory.limit(4095)
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    quantity<-choose(MT-4*limit+3,3)
    assignments<-matrix("A",quantity,MT)
    index1<-1:(MT-4*limit+1)
    index2<-rev(cumsum(index1))
    index.b1.1<-numeric()
    for(it in 1:length(index1)){
      index.b1.1<-c(index.b1.1,(rep((limit+index1[it]),(index2[it]))))
    }
    index.b1.2<-numeric()
    for(itr in index1){
      for(it in (itr-1):(MT-4*limit)){
        index.b1.2<-c(index.b1.2,rep((2*limit+it),(MT-4*limit+1-it)))
      }
    }
    for(it in 1:quantity){
      assignments[it,(index.b1.1[it]:index.b1.2[it])]<-"B"
    }
    indexb2<-numeric()
    for(it in 1:length(index1)){
      indexb2<-c(indexb2,index1[it:length(index1)])
    }
    index.b2<-numeric()
    for(it in 1:length(indexb2)){
      index.b2<-c(index.b2,indexb2[it]:length(index1))
    }
    for(it in 1:quantity){
      assignments[it,(4*limit-limit+index.b2[it]):MT]<-"B"
    }
    if(save=="no"){
      return(assignments)
    }
    if(save=="yes"|save=="check"){
      write.table(assignments,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
      return(assignments)
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
    if(save=="yes"){
      fileASSIGNMENTS<-file.choose(new=FALSE)
    }
    if(save=="no"){
      fileASSIGNMENTS<-tempfile(pattern="assignments",tmpdir=tempdir())
    }
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
    return(assignments)
    if(save=="yes"|save=="check"){
      unlink(fileSTARTPTS,recursive=FALSE)
      unlink(fileCOMBSTARTPOINTS,recursive=FALSE)
    }
    if(save=="no"){
      unlink(fileSTARTPTS,recursive=FALSE)
      unlink(fileCOMBSTARTPOINTS,recursive=FALSE)
      unlink(fileASSIGNMENTS,recursive=FALSE)
    }
  }

}
