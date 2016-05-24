pvalue.random<-function(design,statistic,save="no",number,limit,data=read.table(file.choose(new=FALSE)),starts=file.choose(new=FALSE)){

obs <- data[,2]

  if(design=="CRD"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    file.a<-tempfile(pattern="ascores",tmpdir=tempdir())
    file.b<-tempfile(pattern="bscores",tmpdir=tempdir())
    observed.a<-data[,2][data[,1]=="A"]
    observed.b<-data[,2][data[,1]=="B"]
    repeat{
      index<-sample(data[,1],replace=FALSE)
      scores.a<-data[,2][index=="A"]
      scores.a<-as.matrix(scores.a)
      scores.a<-t(scores.a)
      write.table(scores.a,file=file.a,append=TRUE,col.names=FALSE,row.names=FALSE)
      ascores<-read.table(file.a)
      scores.b<-data[,2][index=="B"]
      scores.b<-as.matrix(scores.b)
      scores.b<-t(scores.b)
      write.table(scores.b,file=file.b,append=TRUE,col.names=FALSE,row.names=FALSE)
      bscores<-read.table(file.b)
      if(nrow(bscores)==number)break
    }
    ascores<-as.matrix(ascores)
    mean.a<-numeric(number)
    for(it in 1:number){
      mean.a[it]<-mean(ascores[it,])
    }
    bscores<-as.matrix(bscores)
    mean.b<-numeric(number)
    for(it in 1:number){
      mean.b[it]<-mean(bscores[it,])
    }
    if(statistic=="A-B"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.a[it]-mean.b[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.a)-mean(observed.b)
    }
    if(statistic=="B-A"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.b[it]-mean.a[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b)-mean(observed.a)
    }
    if(statistic=="|A-B|"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-abs(mean.a[it]-mean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(mean(observed.a)-mean(observed.b))
    }
    test<-distribution>=observed.statistic
    p.value<-sum(test)/number
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE)
    }
    return(p.value)
    unlink(file.a,recursive=FALSE)
    unlink(file.b,recursive=FALSE)
  }
  
  if(design=="RBD"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    file.a<-tempfile(pattern="ascores",tmpdir=tempdir())
    file.b<-tempfile(pattern="bscores",tmpdir=tempdir())
    observed.a<-data[,2][data[,1]=="A"]
    observed.b<-data[,2][data[,1]=="B"]
    MT<-nrow(data)
    ab<-c("A","B")
    repeat{
      index<-numeric()
      repeat{
        index<-c(index,sample(ab,2,replace=FALSE))
        if(length(index)==MT)break
      }
      scores.a<-data[,2][index=="A"]
      scores.a<-as.matrix(scores.a)
      scores.a<-t(scores.a)
      write.table(scores.a,file=file.a,append=TRUE,col.names=FALSE,row.names=FALSE)
      ascores<-read.table(file.a)
      scores.b<-data[,2][index=="B"]
      scores.b<-as.matrix(scores.b)
      scores.b<-t(scores.b)
      write.table(scores.b,file=file.b,append=TRUE,col.names=FALSE,row.names=FALSE)
      bscores<-read.table(file.b)
      if(nrow(bscores)==number)break
    }
    ascores<-as.matrix(ascores)
    mean.a<-numeric(number)
    for(it in 1:number){
      mean.a[it]<-mean(ascores[it,])
    }
    bscores<-as.matrix(bscores)
    mean.b<-numeric(number)
    for(it in 1:number){
      mean.b[it]<-mean(bscores[it,])
    }
    if(statistic=="A-B"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.a[it]-mean.b[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.a)-mean(observed.b)
    }
    if(statistic=="B-A"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.b[it]-mean.a[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b)-mean(observed.a)
    }
    if(statistic=="|A-B|"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-abs(mean.a[it]-mean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(mean(observed.a)-mean(observed.b))
    }
    test<-distribution>=observed.statistic
    p.value<-sum(test)/number
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE)
    }
    return(p.value)
    unlink(file.a,recursive=FALSE)
    unlink(file.b,recursive=FALSE)
  }
  
  if(design=="ATD"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    file.a<-tempfile(pattern="ascores",tmpdir=tempdir())
    file.b<-tempfile(pattern="bscores",tmpdir=tempdir())
    observed.a<-data[,2][data[,1]=="A"]
    observed.b<-data[,2][data[,1]=="B"]
    MT<-nrow(data)
    N<-c(rep(0,MT/2),rep(1,MT/2))
    repeat{
      index<-matrix(0,ncol=MT)
      index<-rbind(sample(N,MT,replace=FALSE))
      check<-numeric()
      for(itr in 1:(MT-limit)){
        check2<-0
        for(it in itr:(itr+limit)){
          check2<-check2+index[,it]
        }
        check<-cbind(check,check2)
      }
      if(sum(check==(limit+1)|check==0)==0){
        for(it in 1:(length(index))){
        if(index[,it]==0){
          index[,it]<-"A"
        }
        else{
          index[,it]<-"B"
        }
      }
      scores.a<-data[,2][index=="A"]
      scores.b<-data[,2][index=="B"]
      scores.a<-as.matrix(scores.a)
      scores.a<-t(scores.a)
      write.table(scores.a,file=file.a,append=TRUE,col.names=FALSE,row.names=FALSE)
      ascores<-read.table(file.a)
      scores.b<-as.matrix(scores.b)
      scores.b<-t(scores.b)
      write.table(scores.b,file=file.b,append=TRUE,col.names=FALSE,row.names=FALSE)
      bscores<-read.table(file.b)
      if(nrow(bscores)==number)break
      }
    }
    ascores<-as.matrix(ascores)
    mean.a<-numeric(number)
    for(it in 1:number){
      mean.a[it]<-mean(ascores[it,])
    }
    bscores<-as.matrix(bscores)
    mean.b<-numeric(number)
    for(it in 1:number){
      mean.b[it]<-mean(bscores[it,])
    }
    if(statistic=="A-B"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.a[it]-mean.b[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.a)-mean(observed.b)
    }
    if(statistic=="B-A"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.b[it]-mean.a[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b)-mean(observed.a)
    }
    if(statistic=="|A-B|"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-abs(mean.a[it]-mean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(mean(observed.a)-mean(observed.b))
    }
    test<-distribution>=observed.statistic
    p.value<-sum(test)/number
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE)
    }
    return(p.value)
    unlink(file.a,recursive=FALSE)
    unlink(file.b,recursive=FALSE)
  }
  
  if(design=="AB"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    observed.a<-data[,2][data[,1]=="A"]
    observed.b<-data[,2][data[,1]=="B"]
    observed<-data[,2]
    MT<-nrow(data)
    quantity<-choose(MT-2*limit+1,1)
    selection<-sample(1:quantity,number,replace=TRUE)
    index.a<-limit:(MT-limit)
    scores.a<-list()
    for(it in 1:number){
      scores.a[[it]]<-c(observed[1:index.a[selection[it]]])
    }
    scores.b<-list()
    for(it in 1:number){
      scores.b[[it]]<-c(observed[(1+index.a[selection[it]]):length(observed)])
    }
    if(statistic=="A-B"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean(scores.a[[it]])-mean(scores.b[[it]])
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.a)-mean(observed.b)
    }
    if(statistic=="B-A"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean(scores.b[[it]])-mean(scores.a[[it]])
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b)-mean(observed.a)
    }
    if(statistic=="|A-B|"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-abs(mean(scores.a[[it]])-mean(scores.b[[it]]))
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(mean(observed.a)-mean(observed.b))
    }
    test<-distribution>=observed.statistic
    p.value<-sum(test)/number
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
      return(p.value)
    }
    if(save=="no"){
      return(p.value)
    }
  }
  
  if(design=="ABA"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    observed.a1<-data[,2][data[,1]=="A1"]
    observed.b1<-data[,2][data[,1]=="B1"]
    observed.a2<-data[,2][data[,1]=="A2"]
    observed.a<-c(observed.a1,observed.a2)
    observed<-data[,2]
    MT<-nrow(data)
    quantity<-choose(MT-3*limit+2,2)
    selection<-sample(1:quantity,number,replace=TRUE)
    index1<-1:(MT-3*limit+1)
    index2<-rev(index1)
    index.a<-numeric()
    for(it in 1:length(index1)){
      index.a<-c(index.a,(rep((limit-1+index1[it]),index2[it])))
    }
    index.b<-numeric()
    for(itr in index1){
      for(it in itr:(MT-3*limit+1)){
        index.b<-c(index.b,2*limit-1+it)
      }
    }
    scores.a1<-list()
    for(it in 1:number){
      scores.a1[[it]]<-c(observed[1:(index.a[selection[it]])])
    }
    scores.b1<-list()
    for(it in 1:number){
      scores.b1[[it]]<-c(observed[(1+index.a[selection[it]]):(index.b[selection[it]])])
    }
    scores.a2<-list()
    for(it in 1:number){
      scores.a2[[it]]<-c(observed[(1+index.b[selection[it]]):(MT)])
    }
    scores.a<-list()
    for(it in 1:number){
      scores.a[[it]]<-c(scores.a1[[it]],scores.a2[[it]])
    }
    mean.a<-numeric(number)
    for(it in 1:number){
      mean.a[it]<-mean(scores.a[[it]])
    }
    mean.b<-numeric(number)
    for(it in 1:number){
      mean.b[it]<-mean(scores.b1[[it]])
    }
    pmean.a<-numeric(number)
    for(it in 1:number){
      pmean.a[it]<-(mean(scores.a1[[it]])+mean(scores.a2[[it]]))/2
    }
    mean.a1<-numeric(number)
    for(it in 1:number){
      mean.a1[it]<-mean(scores.a1[[it]])
    }
    mean.a2<-numeric(number)
    for(it in 1:number){
      mean.a2[it]<-mean(scores.a2[[it]])
    }
    if(statistic=="A-B"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.a[it]-mean.b[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.a)-mean(observed.b1)
    }
    if(statistic=="B-A"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.b[it]-mean.a[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b1)-mean(observed.a)
    }
    if(statistic=="|A-B|"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-abs(mean.a[it]-mean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(mean(observed.a)-mean(observed.b1))
    }
    if(statistic=="PA-PB"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-pmean.a[it]-mean.b[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-((mean(observed.a1)+mean(observed.a2))/2)-mean(observed.b1)
    }
    if(statistic=="PB-PA"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.b[it]-pmean.a[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b1)-((mean(observed.a1)+mean(observed.a2))/2)
    }
    if(statistic=="|PA-PB|"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-abs(pmean.a[it]-mean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(((mean(observed.a1)+mean(observed.a2))/2)-mean(observed.b1))
    }
    if(statistic=="AA-BB"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-(mean.a1[it]+mean.a2[it])-(mean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-(mean(observed.a1)+mean(observed.a2))-(mean(observed.b1))
    }
    if(statistic=="BB-AA"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.b[it]-(mean.a1[it]+mean.a2[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-(mean(observed.b1))-(mean(observed.a1)+mean(observed.a2))
    }
    if(statistic=="|AA-BB|"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-abs((mean.a1[it]+mean.a2[it])-(mean.b[it]))
      }
      distribution<-sort(distribution)
      observed.statistic<-abs((mean(observed.a1)+mean(observed.a2))-(mean(observed.b1)))
    }
    test<-distribution>=observed.statistic
    p.value<-sum(test)/number
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
      return(p.value)
    }
    if(save=="no"){
      return(p.value)
    }
  }
  
  if(design=="ABAB"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    observed.a1<-data[,2][data[,1]=="A1"]
    observed.b1<-data[,2][data[,1]=="B1"]
    observed.a2<-data[,2][data[,1]=="A2"]
    observed.b2<-data[,2][data[,1]=="B2"]
    observed.a<-c(observed.a1,observed.a2)
    observed.b<-c(observed.b1,observed.b2)
    observed<-data[,2]
    MT<-nrow(data)
    quantity<-choose(MT-4*limit+3,3)
    selection<-sample(1:quantity,number,replace=TRUE)
    index1<-1:(MT-4*limit+1)
    index2<-rev(cumsum(index1))
    index.a1<-numeric()
    for(it in 1:length(index1)){
      index.a1<-c(index.a1,(rep((limit+index1[it]-1),index2[it])))
    }
    scores.a1<-list()
    for(it in 1:number){
      scores.a1[[it]]<-c(observed[1:(index.a1[selection[it]])])
    }
    index.b1<-numeric()
    for(itr in index1){
      for(it in (itr-1):(MT-4*limit)){
        index.b1<-c(index.b1,rep((2*limit+it),(MT-4*limit+1-it)))
      }
    }
    scores.b1<-list()
    for(it in 1:number){
      scores.b1[[it]]<-c(observed[(1+index.a1[selection[it]]):index.b1[selection[it]]])
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
    for(it in 1:number){
      scores.a2[[it]]<-c(observed[(1+index.b1[selection[it]]):index.a2[selection[it]]])
    }
    scores.b2<-list()
    for(it in 1:number){
      scores.b2[[it]]<-c(observed[(1+index.a2[selection[it]]):MT])
    }
    scores.a<-list()
    for(it in 1:number){
      scores.a[[it]]<-c(scores.a1[[it]],scores.a2[[it]])
    }
    scores.b<-list()
    for(it in 1:number){
      scores.b[[it]]<-c(scores.b1[[it]],scores.b2[[it]])
    }
    mean.a<-numeric(number)
    for(it in 1:number){
      mean.a[it]<-mean(scores.a[[it]])
    }
    mean.b<-numeric(number)
    for(it in 1:number){
      mean.b[it]<-mean(scores.b[[it]])
    }	
    pmean.a<-numeric(number)
    for(it in 1:number){
      pmean.a[it]<-(mean(scores.a1[[it]])+mean(scores.a2[[it]]))/2
    }
    pmean.b<-numeric(number)
    for(it in 1:number){
      pmean.b[it]<-(mean(scores.b1[[it]])+mean(scores.b2[[it]]))/2
    }
    mean.a1<-numeric(number)
    for(it in 1:number){
      mean.a1[it]<-mean(scores.a1[[it]])
    }
    mean.a2<-numeric(number)
    for(it in 1:number){
      mean.a2[it]<-mean(scores.a2[[it]])
    }
    mean.b1<-numeric(number)
    for(it in 1:number){
      mean.b1[it]<-mean(scores.b1[[it]])
    }
    mean.b2<-numeric(number)
    for(it in 1:number){
      mean.b2[it]<-mean(scores.b2[[it]])
    }
    if(statistic=="A-B"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.a[it]-mean.b[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.a)-mean(observed.b)
    }
    if(statistic=="B-A"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-mean.b[it]-mean.a[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-mean(observed.b)-mean(observed.a)
    }
    if(statistic=="|A-B|"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-abs(mean.a[it]-mean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(mean(observed.a)-mean(observed.b))
    }
    if(statistic=="PA-PB"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-pmean.a[it]-pmean.b[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-((mean(observed.a1)+mean(observed.a2))/2)-((mean(observed.b1)+mean(observed.b2))/2)
    }
    if(statistic=="PB-PA"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-pmean.b[it]-pmean.a[it]
      }
      distribution<-sort(distribution)
      observed.statistic<-((mean(observed.b1)+mean(observed.b2))/2)-((mean(observed.a1)+mean(observed.a2))/2)
    }
    if(statistic=="|PA-PB|"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-abs(pmean.a[it]-pmean.b[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-abs(((mean(observed.a1)+mean(observed.a2))/2)-((mean(observed.b1)+mean(observed.b2))/2))
    }
    if(statistic=="AA-BB"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-(mean.a1[it]+mean.a2[it])-(mean.b1[it]+mean.b2[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-(mean(observed.a1)+mean(observed.a2))-(mean(observed.b1)+mean(observed.b2))
    }
    if(statistic=="BB-AA"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-(mean.b1[it]+mean.b2[it])-(mean.a1[it]+mean.a2[it])
      }
      distribution<-sort(distribution)
      observed.statistic<-(mean(observed.b1)+mean(observed.b2))-(mean(observed.a1)+mean(observed.a2))
    }
    if(statistic=="|AA-BB|"){
      distribution<-numeric(number)
      for(it in 1:number){
        distribution[it]<-abs((mean.a1[it]+mean.a2[it])-(mean.b1[it]+mean.b2[it]))
      }
      distribution<-sort(distribution)
      observed.statistic<-abs((mean(observed.a1)+mean(observed.a2))-(mean(observed.b1)+mean(observed.b2)))
    }
    test<-distribution>=observed.statistic
    p.value<-sum(test)/number
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
    numbers<-numeric(N)
    for(it in 1:N){
      numbers[it]<-length(limits[[it]][[1]])
    }
    fileCOMBSTARTPTS<-tempfile(pattern="startpoints",tmpdir=tempdir())
    repeat{
      startpt<-numeric(N)
      for(it in 1:N){
        if(numbers[it]!=1){
          startpt[it]<-sample(limits[[it]][[1]],1)
        }
        else{
          startpt[it]<-limits[[it]][[1]]
        }
      }
      selectdesign<-sample(startpt,replace=FALSE)
      selectdesign1<-rbind(selectdesign)
      write.table(selectdesign1,file=fileCOMBSTARTPTS,append=TRUE,col.names=FALSE,row.names=FALSE)
      combstartpts<-read.table(fileCOMBSTARTPTS)
      if(nrow(combstartpts)==number)break
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
    for(iter in 1:number){
      ascores<-list()
      for(it in 1:N){
        ascores[[it]]<-data[1:(combstartpts[iter,it]-1),it*2]
      }
      scores.a[[iter]]<-ascores
    }
    scores.b<-list()
    for(iter in 1:number){
      bscores<-list()
      for(it in 1:N){
        bscores[[it]]<-data[combstartpts[iter,it]:MT,it*2]
      }
      scores.b[[iter]]<-bscores
    }
    if(statistic=="A-B"){
      differs<-list()
      for(iter in 1:number){
        differ<-numeric(N)
        for(it in 1:N){
          differ[it]<-mean(scores.a[[iter]][[it]])-mean(scores.b[[iter]][[it]])
        }
        differs[[iter]]<-differ
      }
    }
    if(statistic=="B-A"){
      differs<-list()
      for(iter in 1:number){
        differ<-numeric(N)
        for(it in 1:N){
          differ[it]<-mean(scores.b[[iter]][[it]])-mean(scores.a[[iter]][[it]])
        }
        differs[[iter]]<-differ
      }
    }
    if(statistic=="|A-B|"){
      differs<-list()
      for(iter in 1:number){
        differ<-numeric(N)
        for(it in 1:N){
          differ[it]<-abs(mean(scores.a[[iter]][[it]])-mean(scores.b[[iter]][[it]]))
        }
        differs[[iter]]<-differ
      }
    }
    distribution<-numeric(number)
    for(it in 1:number){
      distribution[it]<-mean(differs[[it]])
    }
    distribution<-sort(distribution)
    test<-distribution>=observed.statistic
    p.value<-sum(test)/number
    if(save=="yes"){
      fileSAVE<-file.choose(new=FALSE)
    }
    if(save=="yes"|save=="check"){
      write.table(distribution,file=fileSAVE,col.names=FALSE,row.names=FALSE,append=FALSE)
      return(p.value)
      unlink(fileCOMBSTARTPTS,recursive=FALSE)
    }
    if(save=="no"){
      return(p.value)
      unlink(fileCOMBSTARTPTS,recursive=FALSE)
    }
  }

}
