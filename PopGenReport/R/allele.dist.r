allele.dist<-function(population, mk.figures=TRUE){
  # package require adegenet and pegas
  if (class(population) != "genind") {
    message("You did not provide a valid genind object! Script stopped!")
    return
  }

  # initial steps...
  numloci<-length(locNames(population))  # this gets the total number of loci across all pops
  numpops<-length(popNames(population)) # this gets the total number of pops
  popnumallele<-population@loc.n.all     # this is a list of the population wide number of alleles at each pop
  lociname<-attributes(popnumallele)[[1]] # this is a list of the locinames (just L01, L02, L03,...)
  subdivpops<-seppop(population)

  # create list of matrices in which to place the numbers from summary
  alleletable<-vector("list",numloci)
  fralleletable<-vector("list",numloci)
  for(i in 1:numloci){
    alleletable[[i]]<-matrix(nrow=popnumallele[[i]],ncol=numpops)
    colnames(alleletable[[i]])<-popNames(population)
    rownames(alleletable[[i]])<-population@all.names[[i]]
    fralleletable[[i]]<-matrix(nrow=popnumallele[[i]],ncol=numpops)
    colnames(fralleletable[[i]])<-popNames(population)
    rownames(fralleletable[[i]])<-population@all.names[[i]]
  }

  # this is going to loop over all populations
  for (i in 1:numpops){
    x<-as.loci(subdivpops[[i]])
    s<-summary(x)
    # this loops over the loci
    for (j in 1:numloci){
    # this is the number of 
      namevec<-(names(s[[j]]$allele))
      numnames<-length(namevec)
      # j<-2
      tablenames<-(rownames(alleletable[[j]]))
      for (k in 1:numnames){
        rownum<-which(tablenames==namevec[k])
        #  message("i = ",i," j = ",j," k = ",k," rownum = ",rownum)
        alleletable[[j]][rownum,i]<-s[[j]]$allele[k]
      }  
    }
  }

  allpops<-as.loci(population)
  numbers<-summary(allpops)
  checkcnts<-matrix(nrow=numloci,ncol=2)
  for (i in 1:numloci){
    checkcnts[i,1]<-sum(numbers[[i]]$allele)
    checkcnts[i,2]<-sum(alleletable[[i]],na.rm=TRUE)
  }

  for (i in 1:numloci){
    for (j in 1:numpops){
      colsum<-sum(alleletable[[i]][,j],na.rm=TRUE)
      fralleletable[[i]][,j]<-round(alleletable[[i]][,j]/colsum, digits=3)
    }
  }
  
  if (mk.figures){
    breaks<-seq(0,1,0.05)
    color.palette  <- colorRampPalette(c("yellow", "red"))(length(breaks) - 1)
    for (i in 1:numloci){
      if(unname(population@loc.n.all[i])>1){
        figlabel<-paste("Loci: ",locNames(population)[i]," List # ",i,sep="")
        dat <- t(fralleletable[[i]])
        dat <- dat[,seq(ncol(dat),1,-1)]
        image( dat, col=color.palette, axes=FALSE, main=figlabel, zlim=c(0,1))
        rn <- rownames(dat)
        cn <- colnames(dat)  
        axis(1, at = seq(0,1,len=nrow(dat)),labels=rn, cex.axis= max(1-nrow(dat)/100,0.5), las=2 )
        axis(2, at = seq(0,1,len=ncol(dat)),labels=cn , las=2, cex.axis= max(1-ncol(dat)/100,0.5))
        box()
        co <- expand.grid(seq(0,1,len=nrow(dat)),seq(0,1,len=ncol(dat)))
        text(co[,1], co[,2],round(dat*100), cex=max(0.5,min(1-nrow(dat)/100, 1-ncol(dat)/100)))
      } else {
        message("Locus ",unname(locNames(population)[i])," has only ",unname(population@loc.n.all[i])," allele, figure not made \n")
      }
    }
  }
  
  # Find private alleles
  locus.private<-list()
  for (i in 1:length(alleletable)){
    counter<-0
    tmpmatrix<-matrix(NA,nrow=dim(alleletable[[i]])[1],ncol=2)
    colnames(tmpmatrix)<-c("Population","Allele")
    for (j in 1:dim(alleletable[[i]])[1]){
      cntpopsallele<-which(alleletable[[i]][j,]>0)
      if(length(cntpopsallele)==1){
        counter<-counter+1
        tmpmatrix[counter,1]<-names(cntpopsallele)
        tmpmatrix[counter,2]<-unname(rownames(alleletable[[i]]))[j]
      }
    }
    if(counter>0){
      locus.private[[i]]<-tmpmatrix[1:counter,]
    } else {
      locus.private[[i]]<-NA
    }
  }
  names(locus.private)<-locNames(population)
  
  alleletables<-list(count=alleletable, frequency=fralleletable, private.alleles=locus.private)
  return(alleletables)
}