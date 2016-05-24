"legalG"<-function(G, A, ped, time_born=NULL, marker.type="MSW"){

  nind<-length(ped[,1])
  nloci<-length(A)
  nall<-unlist(lapply(A, length))
  maxall<-max(nall)
  namesG<-names(G)

  if(is.genotypeD(G[[1]])){
    for(i in 1:length(A)){  
      G[[i]]<-as.matrix(cbind(as.matrix(G[[i]]), as.matrix(G[[i]])))
      ones<-which(G[[i]][,1]=="1")
      G[[i]][,1][ones]<-as.character(rbinom(length(ones), 1, A[[i]][2]))
      G[[i]]<-as.genotype(G[[i]], alleles=c("0","1"))
    }
  }

  #################################### order pedigree and data #############################################

  if(all(is.na(ped))){
    ped[,1]<-1:dim(ped)[1]
  }

  oped<-orderPed(ped, time_born=NULL)

  rearrange_data<-match(oped[,1], ped[,1])
  dam<-match(oped[,2], oped[,1])-1
  dam[which(is.na(dam)==T)]<-nind
  sire<-match(oped[,3], oped[,1])-1
  sire[which(is.na(sire)==T)]<-nind

  G<-lapply(G, function(x){x[rearrange_data]})

  ###########################################################################################################

  legal<-TRUE

  mtype.numeric<-sum(c("MSC", "AFLP", "MSW", "SNP")%in%marker.type*c(1:4))

  output<-.C("legalG",
	as.integer(nind),		 
	as.integer(dam),		
	as.integer(sire),		
	as.integer(nloci),		
	as.integer(nall),		
	as.integer(maxall),		
        as.double(unlist(A)),                   
        as.integer(GtoC(G, (marker.type!="MSC" & marker.type!="MSW"))),                  
        as.logical(legal),
        as.integer(mtype.numeric))


  tmp<-array(output[[8]], c(1+(marker.type=="MSC" | marker.type=="MSW"), length(ped[,1]), length(A)))+as.numeric(marker.type=="MSC" | marker.type=="MSW")

  Gnew<-as.data.frame(matrix(NA, length(ped[,1]), 2*length(A)))

  if(marker.type=="MSC" | marker.type=="MSW"){
    for(i in 1:length(A)){
      Gnew[,c(((i*2)-1):(i*2))]<-t(tmp[,,i])
      Gnew[,(i*2)-1]<-names(A[[i]])[Gnew[,(i*2)-1]]
      Gnew[,(i*2)]<-names(A[[i]])[Gnew[,(i*2)]]
    }
  }else{
    for(i in 1:length(A)){
      Gnew[,i*2]<-tmp[,,i]
      Gnew[,(i*2)-1]<-0
      hom1<-which(Gnew[,i*2]==2)
      if(length(hom1)>0){
        Gnew[,(i*2)-1:0][which(Gnew[,i*2]==2),]<-1
      }
    }
  }
    Gnew<-Gnew[match(ped[,1], oped[,1]),]
    gens <- list()
    for (i in 1:(length(Gnew[1, ])/2)){
       gens[[i]] <- genotype(as.matrix(Gnew[, ((i * 2) - 1):(i *2)]), alleles=names(A[[i]]), reorder="no")
       names(gens)[i] <- namesG[i]
    }

  list(G=gens,valid=output[[9]]) 

}
       
