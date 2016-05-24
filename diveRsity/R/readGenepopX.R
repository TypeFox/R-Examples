################################################################################
# readGenepopX, a function for the generation of basic population parameters   #
################################################################################
readGenepopX <- function (x) {
  #gp=x$gp
  infile=x$infile
  bootstrap=x$bootstrap
  locs=x$locs
  data1 <- fileReader(infile)
  if(is.null(x$gp)){
    rownames(data1) <- NULL
    data1 <- as.matrix(data1)
    # determine genepop format
    p1 <- which(toupper(data1[,1]) == "POP")[1] + 1
    gp <- as.numeric(names(sort(-table(sapply(data1[p1,(-1)], nchar)/2)))[1])
    data1 <- as.data.frame(data1)
  } else {
    gp <- x$gp
  }
  
  if(gp == 3){
    data1[data1==0]<-NA
    data1[data1=="999999"]<-NA
    data1[data1=="000000"]<-NA
    data1[data1=="NANA"]<-NA
  } else if(gp == 2){
    data1[data1==0]<-NA
    data1[data1=="9999"]<-NA
    data1[data1=="0000"]<-NA
    data1[data1=="NA"]<-NA
  }
  raw_data<-data1
  npops<-length(which(toupper(data1[,1]) == "POP"))
  pop_pos<- c(which(toupper(data1[,1]) == "POP"),(nrow(data1)+1))
  pop_sizes<-vector()
  for(i in 1:npops){
    pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
  }
  pop_names<-as.character(data1[(pop_pos[1:npops]+1),1])
  pop_weights<- 1/pop_sizes
  
  n_harmonic<-npops/sum(pop_weights)
  
  N<-pop_sizes
  
  nloci<- (pop_pos[1]-2)
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  pop_list<-list()
  for (i in 1:npops){
    pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                                   2:(nloci+1)])
  }
  # check if all populations have at least some data at loci
  extCheck <- sapply(1:length(pop_list), function(i){
    sum(is.na(pop_list[[i]])) == nloci * pop_sizes[i]
  })
  if (sum(extCheck) > 0){
    npops <- npops - sum(extCheck)
    pop_list <- pop_list[-(which(extCheck == TRUE))]
    pop_sizes <- pop_sizes[-(which(extCheck == TRUE))]
    pop_names <- pop_names[-(which(extCheck == TRUE))]
    pop_weights <- pop_weights[-(which(extCheck == TRUE))]
    N <- N[-(which(extCheck == TRUE))]
    #raw_data fix
    noPop <- which(extCheck == TRUE)
    indexer <- lapply(noPop, function(i){
      (pop_pos[i] + 1):(pop_pos[(i+1)])
    })
    indexer <- unlist(indexer)
    raw_data <- raw_data[-(indexer), ]    
  }  
  if (gp==3) {
    plMake<-function(x){
      out <- matrix(sprintf("%06g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "0000NA"] <- "    NA"
      }
      return(out)
    }
  } else if (gp==2) {
    plMake<-function(x){
      out <- matrix(sprintf("%04g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "00NA"] <- "  NA"
      }
      return(out)
    }
  }
  suppressWarnings(pop_list<-lapply(pop_list, plMake))
  
  
  if (gp == 3){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "    NA"]<-NA
    }
  } else if (gp == 2){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "  NA"] <-NA
    }
  }
  
  
  if(bootstrap == T){
    bs<-function(x){
      return(matrix(x[sample(nrow(x),replace=TRUE), ],ncol=ncol(x)))
    }
    pop_list<-lapply(pop_list, bs)
  }  
  
  ###vectorize loci_pop_sizes################################################
  
  lps<-function(x){#
    lsp_count<-as.vector(colSums(!is.na(x)))#
    return(lsp_count)#
  }#
  pre_loci_pop_sizes<-lapply(pop_list,lps)#
  pls<-matrix(ncol=nloci,nrow=npops)#
  for(i in 1:length(pre_loci_pop_sizes)){#
    pls[i,]<-pre_loci_pop_sizes[[i]]#
  }#
  #convert pls to loci_pop_sizes format
  loci_pop_sizes<-split(pls,col(pls))
  
  
  #vectorized loci_pop_weights##################################################
  
  pre_loc_weights<- 1/pls
  loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
  loci_harm_N<-npops/colSums(pre_loc_weights)
  
  #end vectorized loci_pop_weights##############################################
  
  ###vectorize pop_alleles########################################################
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
      pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
      return(pl)
    }
  }
  pop_alleles<-lapply(pop_list,pl_ss)
  #end vectorize pop_alleles####################################################
  
  #vectorize allele_names#######################################################
  
  alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
    res<-list()
    for(i in 1:ncol(x[[1]])){
      res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
    }
    return(res)
  }
  
  allele_names<-lapply(pop_alleles,alln)
  
  
  loci_combi<-allele_names[[1]]
  for(j in 1:nloci){
    for(i in 2:npops){
      loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
    }
  }
  
  #all_alleles vectorized#######################################################
  
  aaList<-function(x){
    return(sort(unique(x,decreasing=FALSE)))
  }
  all_alleles<-lapply(loci_combi,aaList)
  
  #end all_alleles vectorized###################################################
  
  aa<-all_alleles
  aa<-lapply(aa, FUN=`list`, npops)
  afMatrix<-function(x){
    np<-x[[2]]
    z<-matrix(rep(0,(np*length(x[[1]]))),ncol=np, nrow=length(x[[1]]))
    rownames(z)<-x[[1]]
    return(z)
  }
  allele_freq<-lapply(aa,afMatrix)
  
  
  #combine pop_alleles
  parbind<-function(x){
    rbind(x[[1]],x[[2]])
  }
  pa1<-lapply(pop_alleles, parbind)
  #create a function to tabulate the occurance of each allele
  afTab<-function(x){
    lapply(1:ncol(x), function(i){
      return(table(x[,i]))
    })
  }
  actab<-lapply(pa1, afTab)
  
  afs<-function(x){
    afsint<-function(y){
      length(na.omit(y))/2
    }
    apply(x,2,afsint)
  }
  indtyppop<-lapply(pa1,afs)
  #calculate allele frequencies
  afCalcpop<-lapply(1:length(actab), function(x){
    lapply(1:length(actab[[x]]),function(y){
      actab[[x]][[y]]/(indtyppop[[x]][y]*2)
    })
  })
  #assign allele freqs to frequency matrices
  obs_count<-allele_freq
  for(i in 1:npops){
    for(j in 1:nloci){
      allele_freq[[j]][names(afCalcpop[[i]][[j]]),i]<-afCalcpop[[i]][[j]]
      obs_count[[j]][names(actab[[i]][[j]]),i]<-actab[[i]][[j]]
    }
  }
  
  
  
  indtyp<-list()
  for(i in 1:nloci){
    indtyp[[i]]<-vector()
  }
  for(i in 1:npops){
    for(j in 1:nloci){
      indtyp[[j]][i]<-indtyppop[[i]][j]
    }
  }
  
  if(bootstrap==T){
    ind_vectors<-list()
    for(i in 1:npops){
      ind_vectors[[i]]<-noquote(paste(rep(i,pop_sizes[i]),",",sep=""))
    }
    
    
    pre_data<-matrix(rep("\t",((nloci+1)*(nloci+1))),
                     ncol=(nloci+1))
    pre_data[1,]<-c("Title",rep("\t",nloci))
    for(i in 2:(nloci+1)){
      pre_data[i,1]<-loci_names[(i-1)]
    }
    pop_data<-list()
    for(i in 1:npops){
      pop_data[[i]]<-matrix(rbind(c("POP",as.vector(rep("\t",nloci))),
                                  cbind(ind_vectors[[i]],pop_list[[i]])),
                            ncol=(nloci+1))
    }
    bs_data_file<-matrix(rbind(pre_data,pop_data[[1]]),ncol=(nloci+1))
    for(i in 2:npops){
      bs_data_file<-matrix(rbind(bs_data_file,pop_data[[i]]),ncol=(nloci+1))
    }
    bs_data_file<-data.frame(bs_data_file)
  }
  nalleles<-vector()
  for(i in 1:nloci){
    nalleles[i]<- nrow(allele_freq[[i]])
  }
  ##############################################################################
  if(bootstrap==T){
    list(npops=npops, 
         nloci=nloci, 
         pop_alleles=pop_alleles, 
         pop_list=pop_list,
         loci_names=loci_names, 
         pop_pos=pop_pos, 
         pop_sizes=pop_sizes,
         allele_names=allele_names,
         all_alleles=all_alleles,
         allele_freq=allele_freq,
         raw_data=raw_data,
         loci_harm_N=loci_harm_N,
         n_harmonic=n_harmonic,
         pop_names=pop_names,
         indtyp=indtyp,
         nalleles=nalleles,
         locs=locs,
         bs_file=bs_data_file,
         obs_allele_num=obs_count)
  } else if(bootstrap==F){
    list(npops=npops, 
         nloci=nloci, 
         pop_alleles=pop_alleles, 
         pop_list=pop_list,
         loci_names=loci_names, 
         pop_pos=pop_pos, 
         pop_sizes=pop_sizes,
         allele_names=allele_names,
         all_alleles=all_alleles,
         allele_freq=allele_freq,
         raw_data=raw_data,
         loci_harm_N=loci_harm_N,
         n_harmonic=n_harmonic,
         pop_names=pop_names,
         indtyp=indtyp,
         nalleles=nalleles,
         locs=locs,
         obs_allele_num=obs_count)
  }
}
################################################################################
# readGenepopX end                                                             #
################################################################################