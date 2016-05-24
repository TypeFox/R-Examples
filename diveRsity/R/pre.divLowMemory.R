################################################################################
# pre.divLowMemory, a low memory consumption function for locus bootstrapping  #
################################################################################
pre.divLowMemory <- function(y){
  locs <- y$locs
  fst <- y$fst
  min <- y$min
  if(is.null(min)){
    min = TRUE
  }
  # define all functions first
  # define readGenepopX function
  #############################################################################
  # readGenepopX, a function for the generation of basic population parameters#
  #############################################################################
  readGenepopX <- function (x) {
    infile=x$infile
    #gp=x$gp
    bootstrap=x$bootstrap
    # define file reader
    ###########################################################################
    # Master file reader
    ###########################################################################
    fileReader <- function (infile) {
      if (typeof(infile) == "list") {
        return(infile)
      } else if (typeof(infile) == "character") {
        flForm <- strsplit(infile, split = "\\.")[[1]]
        ext <- flForm[[length(flForm)]]
        if (ext == "arp") {
          convRes <- arp2gen(infile)
          if (!is.null(convRes)) {
            cat("Arlequin file converted to genepop format! \n")
            infile <- paste(flForm[1], ".gen", sep = "")
          } else {
            infile <- paste(flForm[1], ".gen", sep = "")
          }
        }
        dat <- scan(infile, sep = "\n", what = "character", quiet = TRUE)
        if(length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 1){
          locs <- strsplit(dat[2], split = "\\s+")[[1]]
          if(length(locs != 1)){
            locs <- strsplit(dat[2], split = ",")[[1]]
          }
          locs <- as.character(sapply(locs, function(x){
            x <- strsplit(x, split = "")[[1]]
            if(is.element(",", x)){
              x <- x[-(which(x == ","))]
            }
            return(paste(x, collapse = ""))
          }))
          dat <- c(dat[1], locs, dat[-(1:2)])
        }
        
        
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        if (popLoc[1] == 3) {
          locs <- unlist(strsplit(dat[2], split = c("\\,", "\\s+")))
          dat <- c(dat[1], locs, dat[3:length(dat)])
        }
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        dat1 <- sapply(dat, function(x) {
          x <- unlist(strsplit(x, split = "\\s+"))
          if (is.element("", x)) {
            x <- x[-(which(x == ""))]
          }
          if (is.element(",", x)) {
            x <- x[-(which(x == ","))]
          }
          if (length(x) != 1 && length(x) != no_col) {
            x <- paste(x, collapse = "")
          }
          if (length(x) < no_col) {
            tabs <- paste(rep(NA, (no_col - length(x))), 
                          sep = "\t", collapse = "\t")
            line <- paste(x, tabs, sep = "\t")
            line <- unlist(strsplit(line, split = "\t"))
            return(line)
          } else {
            return(x)
          }
        })
      }
      out <- as.data.frame(t(dat1))
      rownames(out) <- NULL
      return(out)
    }
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
    npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                    which(data1[,1]=="pop")))
    pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
                which(data1[,1]=="pop"),(nrow(data1)+1))
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
    
    ###vectorize loci_pop_sizes###############################################
    
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
    
    
    #vectorized loci_pop_weights##############################################
    
    pre_loc_weights<- 1/pls
    loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
    loci_harm_N<-npops/colSums(pre_loc_weights)
    
    #end vectorized loci_pop_weights##########################################
    
    ###vectorize pop_alleles##################################################
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
    #end vectorize pop_alleles################################################
    
    #vectorize allele_names###################################################
    
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
    
    #all_alleles vectorized###################################################
    
    aaList<-function(x){
      return(sort(unique(x,decreasing=FALSE)))
    }
    all_alleles<-lapply(loci_combi,aaList)
    
    #end all_alleles vectorized###############################################
    
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
    ##########################################################################
    list(pop_list = pop_list,
         npops = npops,
         nloci = nloci,
         pop_sizes = pop_sizes,
         pop_alleles = pop_alleles,
         all_alleles = all_alleles,
         allele_freq = allele_freq,
         loci_harm_N = loci_harm_N,
         loci_names = loci_names,
         pop_names = pop_names,
         indtyp = indtyp,
         gp = gp)
  }
  ############################################################################
  # readGenepopX end                                                          #
  ############################################################################
  #
  #
  ############################################################################
  data1 <- readGenepopX(y)
  ############################################################################
  if(fst){
    # define the fst function
    ##########################################################################
    # fstWC: a function co calculate weir and cockerhams fis, fit, and fst
    ##########################################################################
    fstWC<-function(x){
      badData <- sapply(x$indtyp, function(y){
        is.element(0, y)
      })
      if(sum(badData) > 0){
        nl <- x$nloci - (sum(badData))
      } else{
        nl <- x$nloci
      }
      gdData<-which(!badData)
      badData<-which(badData)
      if (nl == 1) {
        all_genot<-x$pop_list[[1]][,gdData]
        if(x$npops > 1){
          for(i in 2:x$npops){
            all_genot <- c(all_genot, x$pop_list[[i]][,gdData])
          }
        }
        all_genot <- matrix(all_genot, ncol = 1)
      } else {
        all_genot<-matrix(x$pop_list[[1]][,gdData], ncol = length(gdData))
        if(x$npops > 1){
          for(i in 2:x$npops){
            all_genot<-rbind(all_genot, x$pop_list[[i]][,gdData])
          }
        }
      }
      genot<-apply(all_genot,2,unique)
      genot<-lapply(genot, function(x){
        if (sum(is.na(x))>0){
          y<-which(is.na(x)==TRUE)
          x_new<-x[-y]
          return(x_new)
        } else {
          return(x)
        }
      })
      #count genotypes
      
      genoCount<-list()
      for(i in 1:ncol(all_genot)){
        genoCount[[i]]<-matrix(0,ncol=length(genot[[i]]))
        for(j in 1:length(genot[[i]])){
          genoCount[[i]][,j]<-length(which(all_genot[,i] == genot[[i]][j]))
        }
        if (x$gp==3){
          colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,3),"/",
                                          substr(genot[[i]],4,6),sep="")
        } else if (x$gp==2){
          colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,2),"/",
                                          substr(genot[[i]],3,4),sep="")
        }
      }
      
      h_sum<-list()
      for(i in 1:ncol(all_genot)){
        h_sum[[i]]<-vector()
        cnSplit<-strsplit(colnames(genoCount[[i]]),"/")
        for(j in 1:length(x$all_alleles[[gdData[i]]])){
          het_id1<-lapply(cnSplit, is.element, x$all_alleles[[gdData[i]]][j])
          het_id2<-lapply(het_id1, sum)
          het_id2<-as.vector(het_id2)
          het_id3<-which(het_id2==1)
          h_sum[[i]][j]<-sum(genoCount[[i]][1,het_id3])
        }
      }
      indtyp_tot<-lapply(x$indtyp, sum)
      kk_hsum <- lapply(1:ncol(all_genot), function(i){
        list(h_sum[[i]], indtyp_tot[[gdData[i]]])
      })
      kk_hbar<-lapply(kk_hsum, function(x){
        return(x[[1]]/x[[2]])
      })
      
      pdat <- lapply(1:ncol(all_genot), function(i){
        list(x$allele_freq[[gdData[i]]], x$indtyp[[gdData[i]]])
      })
      
      kk_p<-lapply(pdat, function(x){
        if(is.null(x[[1]])==FALSE){
          apply(x[[1]], 1, function(y){
            y*(2*x[[2]])
          })
        }
      })
      res<-matrix(0,(x$nloci+1),3)
      colnames(res)<-c("Fis_WC","Fst_WC","Fit_WC")
      rownames(res)<-c(x$loci_names, "All")
      A<-vector()
      a<-vector()
      b<-vector()
      c<-vector()
      for(i in 1:ncol(all_genot)){
        kknbar<-indtyp_tot[[gdData[i]]]/x$npops
        kknC<-(indtyp_tot[[gdData[i]]]-sum(x$indtyp[[gdData[i]]]^2)/
                 indtyp_tot[[gdData[i]]])/(x$npops-1)
        kkptild<-kk_p[[i]]/(2*x$indtyp[[gdData[i]]])
        kkptild[kkptild=="NaN"]<-NA
        kkpbar<-colSums(kk_p[[i]])/(2*indtyp_tot[[gdData[i]]])
        kks2<-colSums(x$indtyp[[gdData[i]]] * (kkptild-rep(kkpbar, each = x$npops))^2)/((x$npops-1)*kknbar)
        kkA <- kkpbar * (1-kkpbar)-(x$npops-1)*kks2/x$npops
        kka<-kknbar*(kks2-(kkA-(kk_hbar[[i]]/4))/(kknbar-1))/kknC
        kkb <- (kknbar/(kknbar - 1))*(kkA-((2*kknbar-1)/(4*kknbar))*kk_hbar[[i]])
        #kkb<-kknbar*(kkA-(2*(kknbar-1))*kk_hbar[[i]]/(4*kknbar))/(kknbar-1)
        kkc<-kk_hbar[[i]]/2
        A[i]<-sum(kkA)
        a[i]<-sum(kka)
        b[i]<-sum(kkb)
        c[i]<-sum(kkc)
        res[gdData[i],"Fis_WC"]<- round(1-sum(kkc)/sum(kkb+kkc),4)
        res[gdData[i],"Fst_WC"]<- round(sum(kka)/sum(kka+kkb+kkc),4)
        res[gdData[i],"Fit_WC"]<- round(1-sum(kkc)/sum(kka+kkb+kkc),4)
      }
      res[res=="NaN"]<-NA
      res[res==0.000]<-NA
      sumA<-sum(na.omit(A))
      suma<-sum(na.omit(a))
      sumb<-sum(na.omit(b))
      sumc<-sum(na.omit(c))
      res[(x$nloci+1),"Fis_WC"]<-round(1-sumc/(sumb+sumc),4)
      res[(x$nloci+1),"Fst_WC"]<-round(suma/(suma+sumb+sumc),4)
      res[(x$nloci+1),"Fit_WC"]<-round(1-sumc/(suma+sumb+sumc),4)
      #res[is.na(res)]<-NaN
      list(Fstats=res,
           multiLoc<-res[(x$nloci+1),])
    }
    ##########################################################################
    # end fstWC
    ##########################################################################
    
    if(locs == TRUE){
      fstats <- fstWC(data1)[[1]]
    }else if (locs == FALSE){
      fstats <- fstWC(data1)[[2]]
    }
  }
  ##############################################################################
  # create 'easy use' objects from data1 (readGenepopX output)
  # pl = pop_list
  pl<-data1$pop_list
  # np = npops
  np<-data1$npops
  # nl = nloci
  nl<-data1$nloci
  # ps = pop sizes
  ps<-data1$pop_sizes
  # pa = pop alleles
  pa<-data1$pop_alleles
  # ant = allele names total
  ant<-data1$all_alleles
  # af = allele frequencies
  af<-data1$allele_freq
  # lnharm = locus harmonic sample size
  lnharm<-round(as.numeric(data1$loci_harm_N), 4)
  # ln = locus names
  ln<-data1$loci_names
  # pn = population names
  pn<-data1$pop_names
  # ntpl = number (of individuals) typed per locus
  nt<-data1$indtyp
  
  # remove data1 to save ram
  rm(data1)
  # garbage collect
  zz <- gc(reset = TRUE)
  rm(zz)
  ##############################################################################
  #observed heterozygosity count vectorize######################################
  
  ohcFUN<-function(x){
    lapply(1:ncol(x[[1]]), function(y){
      (x[[1]][,y]!=x[[2]][,y])*1 #multiply by 1 to conver logical to numeric
    })
  }
  ohc_data<-lapply(pa, ohcFUN)
  ohcConvert<-function(x){
    matrix(unlist(x),nrow=length(x[[1]]))
  }
  ohc<-lapply(ohc_data,ohcConvert)
  
  
  #end observed heterozygosity count vectorize##################################
  #exhmf & exhtf vectorize######################################################
  
  # calculate Heterozygosity exp
  tsapply <- function(...){t(sapply(...))}
  exhtf <- tsapply(af, function(x){
    apply(x, 2, function(y){
      1 - (sum(y^2))
    })
  })
  
  #end exhmf & exhtf vectorize##################################################
  #mean frequency vectorize#####################################################
  
  mf<-lapply(af,function(x){
    rowSums(x)/np  
  })
  ht<-sapply(mf, function(x){
    1- sum(x^2)
  })
  ht[ht=="NaN"]<-NA
  
  #end mean frequency vectorize#################################################
  
  
  ###end locus stats legacy code
  #locus stats vectorize########################################################
  
  hs<-round(rowSums(exhtf)/np,4)
  hs_est<-round(hs*((2*lnharm)/((2*lnharm)-1)),4)
  ht_est<-round((ht + (hs_est/(2*lnharm*np))),4)
  ht_est[ht_est=="NaN"]<-NA
  hst<-(ht-hs)/(1-hs)
  dst<-ht-hs
  gst<-dst/ht
  djost<-((ht-hs)/(1-hs))*(np/(np-1))
  hst_est<-(ht_est-hs_est)/(1-hs_est)
  dst_est<-ht_est- hs_est
  gst_est<-(ht_est-hs_est)/ht_est
  gst_max<-((np-1)*(1-hs))/(np-1+hs)
  gst_est_max<-(((np-1)*(1-hs_est))/(np-1+hs_est))
  gst_hedrick<-gst/gst_max
  gst_est_hedrick<-gst_est/gst_est_max
  gst_est_hedrick[gst_est_hedrick > 1] <- 1
  djost_est<-(np/(np-1))*((ht_est-hs_est)/(1 - hs_est))
  
  #end locus stats vectorize####################################################
  # Across all loci stats #
  ht_mean<-round(mean(na.omit(ht)),4)
  hs_mean<-round(mean(hs),4)
  gst_all<-round((ht_mean-hs_mean)/ht_mean,4)
  gst_all_max<-round(((np-1)*(1-hs_mean))/(np-1+hs_mean),4)
  gst_all_hedrick<-round(gst_all/gst_all_max,4)
  djost_all<-round(((ht_mean-hs_mean)/(1-hs_mean))*(np/(np-1)),4)
  ##############################################################################
  # Across all loci estimated stats #
  hs_est_mean<-round(mean(hs_est),4)
  ht_est_mean<-round(mean(na.omit(ht_est)),4)
  gst_est_all<-round((ht_est_mean-hs_est_mean)/ht_est_mean,4)
  gst_est_all_max<-round((((np-1)*(1-hs_est_mean))/(np-1+hs_est_mean)),4)
  gst_est_all_hedrick<-round(gst_est_all/gst_est_all_max,4)
  gst_est_all_hedrick[gst_est_all_hedrick > 1] <- 1
  #djost_est_all<-round((np/(np-1))*((ht_est_mean-hs_est_mean)/
  #(1 - hs_est_mean)),4)
  if (nl == 1){
    djost_est_all <- round(djost_est,4)
  } else {
    djost_est_all<-round(1/((1/mean(na.omit(djost_est))+(var(na.omit(djost_est))*
                                                           ((1/mean(na.omit(djost_est)))^3)))),4)
  }
  djost_est[djost_est==0]<-NaN
  djost[djost==0]<-NaN
  ##############################################################################
  if(fst == TRUE){
    if(locs == TRUE && min == FALSE){  
      list(hs=hs,
           hst=hst,
           dst=dst,
           gst=gst,
           djost=djost,
           hs_est=hs_est,
           ht_est=ht_est,
           hst_est=hst_est,
           dst_est=dst_est,
           gst_est=gst_est,
           djost_est=djost_est,
           gst_max=gst_max,
           gst_est_max=gst_est_max,
           gst_hedrick=gst_hedrick,
           gst_est_hedrick=gst_est_hedrick,
           ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn,
           fstats=fstats)
    } else if (locs == TRUE && min == TRUE){
      list(hs=hs,
           hst=hst,
           dst=dst,
           gst=gst,
           djost=djost,
           hs_est=hs_est,
           ht_est=ht_est,
           hst_est=hst_est,
           dst_est=dst_est,
           gst_est=gst_est,
           djost_est=djost_est,
           gst_max=gst_max,
           gst_est_max=gst_est_max,
           gst_hedrick=gst_hedrick,
           gst_est_hedrick=gst_est_hedrick,
           ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           #pop_list=pl,
           pop_names=pn,
           fstats=fstats)
    } else if (locs == FALSE && min == FALSE){
      list(ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn,
           fstats=fstats)
    } else if(locs == FALSE && min == TRUE){
      list(ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           #pop_list=pl,
           pop_names=pn,
           fstats=fstats)
    }
  } else {
    if(locs==T && min == FALSE){  
      list(hs=hs,
           hst=hst,
           dst=dst,
           gst=gst,
           djost=djost,
           hs_est=hs_est,
           ht_est=ht_est,
           hst_est=hst_est,
           dst_est=dst_est,
           gst_est=gst_est,
           djost_est=djost_est,
           gst_max=gst_max,
           gst_est_max=gst_est_max,
           gst_hedrick=gst_hedrick,
           gst_est_hedrick=gst_est_hedrick,
           ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn)
    } else if (locs == TRUE && min == TRUE){
      list(hs=hs,
           hst=hst,
           dst=dst,
           gst=gst,
           djost=djost,
           hs_est=hs_est,
           ht_est=ht_est,
           hst_est=hst_est,
           dst_est=dst_est,
           gst_est=gst_est,
           djost_est=djost_est,
           gst_max=gst_max,
           gst_est_max=gst_est_max,
           gst_hedrick=gst_hedrick,
           gst_est_hedrick=gst_est_hedrick,
           ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           #pop_list=pl,
           pop_names=pn)
    } else if (locs == FALSE && min == FALSE){
      list(ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn)
    } else if(locs == FALSE && min == TRUE){
      list(ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           #pop_list=pl,
           pop_names=pn)
    }
  }
}
################################################################################
# end pre.divLowMemory                                                         #
################################################################################