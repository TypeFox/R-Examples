single.snp.test.families <-
function ( snps, trait, adj.var=NULL, famid , patid ,
    fid , mid , prt=T  ) {


  snps <- as.matrix(snps)

  if ( !all(trait[!is.na(trait)] == 0 | trait[!is.na(trait)] == 1) )
    stop("trait should be 0 for unaffected or 1 for affected children.")
    
  # remove unaffected children
  not.rm.child <- ifelse ( (fid!=0) & (mid!=0) & (trait==0) , F , T )

  if ( sum(!not.rm.child)>0 ) {
  
    cat ( sum(!not.rm.child) , " children removed because they were unaffected.\n",sep="")

  }

  famid <- famid[not.rm.child]
  patid <- patid[not.rm.child]
  fid   <- fid[not.rm.child]
  mid   <- mid[not.rm.child]
  snps  <- snps[not.rm.child,,drop=FALSE]
  trait <- trait[not.rm.child]
  
  if ( length(famid)<1 ) {
    stop ("no families for TDT observed.")
  }


    # exclusion of nuclear families without two parents.
  excl.fam <- NULL
  #i <- unique(famid)[1]
  for ( i in unique(famid) ) {
    selfam <- famid==i
    selfid <- unique(fid[selfam])
    selfid <- selfid[ (!is.na(selfid)) & (selfid!=0) ]
    selmid <- unique(mid[selfam])
    selmid <- selmid[ (!is.na(selmid)) & (selmid!=0) ]
    if ( length(selfid)!=1  ) {
      excl.fam <- c(excl.fam,i)
    } else {
        if ( ! any(patid[selfam] %in% selfid) ) {
          excl.fam <- c(excl.fam,i)
        }
    }

    if ( length(selmid)!=1  ) {
      excl.fam <- c(excl.fam,i)
    } else {
        if ( ! any(patid[selfam] %in% selmid) ) {
          excl.fam <- c(excl.fam,i)
        }
    }
  }
  excl.fam <- unique(excl.fam)
  if ( !is.null(excl.fam) ) {
    print ( paste ( "Exclusion of nuclear families without two parents: " , 
                    paste(excl.fam,collapse=" ") , sep="") )
    selcond <-  !(famid %in% excl.fam)
    famid <- famid [ selcond   ]
    patid <- patid [ selcond  ]
    fid   <- fid   [ selcond  ]
    mid   <- mid   [ selcond  ]
    snps  <- snps  [ selcond  ,,drop=FALSE]
    trait <- trait [ selcond  ]
  }
  if ( length(famid)<1 ) {
    stop ("no families for TDT observed.")
  }
  
  ordfam <- order.families ( famid, patid, fid, mid )
  snps  <- snps[ordfam,,drop=FALSE]
  trait <- trait[ordfam]
  famid <- famid[ordfam]
  patid <- patid[ordfam]
  fid   <- fid[ordfam]
  mid   <- mid[ordfam]

  if ( !all(is.null(adj.var) ) ) {
     stop("Error in single.snp.test.families: no adjustment variables are allowed by using wTDT.")
  }
  
  freq.allele1 <- trans.allele.1 <- nontrans.allele1 <- pval <- rep(NA,dim(snps)[2])
  
  for ( i in 1:dim(snps)[2] ) {

    sht <- single.haplotype.test ( snps[,i], trait, famid , patid , fid , mid ,
    adj.var=adj.var , type = "families" , lim=1e-06 , prt=F )

    if ( dim(sht$haplotypes)[1] == 3 ) {
      freq.allele1 [i]       <- as.numeric(sht$haplotypes[ sht$haplotypes[,"Hap"]=="1" , "Freq" ])
      trans.allele.1 [i]    <- as.numeric(sht$haplotype.i [ "1" , "trans" ])
      nontrans.allele1 [i] <- as.numeric(sht$haplotype.i [ "1" , "non-trans" ])
      pval [i]               <- as.numeric(sht$global.test["pvalue"])
    } else {
      freq.allele1 [i]       <- NA
      trans.allele.1 [i]    <- NA
      nontrans.allele1 [i] <- NA
      pval [i]               <- NA

    }
  }
  
  res <- data.frame ( SNP=1:dim(snps)[2] , type="families", freq.allele1=freq.allele1 ,
      trans.allele.1=trans.allele.1, nontrans.allele.1=nontrans.allele1,
      OR=trans.allele.1/nontrans.allele1,test="wTDT",
      p.value=pval, stringsAsFactors=F)

  return ( res )
}
