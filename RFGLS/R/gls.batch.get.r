#Adapted from Xiang Li's gls.batch(), by Rob Kirkpatrick.
gls.batch.get <-
  function(phenfile,genfile,pedifile,covmtxfile.in=NULL,theta=NULL, #Input arguments
           snp.names=NULL,input.mode=c(1,2,3),pediheader=FALSE, 
           pedicolname=c("FAMID","ID","PID","MID","SEX"),sep.phe=" ",sep.gen=" ",sep.ped=" ",
           phen,covars=NULL,med=c("UN","VC"), #Arguments for what to do with the input
           outfile,col.names=TRUE,return.value=FALSE, #Output arguments
           covmtxfile.out=NULL,
           covmtxparams.out=NULL,
           sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL #Optional arguments that trigger checks.
  ){
  ###########################################################
    
  #Idiot check for genfile: 
  if(missing(genfile)){ #genfile is also too important to just default to NULL.
    stop("Argument 'genfile' is missing, with no default (though value of NULL is accepted).")
  }
  
  #Miscellaneous input checks:
  if(length(input.mode)>1){input.mode <- input.mode[1]}
  if( !(input.mode %in% c(1,2,3)) ){
    warning("Argument 'input.mode' must be either 1, 2, or 3; coercing value to 1.")
    input.mode <- 1
  }
  #if(length(med)>1){med <- med[1]}
  #if( !(med %in% c("UN","VC")) ){
  #  stop("Argument 'med' must be a character string, and either 'UN' or 'VC'.")
  #}
  
  print("Reading in data.")
  
  #Phenotype file:
  if(is.character(phenfile)){phen.dat=read.table(phenfile,header=TRUE,sep=sep.phe,as.is=T)}
  else{phen.dat <- phenfile}
  if(all(c("FTYPE","INDIV") %in% colnames(phen.dat))){FSV.in.phenfile <- TRUE}
  else{
    FSV.in.phenfile <- FALSE
    if(input.mode==1){
      warning("Columns 'FTYPE' and 'INDIV' not found in phenotype file; input.mode coerced to 2.")
      input.mode <- 2
    }}     
  #Idiot check on phenotype file--is it sorted correctly?:
  if(input.mode==1){
    if( !(all(phen.dat[,c("FAMID","INDIV")]==phen.dat[order(phen.dat$FAMID, phen.dat$INDIV),c("FAMID","INDIV")])) ){
      warning("Phenotype file not sorted by FAMID, and by INDIV within FAMID; doing sorting now.")
      phen.dat <- phen.dat[order(phen.dat$FAMID, phen.dat$INDIV),]
    }}
  
  #Pedigree file:
  if(all(c("FTYPE","INDIV") %in% pedicolname)){FSV.in.pedifile <- TRUE}
  else{
    FSV.in.pedifile <- FALSE
    if(input.mode==2){
      warning("Columns 'FTYPE' and 'INDIV' not found in pedigree file; input.mode coerced to 3.")
      input.mode <- 3
    }}
  if(is.character(pedifile)){
    pedi.dat <- read.table(pedifile,header=pediheader,sep=sep.ped)#[,1:length(pedicolname)]
  }
  else{pedi.dat <- pedifile}
  if(length(pedicolname)>=ncol(pedi.dat)){pedi.dat <- pedi.dat[,1:length(pedicolname)]}
  else{stop("Length of argument 'pedicolname' is greater than the number of columns in pedigree file.")}
  names(pedi.dat) <- pedicolname
  if(input.mode==2){
    pedi.dat2 <- pedi.dat[,c("FTYPE","INDIV")]
    rownames(pedi.dat2) <- pedi.dat$ID
    phen.dat$FTYPE <- pedi.dat2[as.character(phen.dat$ID),"FTYPE"]
    phen.dat$INDIV <- pedi.dat2[as.character(phen.dat$ID),"INDIV"]
    rm(pedi.dat2)
    if( !(all(phen.dat[,c("FAMID","INDIV")]==phen.dat[order(phen.dat$FAMID, phen.dat$INDIV),c("FAMID","INDIV")])) ){
      phen.dat <- phen.dat[order(phen.dat$FAMID, phen.dat$INDIV),]
    }}
  
  #Genotype file:
  if(is.null(genfile)){test.dat <- NULL}
  else{
    if(is.character(genfile)){
      test.dat <- read.table(genfile,header=F,na.strings="NA",sep=sep.gen)
      test.dat <- data.frame(t(test.dat))
    } else {
      test.dat <- genfile
      rm(genfile)
    }
    #rownames(test.dat) <- pedi.dat$ID
  }
  gc()
  
  #SNP names:
  if(!is.null(test.dat)){
    if(is.character(snp.names)){
      if(length(snp.names)==length(names(test.dat))){names(test.dat) <- snp.names}
      else{
        warning("Mismatch between length of 'snp.names' and number of SNPs in 'genfile'; using generic snp.names instead.")
        names(test.dat) <- paste("snp",1:dim(test.dat)[2],sep=".")
      }
    }
    if(is.null(snp.names)){names(test.dat) <- paste("snp",1:dim(test.dat)[2],sep=".")}
    if(!is.null(snp.names) & !is.character(snp.names)){
      warning("Non-NULL value for argument 'snp.names' must be a character vector; using generic snp.names instead.")
    }
  }
  snplist <- names(test.dat) #<--Will be NULL if test.dat is also, which happens if genfile is NULL.
  
  if(input.mode==3){ #If family structure variables must be inferred from pedigree file.
    phen.dat <- FSV.frompedi(pedi.dat=pedi.dat,phen.dat=phen.dat)
    if( !(all(phen.dat[,c("FAMID","INDIV")]==phen.dat[order(phen.dat$FAMID, phen.dat$INDIV),c("FAMID","INDIV")])) ){
      phen.dat <- phen.dat[order(phen.dat$FAMID, phen.dat$INDIV),]
    }}
  
  print("Done reading in data.")
  
  #Data merging:
  if(!is.null(test.dat)){test.dat$ID <- pedi.dat$ID}
  idlab <- "ID" #"iid"
  result <- NULL
  famid <- "FAMID" #"fid"
  famtype <- "FTYPE" #"ftype"
  sid <- "INDIV"
  if(!is.null(covars)){phen.dat <- na.omit(phen.dat[,c(idlab,famid,famtype,sid,phen,covars)])}
  else{phen.dat <- na.omit(phen.dat[,c(idlab,famid,famtype,sid,phen)])}
  if(is.null(test.dat)){test.dat <- phen.dat}
  else{test.dat <- merge(phen.dat,test.dat,by="ID",sort=F)}
  print("Done merging data and trimming out incomplete cases.")
  
  #Idiot checks pertaining to family types:
  if( !is.null(Mz) ){
    if(Mz!=any(test.dat$FTYPE==1)){
      warning("Value supplied for 'Mz' not consistent with family types actually present in trimmed & merged data.")
    }}
  if( !is.null(Bo) ){
    if( Bo!=any(test.dat$FTYPE %in% c(2,4)) ){
      warning("Value supplied for 'Bo' not consistent with family types actually present in trimmed & merged data.")
    }}
  if( !is.null(Ad) ){
    if( Ad!=any(test.dat$FTYPE==3) ){
      warning("Value supplied for 'Ad' not consistent with family types actually present in trimmed & merged data.")
    }}
  if( !is.null(Mix) ){
    if( Mix!=any(test.dat$FTYPE==5) ){
      warning("Value supplied for 'Mix' not consistent with family types actually present in trimmed & merged data.")
    }}
  if( !is.null(indobs) ){
    if( indobs!=any(test.dat$FTYPE==6) ){
      warning("Value supplied for 'indobs' not consistent with family types actually present in trimmed & merged data.")
    }}
  #End of these idiot checks.
  
  #create famsize column
  test.dat$famsize = 1
  test.dat$famsize[test.dat$FTYPE!=6]=ave(test.dat$FAMID[test.dat$FTYPE!=6],test.dat$FAMID[test.dat$FTYPE!=6],FUN=length)
  #create unisid column, c-mz twin, b-bio-offspring, a-adopted offspring, f-father, m-mother
  test.dat$unisid=NULL
  test.dat$unisid[test.dat$INDIV==4 & test.dat$FTYPE!=6]="f"
  test.dat$unisid[test.dat$INDIV==3 & test.dat$FTYPE!=6]="m"
  test.dat$unisid[test.dat$FTYPE==1 & test.dat$INDIV==1]="c"
  test.dat$unisid[test.dat$FTYPE==1 & test.dat$INDIV==2]="c"
  test.dat$unisid[test.dat$FTYPE==2 & test.dat$INDIV==1]="b"
  test.dat$unisid[test.dat$FTYPE==2 & test.dat$INDIV==2]="b"
  test.dat$unisid[test.dat$FTYPE==4 & test.dat$INDIV==1]="b"
  test.dat$unisid[test.dat$FTYPE==4 & test.dat$INDIV==2]="b"
  test.dat$unisid[test.dat$FTYPE==3 & test.dat$INDIV==1]="a"
  test.dat$unisid[test.dat$FTYPE==3 & test.dat$INDIV==2]="a"
  test.dat$unisid[test.dat$FTYPE==5 & test.dat$INDIV==1]="b"
  test.dat$unisid[test.dat$FTYPE==5 & test.dat$INDIV==2]="a"
  #test.dat$unisid[test.dat$INDIV==5] <- "u"
  #create fam labs
  test.dat$famlab="INDPT"
  test.dat$famlab[test.dat$FTYPE!=6] = ave(test.dat$unisid[test.dat$FTYPE!=6],test.dat$FAMID[test.dat$FTYPE!=6],
                                           FUN=function(x) do.call("paste",c(data.frame(matrix(x,1,length(x))),sep="")))
  #get tlist and famsize list; tlist is the list of family labels, and famsize is the list of family sizes
  tlist = tapply(test.dat$famlab[test.dat$famlab!="INDPT"],test.dat$FAMID[test.dat$famlab!="INDPT"],FUN=function(x) unique(x))
  names=as.character(unique(test.dat$FAMID[test.dat$famlab!="INDPT"]))
  tlist=tlist[names]
  tlist = c(tlist,rep("INDPT",sum(test.dat$famlab=="INDPT")))
  sizelist = tapply(test.dat$famsize[test.dat$famlab!="INDPT"],test.dat$FAMID[test.dat$famlab!="INDPT"],FUN=function(x) unique(x))
  sizelist = sizelist[names]
  sizelist = c(sizelist ,rep(1,sum(test.dat$famlab=="INDPT")))
  test.dat <- rbind(test.dat[test.dat$famlab != "INDPT",], test.dat[test.dat$famlab=="INDPT",])
  id<-test.dat[,idlab]
  
  #Idiot checks pertaining to family sizes:
  if( !is.null(sizeLab)  ){
    if( !is.character(sizeLab) ){
      warning("Value provided for 'sizeLab' is not a character string, and will be ignored.")
    }
    else{
      warning("Use of argument 'sizeLab' is optional, and may be eliminated in future package versions.")
      if(nchar(sizeLab)>max(sizelist)){
        warning("Value provided for 'sizeLab' specifies a larger maximum family size than appears in the trimmed & merged data.")
      }
      if(nchar(sizeLab)<max(sizelist)){
        warning("Value provided for 'sizeLab' specifies a smaller maximum family size than appears in the trimmed & merged data.")
      }
    }}
  #End of these idiot checks.
  
  getout <- list(test.dat, tlist, sizelist)
  names(getout) <- c("test.dat", "tlist", "sizelist")
  return(getout)
}