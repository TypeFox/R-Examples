## This function was developed based on the lme.batch function in GWAF package
## by  Qiong Yang and Ming-Huei Chen.
## Function written by Xiang Li (last update 5/18/11) and Saonli Basu (last update 7/18/12), Rob Kirkpatrick (last update May 2013)
###############################################################################
gls.batch <- 
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
    
  #Idiot checks for outfile:
  if( missing(outfile) | (!is.null(outfile) & !is.character(outfile)) ){ #outfile is too important to default to NULL.
    stop("Argument 'outfile' must be a character string, and has no default (though value of NULL is accepted).")
  }
  if(is.null(outfile) & return.value==FALSE){
    stop("Arguments 'outfile=NULL' and 'return.value=FALSE' are mutually exclusive.")
  }
  if(is.character(outfile)){
    #Idiot check for nonexistent outfile directory (in  which case this writing will fail):
    write.table(x=" ", file=outfile, quote=F, row.names=F, col.names=F) 
  }

  #Idiot checks for covmtxfile.out:
  if( !is.null(covmtxfile.out) ){
    if( !is.character(covmtxfile.out) ){
      warning("Argument 'covmtxfile.out' not a character string; residual covariance matrix will not be written to disk.")
    }
    #Idiot check for nonexistent covmtxfile.out directory (in  which case this writing will fail):
    else{write.table(x=" ", file=covmtxfile.out, quote=F, row.names=F, col.names=F)}
  }

  #Idiot checks for covmtxparams.out:
  if( !is.null(covmtxparams.out) ){
    if( !is.character(covmtxparams.out) ){
      warning("Argument 'covmtxparams.out' not a character string; residual covariance parameters will not be written to disk.")
    }
    #Idiot check for nonexistent covmtxparams.out directory (in  which case this writing will fail):
    else{write.table(x=" ", file=covmtxparams.out, quote=F, row.names=F, col.names=F)}
  }
  
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
  if(length(med)>1){med <- med[1]}
  if( !(med %in% c("UN","VC")) ){
    stop("Argument 'med' must be a character string, and either 'UN' or 'VC'.")
  }
  
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
  
  #Idiot checks pertaining to previously estimated residual covariance matrix:
  if( !is.null(covmtxfile.in) & !is.null(theta) ){stop("At least one of arguments 'covmtxfile.in' and 'theta' must be NULL.")}
  if( !is.null(theta) & (med=="UN" & length(theta)!=12) ){
    theta <- NULL
    warning("Non-NULL value to argument 'theta' must be a numerical vector of length 12 (NA's are accepted) for med='UN'; coercing 'theta' to NULL.")
  }
  if( !is.null(theta) & (med=="VC" & length(theta)!=3) ){
    theta <- NULL
    warning("Non-NULL value to argument 'theta' must be a numerical vector of length 3 or 12 (NA's are accepted) for med='VC'; coercing 'theta' to NULL.")
  }
  
  #Load residual covariance matrix (if specified):
  vmat <- NULL
  varfile <- NULL
  if(!is.null(covmtxfile.in)){
    if(is.character(covmtxfile.in)){
      if(file.exists(covmtxfile.in)==T){
        varfile <- read.csv(covmtxfile.in,header=T,colClasses="numeric")
      }
      else{
        varfile <- NULL
        warning("Specified file for 'covmtxfile.in' not found; will (re-)estimate residual covariance matrix.")
      }
    }
    else{vmat <- covmtxfile.in}
  }
  
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
  if(!is.null(covars)){
 	  phen.dat <- na.omit(phen.dat[,c(idlab,famid,famtype,sid,phen,covars)])
#     if( any( sapply(phen.dat[,covars],is.factor) + sapply(phen.dat[,covars],is.character) ) ){
#       stop("Variables named in argument 'covars' cannot be of mode factor or character.")
#     }
  }
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
  
  #Build matrix if parameters were supplied as 'theta':
  if( !is.null(theta) ){
    if(med=="VC"){blocks <- getblocks.ACE(theta=theta,tlist=tlist,sizelist=sizelist)}
    else{blocks <- getblocks(theta=theta,tlist=tlist,sizelist=sizelist,pad=FALSE)}
    vmat <- bdsmatrix(sizelist,blocks=blocks$blocks,dimnames=list(id,id))
  }

  #Estimate residual covariance matrix from covariates-only (or intercept-only) model:
  if(is.null(vmat)){
    if(is.null(varfile)){
	    if(is.null(covars)){
        lme.out <- try(fgls(test.dat[,phen]~1,data=test.dat,tlist=tlist, sizelist=sizelist, med=med))
        vmat <- bdsmatrix(sizelist,lme.out$sigma@blocks,dimnames=list(id,id)) 
        if(is.character(covmtxfile.out)){write.csv(vmat@blocks,file=covmtxfile.out,quote=F,row.names=FALSE)}
        if(is.character(covmtxparams.out)){write.table(lme.out$estimates,file=covmtxparams.out,row.names=F,col.names=F)}
      }
      else{
        if(length(covars)==1){x.covar <- model.matrix(as.formula(paste("~",covars)),data=test.dat)[,-1]}
        else{x.covar <- model.matrix(as.formula(
          paste("~",covars[1],paste("+",covars[2:length(covars)],sep="",collapse=""))),data=test.dat)[,-1]}  
        #x.covar <- as.matrix(test.dat[,covars])
    	  lme.out <- try(fgls(test.dat[,phen]~1+x.covar,data=test.dat,tlist=tlist, sizelist=sizelist, med=med))
        vmat <- bdsmatrix(sizelist,lme.out$sigma@blocks,dimnames=list(id,id))
        if(is.character(covmtxfile.out)){write.csv(vmat@blocks,file=covmtxfile.out,quote=F,row.names=FALSE)}
        if(is.character(covmtxparams.out)){write.table(lme.out$estimates,file=covmtxparams.out,row.names=F,col.names=F)}
	  }}
    else{vmat <- bdsmatrix(sizelist,c(varfile[,1]),dimnames=list(id,id))}
	} #By this point, vmat should be defined as other than NULL; possibly, varfile might not be.
  gc()

  #If there are no SNPs (because genfile=NULL), the program can terminate normally:
  if(is.null(snplist)){
    if(return.value==TRUE){
      if(exists("lme.out")){
        return(lme.out) #If a function value is requested, it will be the fgls() regresssion object.
      }
      if(is.null(covars)){return(try(fgls(
        test.dat[,phen]~1,data=test.dat,tlist=tlist, sizelist=sizelist,vmat=vmat,med=med)))
      }
      else{
        if(length(covars)==1){x.covar <- model.matrix(as.formula(paste("~",covars)),data=test.dat)[,-1]}
        else{x.covar <- model.matrix(as.formula(
          paste("~",covars[1],paste("+",covars[2:length(covars)],sep="",collapse=""))),data=test.dat)[,-1]} 
        return(try(fgls(test.dat[,phen]~1+x.covar,data=test.dat,tlist=tlist, sizelist=sizelist,vmat=vmat,med=med)))
      }
    }
    else{return()}
  }
  
  #RMK Jun'13: We can make the Cholesky factor of the inverted residual covariance matrix now, and re-use it for SNPs that
  #have complete data; this way, we're not re-inverting and re-factoring it on every SNP:
  list.vmat<-listbdsmatrix(vmat,diag=T,id=F)
  vmat1<-sparseMatrix(list.vmat[,2],list.vmat[,1],x=list.vmat[,3],symmetric=T)
  vmat.Inv<-as(solve(vmat1,full=T),"sparseMatrix")
  vmat.Inv<-forceSymmetric(vmat.Inv)
  gkmat<-as(chol(vmat.Inv),"sparseMatrix")
  rm(list.vmat,vmat1,vmat.Inv); gc()
  ################# begins single snp analysis ###################

#   #RMK May'13: Writing a function that does a 1-SNP-at-a-time analysis and applying it to the columns of the matrix
#   #of genotypes is faster than a 'for' loop.
#   if(is.null(covars)){
#     result <- apply(X=test.dat[,snplist],MARGIN=2,FUN=gwasfunc,pheno=test.dat[,phen],
#                     id=id, x.covar.full=NULL, vmat=vmat, tlist=tlist, sizelist=sizelist, med=med)
#   }
#   else{
#     result <- apply(X=test.dat[,snplist],MARGIN=2,FUN=gwasfunc,pheno=test.dat[,phen],
#                     id=id, x.covar.full=test.dat[,covars], vmat=vmat, tlist=tlist, sizelist=sizelist, med=med)
#   }
#   result <- data.frame(snp=snplist,coef=as.numeric(result[2,]),se=as.numeric(result[3,]),
#                        t.stat=as.numeric(result[4,]), df=as.integer(result[5,]), pval=as.numeric(result[6,])) 
  
  #RMK Jun'13: by creating objects first and then filling in their elements as the loop progresses, those objects' 
  #memory-addressing only has to be done once.
  m <- length(snplist)
  result <- data.frame(snp=snplist,coef=as.numeric(rep(NA,m)),se=as.numeric(rep(NA,m)),
                       t.stat=as.numeric(rep(NA,m)), df=as.integer(rep(NA,m)),
                       pval=as.numeric(rep(NA,m)), stringsAsFactors=F)
  if(is.null(covars)){X <- as.matrix(cbind(rep(1,nrow(test.dat)),rep(NA,nrow(test.dat))))}
  else{
    if(length(covars)==1){mm <- model.matrix(as.formula(paste("~",covars)),data=test.dat)[,-1]}
    else{mm <- model.matrix(as.formula(
      paste("~",covars[1],paste("+",covars[2:length(covars)],sep="",collapse=""))),data=test.dat)[,-1]}  
    #X <- model.matrix(terms(model.frame(lm(test.dat[,phen] ~ test.dat[,covars]))))
    X <- as.matrix(cbind(1,NA,mm))
  }
  Y <- test.dat[,phen]
  for(i in 1:m){ #start snplist loop
    if(i%%100==0){gc()}
    X[,2] <- test.dat[,snplist[i]]
    if(length(table(X[,2]))==1){next} #If SNP is monomorphic.
    if(any(is.na(X[,2]))){ #Check if anyone is missing the current SNP; if so, cut them out of the data for this iteration.
      na.rows <- which(is.na(X[,2]))
      X0 <- X[-na.rows,]
      Y0 <- test.dat[-na.rows,phen]
      vmat0 <- vmat[-na.rows,-na.rows]
      lme.out <- try(bare_fgls(Y=Y0, X=X0, vmat=vmat0, gkmat=NULL))
    }
    else{lme.out <- try(bare_fgls(Y=Y, X=X, vmat=vmat, gkmat=gkmat))}    
    if(class(lme.out)=="try-error"){next} #If bare_fgls() fails for some reason.
    else{
      result[i,2:6] <- c(lme.out$ctable[2,1],lme.out$ctable[2,2],lme.out$ctable[2,3],
                         lme.out$df.residual,lme.out$ctable[2,4])
}} #end of snplist loop
  
  rm(test.dat)
   result <- data.frame(snp=as.character(result[,1]),coef=as.numeric(result[,2]),se=as.numeric(result[,3]),
                          t.stat=as.numeric(result[,4]), df=as.integer(result[,5]), pval=as.numeric(result[,6]), 
                        stringsAsFactors=F) 
  if(!is.null(outfile)){
    write.table(result, outfile, quote=F,row.names=F, col.names=col.names,sep=" ",append=F)
  }
  
  if(return.value==T){return(result)}
  else{ return() }
}
