BLUP=function(trait="yield",family="all",env="all",
              MAF=0.05,use.check=TRUE,impute="FM",rm.rep=TRUE){
    
    # Line added for debugging purpose
    gen.raw=data.check=data.line=matrix(NA,2,2)
  
    # Load data
    data(soynam,envir=environment(),package="SoyNAM")
        
    # Genotypic matrix of lines 
    geno = gen.raw[grep('DS1',rownames(gen.raw)),]
        
    # FAM
    fam=rownames(geno)
    fam=gsub('DS1.-','',fam)
    fam=gsub('...$','',fam,perl = T)
    fam=as.numeric(fam)
    
    # CHR
    chr=rep(NA,20)
    for(i in 1:20) chr[i]=length(grep(
      paste("Gm",sprintf("%02d",i),sep=''),
      colnames(geno)));rm(i)
    
    # Subsetting
    if(is.numeric(family)) data.line = data.line[data.line$family%in%family,]
    
    if(is.numeric(env))  if(length(env)==1) stop("At least two environments where the trait was measured are required")
    
    if(is.numeric(env)){
      E1 = as.numeric(data.line$environ)
      data.line = data.line[E1%in%env,]
      E2 = as.numeric(data.check$environ)
      data.check = data.check[E2%in%env,]
    }
    
    # Check function
    CHECK=function(trait){ test=dcast(data.check,environ+spot~strain,value.var=trait,mean)
                           rownames(test)=test[,2];E=test[,1];test=test[,-c(1,2)];test=data.matrix(test);test[is.nan(test)]=NA;
                           X=function(X) unlist(tapply(X,E,FUN=function(x){m=mean(x,na.rm=T);SD=sd(x,na.rm=T);return((x-m)/SD)}))
                           MEAN=apply(test,2,X);C=rowMeans(MEAN,na.rm=T);names(C)=rownames(test);C[is.nan(C)]=0;return(C)}
    
    # Model terms
    Y = data.line[,trait]
    G = data.line[,"strain"]
    E = data.line[,"environ"]
    
    # BLUP
    if(use.check){
      cat('solving BLUE of checks\n')
      check = CHECK(trait);set = as.character(data.line[,"spot"])
      C = check[set]
      cat('solving BLUP of phenotypes\n')
      blup=lmer(Y~C+(1|E)+(1|G))
    }else{
      cat('solving BLUP of phenotypes\n')
      blup=lmer(Y~(1|E)+(1|G))}
    BV = coef(blup)$G[,1]
    names(BV) = rownames(coef(blup)$G)
    BV = BV[rownames(geno)]
    
    if(is.numeric(family)){
      BV = BV[fam%in%family]
      geno = geno[fam%in%family,]
      fam = fam[fam%in%family]
    }
    
    geno = snpQC(geno,MAF=MAF)
    
    cat('Removing markers with more than 80% missing values\n')
    mi=apply(geno,2,function(q)mean(is.na(q)))
    geno = geno[,mi<.8]
    
    for(i in 1:20) chr[i]=length(grep(
      paste("Gm",sprintf("%02d",i),sep=''),
      colnames(geno)));rm(i)
    
    if(anyNA(geno)){
      
      if(impute=="RF"){
        CHR = as.numeric(substring(colnames(geno),1,1))
        for(i in 1:length(chr)){
          w = which(CHR==i)
          geno[,w] = snpQC(geno[,w],1,0,F,T)
          cat('DONE IMPUTING CHROMOSOME',i,'\n')
        }  
      }
      
      if(impute=="FM"){
        cat('Imputing with expectation (based on transition prob) \n')
        geno = markov(geno,chr)
      }
      
    }
        
    # CHR
    for(i in 1:20) chr[i]=length(grep(
      paste("Gm",sprintf("%02d",i),sep=''),
      colnames(geno)));rm(i)
    
    if(!anyNA(geno)){
      if(rm.rep){
        cat('removing repeated genotypes\n')
        d=cleanREP(y = cbind(BV,BV),fam = fam,gen = geno)
        BV=d$y[,1]
        geno=d$gen
        Chrom=chr
        fam=d$fam
      }
      rownames(geno) = names(BV)
    }
  
    LIST = list('Phen'=BV,'Gen'=geno,'Chrom'=chr,'Fam'=fam)
    
    return(LIST)
  }

