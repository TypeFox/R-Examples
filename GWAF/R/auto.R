
auto <- function(genopath,phenfile,pedfile,outfile,phen,covars=NULL,cov.int,sub="N",analysis,lib.loc,model=NULL,kinmat=NULL,col.names=F,sep.ped=",",sep.phe=",",sep.gen=","){

gfiles <- list.files(path=genopath,full.names=T)
gfs <- list.files(path=genopath,full.names=F)

if (length(gfs)==0) stop(paste("No genotype files in ",genopath,"!!",sep=""))
if (!analysis %in% c("lmepack","lmepack.imputed","lmeVpack.imputed","glmm","geepack","geepack.imputed","geepack.quant","geepack.quant.imputed","lmepack.int","lmepack.int.imputed","geepack.int","geepack.int.imputed","geepack.quant.int","geepack.quant.int.imputed")) 
   stop("Please choose one appropriate option from lmepack, lmepack.imputed, geepack, geepack.imputed, geepack.quant, geepack.quant.imputed, 
       lmepack.int, lmepack.int.imputed, geepack.int, geepack.int.imputed, geepack.quant.int, geepack.quant.int.imputed for analysis")
if (!file.exists(phenfile)) stop(paste(phenfile," does not exist!",sep=""))
if (!file.exists(pedfile)) stop(paste(pedfile," does not exist!",sep=""))
if (file.exists(outfile)) stop(paste(outfile," already exists!",sep=""))
if (!is.null(kinmat)) {
   trykin <- try(load(kinmat))
   if (inherits(trykin,"try-error")) stop(paste('kinship matrix does not exist at ',kinmat))
}

pheno <- read.table(phenfile,as.is=T,header=T,sep=sep.phe)
if (!phen %in% names(pheno)) stop(paste(phen," does not exist in ",phenfile,"!!",sep=""))
if (!is.null(covars) & sum(covars %in% names(pheno))!=length(covars)) stop(paste("Not all covariates exist in ",phenfile,"!!",sep=""))
if (analysis%in%c("lmepack.int","lmepack.int.imputed","geepack.int","geepack.int.imputed","geepack.quant.int","geepack.quant.int.imputed") && !cov.int%in%covars) stop("covariate for interaction has to be one of the user-specified covariates")
if (analysis%in%c("lmepack.int","lmepack.int.imputed","geepack.int","geepack.int.imputed","geepack.quant.int","geepack.quant.int.imputed") && length(cov.int)!=1) stop("Please specify one covariate for ineteraction!")
if (analysis%in%c("lmepack.int","lmepack.int.imputed","geepack.int","geepack.int.imputed","geepack.quant.int","geepack.quant.int.imputed") && length(model)!=0) stop("No model argument for ineteraction analysis!")

when <- format(Sys.time(), "%Y-%b-%d-%H-%M")

if (analysis=="lmepack") {
   cmds <- character(8)
   if (is.null(model)) model <- "a"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("lmepack.batch(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       if (is.null(covars)) cmds[5] <- paste("kinmat='",kinmat,"',",sep="") else 
           cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),kinmat='",kinmat,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',model='",model,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".lmepack.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".lmepack.",when,".",j,".R ",phen,".lmepack.",when,".",j,".sh',sep=''))",sep="")
       write(cmds,paste(phen,".lmepack.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".lmepack.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".lmepack.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".lmepack.",when,".",j,".R --no-save",sep=""),file=paste(phen,".lmepack.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".lmepack.",when,".",j,".sh",sep=""),paste(phen,".lmepack.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   if (model %in% c("a","d","r")){
  	write(c("phen","snp","n0","n1","n2","h2q","beta","se","chisq","df","model","pval"),outfile,sep=",",ncolumns=12)
   } else write(c("phen","snp","n0","n1","n2","h2q","beta10","beta20","beta21","se10","se20","se21","chisq","df","model","pval"),outfile,ncolumns=16,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".lmepack.",when,".lst",sep=""))

}  else 
if (analysis=="lmepack.imputed") {
   cmds <- character(8)
   if (!is.null(model)) stop(paste("No model option for imputed analysis!!"))
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("lmepack.batch.imputed(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       if (is.null(covars)) cmds[5] <- paste("kinmat='",kinmat,"',",sep="") else 
          cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),kinmat='",kinmat,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".lmepack.imputed.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".lmepack.imputed.",when,".",j,".R ",phen,".lmepack.imputed.",when,".",j,".sh',sep=''))",sep="")
       write(cmds,paste(phen,".lmepack.imputed.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".lmepack.imputed.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".lmepack.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".lmepack.imputed.",when,".",j,".R --no-save",sep=""),file=paste(phen,".lmepack.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".lmepack.imputed.",when,".",j,".sh",sep=""),paste(phen,".lmepack.imputed.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   write(c("phen","snp","N","AF","h2q","beta","se","pval"),outfile,ncolumns=8,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".lmepack.imputed.",when,".lst",sep=""))

} else
if (analysis=="lmeVpack.imputed") {
   cmds <- character(8)
   if (!is.null(model)) stop(paste("No model option for imputed analysis!!"))
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("lmeVpack.batch.imputed(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       if (is.null(covars)) cmds[5] <- paste("kinmat='",kinmat,"',",sep="") else 
          cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),kinmat='",kinmat,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".lmeVpack.imputed.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".lmeVpack.imputed.",when,".",j,".R ",phen,".lmeVpack.imputed.",when,".",j,".sh',sep=''))",sep="")
       write(cmds,paste(phen,".lmeVpack.imputed.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".lmeVpack.imputed.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".lmeVpack.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".lmeVpack.imputed.",when,".",j,".R --no-save",sep=""),file=paste(phen,".lmeVpack.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".lmeVpack.imputed.",when,".",j,".sh",sep=""),paste(phen,".lmeVpack.imputed.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   write(c("phen","snp","N","AF","h2q","beta","se","pval"),outfile,ncolumns=8,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".lmeVpack.imputed.",when,".lst",sep=""))

} else
if (analysis=="geepack") {       
   cmds <- character(8)
   if (is.null(model)) model <- "a"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("geepack.lgst.batch(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       if (is.null(covars)) cmds[5] <- paste("',model='",model,"',",sep="") else 
           cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),model='",model,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".geepack.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".geepack.",when,".",j,".R ",phen,".geepack.",when,".",j,".sh',sep=''))",sep="")
       write(cmds,paste(phen,".geepack.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".geepack.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".geepack.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".geepack.",when,".",j,".R --no-save",sep=""),file=paste(phen,".geepack.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".geepack.",when,".",j,".sh",sep=""),paste(phen,".geepack.",when,".lst",sep=""),ncolumns=1,append=T)
   }
  if (model %in% c("a","d","r")) {  	
	write(c("phen","snp","n0","n1","n2","nd0","nd1","nd2","miss.0","miss.1","miss.diff.p","beta","se","chisq","df","model","remark","pval"),
             outfile,sep=",",ncolumns=18)
  } else
	write(c("phen","snp","n0","n1","n2","nd0","nd1","nd2","miss.0","miss.1","miss.diff.p","beta10","beta20","beta21",
			"se10","se20","se21","chisq","df","model","remark","pval"),outfile,sep=",",ncolumns=22)
   print(paste("Quit R and submit all jobs by using ksh ",phen,".geepack.",when,".lst",sep=""))

}  else 
if (analysis=="glmm") {       
   cmds <- character(8)
   if (is.null(model)) model <- "a"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("glmm.lgst.batch(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       if (is.null(covars)) cmds[5] <- paste("',model='",model,"',",sep="") else 
           cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),model='",model,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".glmm.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".glmm.",when,".",j,".R ",phen,".glmm.",when,".",j,".sh',sep=''))",sep="")
       write(cmds,paste(phen,".glmm.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".glmm.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".glmm.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".glmm.",when,".",j,".R --no-save",sep=""),file=paste(phen,".glmm.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".glmm.",when,".",j,".sh",sep=""),paste(phen,".glmm.",when,".lst",sep=""),ncolumns=1,append=T)
   }
  if (model %in% c("a","d","r")) {  	
	write(c("phen","snp","n0","n1","n2","nd0","nd1","nd2","miss.0","miss.1","miss.diff.p","beta","se","chisq","df","model","remark","pval"),
             outfile,sep=",",ncolumns=18)
  } else
	write(c("phen","snp","n0","n1","n2","nd0","nd1","nd2","miss.0","miss.1","miss.diff.p","beta10","beta20","beta21",
			"se10","se20","se21","chisq","df","model","remark","pval"),outfile,sep=",",ncolumns=22)
   print(paste("Quit R and submit all jobs by using ksh ",phen,".glmm.",when,".lst",sep=""))

}  else 
if (analysis=="geepack.imputed") {       
   cmds <- character(7)
   if (!is.null(model)) stop(paste("No model option for imputed analysis!!"))
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("geepack.lgst.batch.imputed(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F,",sep="")
       if (is.null(covars)) cmds[5] <- paste("phen='",phen,"')",sep="") else 
           cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),phen='",phen,"')",sep="")                   
       cmds[6] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".geepack.imputed.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[7] <- paste("system(paste('rm ",phen,".geepack.imputed.",when,".",j,".R ",phen,".geepack.imputed.",when,".",j,".sh',sep=''))",sep="") 
       write(cmds,paste(phen,".geepack.imputed.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".geepack.imputed.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".geepack.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".geepack.imputed.",when,".",j,".R --no-save",sep=""),file=paste(phen,".geepack.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".geepack.imputed.",when,".",j,".sh",sep=""),paste(phen,".geepack.imputed.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   write(c("phen","snp","N","Nd","AF","AFd","beta","se","remark","pval"),outfile,ncolumns=10,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".geepack.imputed.",when,".lst",sep=""))
}  else
if (analysis=="geepack.quant") {       
   cmds <- character(8)
   if (is.null(model)) model <- "a"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("geepack.quant.batch(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       if (is.null(covars)) cmds[5] <- paste("model='",model,"',",sep="") else 
           cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),model='",model,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".geepack.quant.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".geepack.quant.",when,".",j,".R ",phen,".geepack.quant.",when,".",j,".sh',sep=''))",sep="")
       write(cmds,paste(phen,".geepack.quant.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".geepack.quant.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".geepack.quant.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".geepack.quant.",when,".",j,".R --no-save",sep=""),file=paste(phen,".geepack.quant.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".geepack.quant.",when,".",j,".sh",sep=""),paste(phen,".geepack.quant.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   if (model %in% c("a","d","r")){
  	write(c("phen","snp","n0","n1","n2","beta","se","chisq","df","model","pval"),outfile,sep=",",ncolumns=11)
   } else write(c("phen","snp","n0","n1","n2","beta10","beta20","beta21","se10","se20","se21","chisq","df","model","pval"),outfile,ncolumns=15,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".geepack.quant.",when,".lst",sep=""))

}  else 
if (analysis=="geepack.quant.imputed") {       
   cmds <- character(7)
   if (!is.null(model)) stop(paste("No model option for imputed analysis!!"))
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("geepack.quant.batch.imputed(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F,",sep="")
       if (is.null(covars)) cmds[5] <- paste("phen='",phen,"')",sep="") else 
           cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),phen='",phen,"')",sep="")                   
       cmds[6] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".geepack.quant.imputed.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[7] <- paste("system(paste('rm ",phen,".geepack.quant.imputed.",when,".",j,".R ",phen,".geepack.quant.imputed.",when,".",j,".sh',sep=''))",sep="") 
       write(cmds,paste(phen,".geepack.quant.imputed.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".geepack.quant.imputed.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".geepack.quant.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".geepack.quant.imputed.",when,".",j,".R --no-save",sep=""),file=paste(phen,".geepack.quant.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".geepack.quant.imputed.",when,".",j,".sh",sep=""),paste(phen,".geepack.quant.imputed.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   write(c("phen","snp","N","AF","beta","se","pval"),outfile,ncolumns=7,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".geepack.quant.imputed.",when,".lst",sep=""))
}  else
if (analysis=="lmepack.int") {
   cmds <- character(8)
   if (is.null(sub)) sub <- "N"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("lmepack.int.batch(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),kinmat='",kinmat,"',cov.int='",cov.int,"',sub='",sub,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".lmepack.int.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".lmepack.int.",when,".",j,".R ",phen,".lmepack.int.",when,".",j,".sh',sep=''))",sep="")
       write(cmds,paste(phen,".lmepack.int.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".lmepack.int.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".lmepack.int.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".lmepack.int.",when,".",j,".R --no-save",sep=""),file=paste(phen,".lmepack.int.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".lmepack.int.",when,".",j,".sh",sep=""),paste(phen,".lmepack.int.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   if (sub=="N"){
  	write(c("phen","snp","covar_int","n","AF","cov_beta_snp_beta_int","model","beta_snp","se_snp","pval_snp","beta_int","se_int","pval_int"),outfile,sep=",",ncolumns=13)
   } else write(c("phen","snp","covar_int","n","AF","model","beta_snp","se_snp","pval_snp","beta_snp_cov0","se_snp_cov0","pval_snp_cov0","beta_snp_cov1","se_snp_cov1","pval_snp_cov1","beta_int","se_int","pval_int"),outfile,ncolumns=18,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".lmepack.int.",when,".lst",sep=""))

}  else 
if (analysis=="lmepack.int.imputed") {
   cmds <- character(8)
   if (is.null(sub)) sub <- "N"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("lmepack.int.batch.imputed(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),kinmat='",kinmat,"',cov.int='",cov.int,"',sub='",sub,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".lmepack.int.imputed.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".lmepack.int.imputed.",when,".",j,".R ",phen,".lmepack.int.imputed.",when,".",j,".sh',sep=''))",sep="")
       write(cmds,paste(phen,".lmepack.int.imputed.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".lmepack.int.imputed.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".lmepack.int.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".lmepack.int.imputed.",when,".",j,".R --no-save",sep=""),file=paste(phen,".lmepack.int.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".lmepack.int.imputed.",when,".",j,".sh",sep=""),paste(phen,".lmepack.int.imputed.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   if (sub=="N"){
  	write(c("phen","snp","covar_int","n","AF","cov_beta_snp_beta_int","model","beta_snp","se_snp","pval_snp","beta_int","se_int","pval_int"),outfile,sep=",",ncolumns=13)
   } else write(c("phen","snp","covar_int","n","AF","model","beta_snp","se_snp","pval_snp","beta_snp_cov0","se_snp_cov0","pval_snp_cov0","beta_snp_cov1","se_snp_cov1","pval_snp_cov1","beta_int","se_int","pval_int"),outfile,ncolumns=18,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".lmepack.int.imputed.",when,".lst",sep=""))

} else
if (analysis=="geepack.int") {       
   cmds <- character(8)
   if (is.null(sub)) sub <- "N"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("geepack.lgst.int.batch(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),cov.int='",cov.int,"',sub='",sub,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".geepack.int.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".geepack.int.",when,".",j,".R ",phen,".geepack.int.",when,".",j,".sh',sep=''))",sep="")
       write(cmds,paste(phen,".geepack.int.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".geepack.int.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".geepack.int.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".geepack.int.",when,".",j,".R --no-save",sep=""),file=paste(phen,".geepack.int.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".geepack.int.",when,".",j,".sh",sep=""),paste(phen,".geepack.int.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   if (sub=="N"){
  	write(c("phen","snp","covar_int","n","AF","nd","AFd","cov_beta_snp_beta_int","model","beta_snp","se_snp","pval_snp","beta_int","se_int","remark","pval_int"),outfile,sep=",",ncolumns=16)
   } else write(c("phen","snp","covar_int","n","AF","nd","AFd","model","beta_snp","se_snp","pval_snp","beta_snp_cov0","se_snp_cov0","pval_snp_cov0","beta_snp_cov1","se_snp_cov1","pval_snp_cov1","beta_int","se_int","remark","pval_int"),outfile,ncolumns=21,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".geepack.int.",when,".lst",sep=""))

}  else 
if (analysis=="geepack.int.imputed") {       
   cmds <- character(7)
   if (is.null(sub)) sub <- "N"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("geepack.lgst.int.batch.imputed(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F,",sep="")
       cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),phen='",phen,"',cov.int='",cov.int,"',sub='",sub,"')",sep="")                   
       cmds[6] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".geepack.int.imputed.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[7] <- paste("system(paste('rm ",phen,".geepack.int.imputed.",when,".",j,".R ",phen,".geepack.int.imputed.",when,".",j,".sh',sep=''))",sep="") 
       write(cmds,paste(phen,".geepack.int.imputed.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".geepack.int.imputed.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".geepack.int.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".geepack.int.imputed.",when,".",j,".R --no-save",sep=""),file=paste(phen,".geepack.int.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".geepack.int.imputed.",when,".",j,".sh",sep=""),paste(phen,".geepack.int.imputed.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   if (sub=="N"){
  	write(c("phen","snp","covar_int","n","AF","nd","AFd","cov_beta_snp_beta_int","model","beta_snp","se_snp","pval_snp","beta_int","se_int","remark","pval_int"),outfile,sep=",",ncolumns=16)
   } else write(c("phen","snp","covar_int","n","AF","nd","AFd","model","beta_snp","se_snp","pval_snp","beta_snp_cov0","se_snp_cov0","pval_snp_cov0","beta_snp_cov1","se_snp_cov1","pval_snp_cov1","beta_int","se_int","remark","pval_int"),outfile,ncolumns=21,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".geepack.int.imputed.",when,".lst",sep=""))
}  else
if (analysis=="geepack.quant.int") {       
   cmds <- character(8)
   if (is.null(sub)) sub <- "N"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("geepack.quant.int.batch(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),cov.int='",cov.int,"',sub='",sub,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".geepack.quant.int.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".geepack.quant.int.",when,".",j,".R ",phen,".geepack.quant.int.",when,".",j,".sh',sep=''))",sep="")
       write(cmds,paste(phen,".geepack.quant.int.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".geepack.quant.int.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".geepack.quant.int.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".geepack.quant.int.",when,".",j,".R --no-save",sep=""),file=paste(phen,".geepack.quant.int.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".geepack.quant.int.",when,".",j,".sh",sep=""),paste(phen,".geepack.quant.int.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   if (sub=="N"){
  	write(c("phen","snp","covar_int","n","AF","cov_beta_snp_beta_int","model","beta_snp","se_snp","pval_snp","beta_int","se_int","pval_int"),outfile,sep=",",ncolumns=13)
   } else write(c("phen","snp","covar_int","n","AF","model","beta_snp","se_snp","pval_snp","beta_snp_cov0","se_snp_cov0","pval_snp_cov0","beta_snp_cov1","se_snp_cov1","pval_snp_cov1","beta_int","se_int","pval_int"),outfile,ncolumns=18,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".geepack.quant.int.",when,".lst",sep=""))

}  else 
if (analysis=="geepack.quant.int.imputed") {       
   cmds <- character(7)
   if (is.null(sub)) sub <- "N"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("geepack.quant.int.batch.imputed(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F,",sep="")
       cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),phen='",phen,"',cov.int='",cov.int,"',sub='",sub,"')",sep="")                   
       cmds[6] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".geepack.quant.int.imputed.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[7] <- paste("system(paste('rm ",phen,".geepack.quant.int.imputed.",when,".",j,".R ",phen,".geepack.quant.int.imputed.",when,".",j,".sh',sep=''))",sep="") 
       write(cmds,paste(phen,".geepack.quant.int.imputed.",when,".",j,".R",sep=""),ncolumns=1)
       write("#$ -o /dev/null",file=paste(phen,".geepack.quant.int.imputed.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".geepack.quant.int.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".geepack.quant.int.imputed.",when,".",j,".R --no-save",sep=""),file=paste(phen,".geepack.quant.int.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".geepack.quant.int.imputed.",when,".",j,".sh",sep=""),paste(phen,".geepack.quant.int.imputed.",when,".lst",sep=""),ncolumns=1,append=T)
   }
   if (sub=="N"){
  	write(c("phen","snp","covar_int","n","AF","cov_beta_snp_beta_int","model","beta_snp","se_snp","pval_snp","beta_int","se_int","pval_int"),outfile,sep=",",ncolumns=13)
   } else write(c("phen","snp","covar_int","n","AF","model","beta_snp","se_snp","pval_snp","beta_snp_cov0","se_snp_cov0","pval_snp_cov0","beta_snp_cov1","se_snp_cov1","pval_snp_cov1","beta_int","se_int","pval_int"),outfile,ncolumns=18,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".geepack.quant.int.imputed.",when,".lst",sep=""))

}

}