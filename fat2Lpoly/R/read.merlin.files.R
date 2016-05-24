
# Jordie Croteau
# 17 aout 2012

# correspond au fichier read_merlin_files_v3.R dans le dossier "programmes"

# modifie le 21 aout pour avoir la possibilite de prendre comme argument des vecteurs d'alleles mineures au lieu de fichiers de frequences d'alleles.

# modifie le 23 aout pour integrer la lecture des donnees d'IBD.

# 4 avril 2013: valeur NULL par defaut ajoutee pour l'argument ibdfilenames (dans ce cas, aucune lecture des donnees d'IBD)

read.merlin.files=function(pedfilenames,datfilenames,freq.data,ibdfilenames=NULL)
{
###################### Definition des arguments #####################################################################################
# pedfilenames : vecteur des noms de fichiers ped (un fichier par locus).  Les sujets inclus peuvent etre un sous-ensemble de ceux 
#                inclus dans les fichiers d'IDB.
# datfilenames : vecteur des noms de fichiers dat (un fichier par locus). 
# freqfilenames: vecteur des noms de fichiers freq (un fichier par locus). 
# ibdfilenames: 
# Tous ces fichiers doivent etre en format Merlin (voir http://www.sph.umich.edu/csg/abecasis/Merlin/tour/input_files.html pour la description detaillee de ce format)
#####################################################################################################################################

############### lecture des fichiers datfile pour obtenir les noms de SNPs et de phenotype #######################
dat1=read.table(datfilenames[1],as.is=TRUE)
if(dat1[1,1]!="A") stop(paste("The first line of",datfilenames[1],"should be starting by 'A'."))

# si 2 noms d'affection status sont fournis, on suppose que le 1er est l'endophenotype et le 2e, le phenotype.
if(sum(dat1[,1]=="A")==2){ 
nb.pheno=2
pheno.name=dat1[2,2]
endo.name=dat1[1,2]
cat("\n")
cat("Y1 data extracted from input files:   ",endo.name,"\n")
cat("Y2 data extracted from input files:   ",pheno.name,"\n")
cat("\n")
if(!all(dat1[3:nrow(dat1),1]=="M")) stop(paste("Lines number 3 to",nrow(dat1),"of",datfilenames[1],"should all be starting by 'M'."))
}
# si seulement 1 nom d'affection status est fourni, on met pheno.name=NULL et on travaille avec endo.name (bien que ce dernier puisse representer un phenotype (pas endo))
else{
nb.pheno=1
endo.name=dat1[1,2]
pheno.name=NULL
cat("\n")
cat("Analysis will be for dichotomous phenotype:  ",endo.name,"\n")
cat("\n")
if(!all(dat1[2:nrow(dat1),1]=="M")) stop(paste("Lines number 2 to",nrow(dat1),"of",datfilenames[1],"should all be starting by 'M'."))
}

snp.names.dat=dat1[dat1[,1]=="M",2]

n.loc=length(pedfilenames)
if(n.loc>1)
 {
  for(loc.num in 2:n.loc)
   {
    dat.tmp=read.table(datfilenames[loc.num],as.is=TRUE)
	if(dat.tmp[1,1]!="A") stop(paste("The first line of",datfilenames[loc.num],"should be starting by 'A'."))
	if(nb.pheno==2){ if(!all(dat.tmp[3:nrow(dat.tmp),1]=="M")) stop(paste("Lines number 3 to",nrow(dat.tmp),"of",datfilenames[loc.num],"should all be starting by 'M'."))}
    else{ if(!all(dat.tmp[2:nrow(dat.tmp),1]=="M")) stop(paste("Lines number 2 to",nrow(dat.tmp),"of",datfilenames[loc.num],"should all be starting by 'M'."))}
    snp.names.dat=c(snp.names.dat,dat.tmp[dat.tmp[,1]=="M",2])
   }
 }
################################################################################################################################


# extraction des fam.id, subject.ids, y1 et y2, a partir du premier fichier ped.
ped1.tmp=read.table(pedfilenames[1],header=FALSE,as.is=TRUE)
if(any(apply(ped1.tmp,2,is.character))) stop(paste(pedfilenames[1],"contains letters in some fields.  All fields must be numeric."))
fam.id=ped1.tmp[,1]
subject.ids=ped1.tmp[,2]
if(nb.pheno==2){ # 2 affection status fournis
y1=ped1.tmp[,6]
y2=ped1.tmp[,7]
y1[is.na(y1)]=0
y2[is.na(y2)]=0
}
else{ # 1 seul affection status fourni
y1=ped1.tmp[,6]
y1[is.na(y1)]=0
}

################ extraction des genotypes de tous les locus #####################################
if(!is.null(pheno.name)) ped=data.frame(fam.id,subject.ids,ped1.tmp[,3:4],y1,y2,ped1.tmp[,8:ncol(ped1.tmp)])
else ped=data.frame(fam.id,subject.ids,ped1.tmp[,3:4],y1,ped1.tmp[,7:ncol(ped1.tmp)])

if(n.loc>1)
 {
  for(loc.num in 2:n.loc)
   {
    ped.tmp=read.table(pedfilenames[loc.num],header=FALSE,as.is=TRUE)
	if(any(apply(ped.tmp,2,is.character))) stop(paste(pedfilenames[loc.num],"contains letters in some fields.  All fields must be numeric."))
	ped.tmp=ped.tmp[,c(1,2,(6+nb.pheno):ncol(ped.tmp))] 
	colnames(ped.tmp)[1:2]=c("fam.id","subject.ids")
	ped=merge(ped,ped.tmp,by=c("fam.id","subject.ids"),all.x=FALSE,all.y=FALSE,sort=FALSE)
   }
 }
pere=ped[,3]
mere=ped[,4]
ped=ped[,-(3:4)]
if(nrow(ped)<nrow(ped1.tmp)) warning(paste("Subjects from 1 or more ped files differ from those of other ped files. Only subjects found in the",n.loc,"ped files at the same time are kept for analysis."))
##################################################################################################
 

################################################ obtention des alleles mineures ##########################################################
if(is.list(freq.data)) MA=unlist(freq.data)
else
{ 
######################### lecture des fichiers freq pour obtenir les alleles mineures  ###########################################

# il y a 2 formats de fichiers freq possibles dans merlin (avec frequences d'alleles une sous l'autre ou une a cote de l'autre).
# Donc, pour chaque fichier freq, il faut d'abord detecter quel est le format.
freq=readLines(freq.data[1])
n.lines=length(freq)
greps.M=grep("M",freq)
format.freq=as.numeric(length(greps.M)==(n.lines/2))

name.lines=freq[greps.M]
write(name.lines,"name.lines_asdhskjhfdak.txt")
snp.names.freq=as.character(read.table("name.lines_asdhskjhfdak.txt",as.is=TRUE)[,2])
unlink("name.lines_asdhskjhfdak.txt")

if(format.freq==0)
 {
  freq.col=read.table(freq.data[1],as.is=TRUE)[,2]
  freqs=vector(mode="list",length=length(greps.M)) 
  greps.M.mod=c(greps.M,length(freq.col)+1)
  for(j in 1:(length(greps.M.mod)-1)) freqs[[j]]=as.numeric(freq.col[(greps.M.mod[j]+1):(greps.M.mod[j+1]-1)])
 }
 
if(format.freq==1) freqs=lapply(strsplit(freq[-greps.M],split=" "),function(x){ x=as.numeric(x); return(x[!is.na(x)])})

# alleles mineures
MA=unlist(lapply(freqs,function(x) {y=x[x>0]; ifelse(length(y)==2,which(x==min(y)),NA)}))
# alleles majeures
MjA=unlist(lapply(freqs,which.max))

if(n.loc>1)
 {
  for(loc.num in 2:n.loc)
   {
    freq=readLines(freq.data[loc.num])
    n.lines=length(freq)
    greps.M=grep("M",freq)
    format.freq=as.numeric(length(greps.M)==(n.lines/2))

    name.lines=freq[greps.M]
    write(name.lines,"name.lines_asdhskjhfdak.txt")
    snp.names.freq.tmp=as.character(read.table("name.lines_asdhskjhfdak.txt",as.is=TRUE)[,2])
    snp.names.freq=c(snp.names.freq,snp.names.freq.tmp)

    unlink("name.lines_asdhskjhfdak.txt")

    if(format.freq==0)
     {
      freq.col=read.table(freq.data[loc.num],as.is=TRUE)[,2]
      freqs=vector(mode="list",length=length(greps.M)) 
      greps.M.mod=c(greps.M,length(freq.col)+1)
      for(j in 1:(length(greps.M.mod)-1)) freqs[[j]]=as.numeric(freq.col[(greps.M.mod[j]+1):(greps.M.mod[j+1]-1)])
     }
 
    if(format.freq==1) freqs=lapply(strsplit(freq[-greps.M],split=" "),function(x){ x=as.numeric(x); return(x[!is.na(x)])})
	
	# alleles mineures
    MA=c(MA,unlist(lapply(freqs,function(x) {y=x[x>0]; ifelse(length(y)==2,which(x==min(y)),NA)})))
    # alleles majeures
    MjA=c(MjA,unlist(lapply(freqs,which.max)))
   }
 }
if(any(is.na(MA))) stop("one or more SNP(s) with unknown minor allele, probably because of nul minor allele frequency(ies)")
 
if(!all(snp.names.freq==snp.names.dat)) stop("SNP names or order in freq files not the same as in the data files")

# verifier que toutes les colonnes de genotypes ne contiennent pas d'autres alleles que les alleles mineures et majeures obtenues 
# a partir des fichiers freq.
for(j in 1:length(MA))
 {
  genos.tmp=as.vector(ped[,c(1+nb.pheno+2*j,2+nb.pheno+2*j)])
  if(!all(is.na(genos.tmp)|genos.tmp==0))
   {
    genos.tmp=genos.tmp[!(is.na(genos.tmp)|genos.tmp==0)]
    if(sum(!(genos.tmp%in%c(MA[j],MjA[j])))!=0) stop(paste("one or more alleles differ from minor and major alleles expected from freq files, for SNP",snp.names.dat[j]))
   }
 }
################################################################################################################################
}

# conversion des genotypes en valeurs 0,0.5,1 (mode="allelic")
x.all=alleles2sums(ped[,(3+nb.pheno):ncol(ped)],MA.vec=MA,snp.names=snp.names.dat,mode="allelic")

MA.table=data.frame(snp.names.dat,MA)

# lecture des fichiers d'IBD
ibd.dat.list=vector(mode="list",length=n.loc)
if(!is.null(ibdfilenames)) for(loc.num in 1:n.loc) ibd.dat.list[[loc.num]]=as.data.frame(read.table(ibdfilenames[loc.num],header=TRUE,as.is=TRUE))

list(ped=ped[,1:(2+nb.pheno)],x.all=x.all,MA.table=MA.table,ibd.dat.list=ibd.dat.list,y1.name=endo.name,y2.name=pheno.name,ibdfilenames=ibdfilenames,pere=pere,mere=mere)
}                      