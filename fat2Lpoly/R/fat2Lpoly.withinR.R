
# Jordie Croteau
# 23 aout 2012

# correspond au fichier fonction_fat2Lpoly.withinR_v2.R dans le dossier "programmes"


# cette fonction est identique a l'ancienne version de fat2Lpoly (fichier fonction_globale_donnees_SPAP_v2.R), excepte qu'au lieu d'appeler 
# read.merlin.files, fat2Lpoly.withinR prend comme argument un objet dont le format est celui de la sortie de la fonction read.merlin.files.
# De plus, fat2Lpoly.withinR s'arrete juste avant l'appel de get.scores.pvalues dans l'ancienne version de fat2Lpoly.

# petites modifs apportees le 28 aout 2012 par Jordie.

# ajout de code enlevant les familles n'ayant pas plus d'une categorie representee, le 17 octobre 2012

# 4 avril 2013: ajout du calcul des coefficients de kinship (a priori) au lieu des IBD, dans le cas ou cette derniere information n'est pas fournie.

fat2Lpoly.withinR=function(ped.x.all,snp.names.mat,ibd.loci=NULL,contingency.file=FALSE,design.constraint,par.constrained,constraints,pairweights=calcule.poids.alphafixe,lc=NULL,alpha=NULL)
{
ped=ped.x.all$ped
x.all=ped.x.all$x.all
c.names=colnames(x.all) 
ibd.dat.list=ped.x.all$ibd.dat.list
n.loc=length(ibd.dat.list)
ibdfilenames=ped.x.all$ibdfilenames

res=vector(mode="list",length=nrow(snp.names.mat))

if(contingency.file)
 {
  date.heure=Sys.time()
  date.heure.mod=paste(substr(date.heure,1,10),"_",substr(date.heure,12,13),".",substr(date.heure,15,16),".",substr(date.heure,18,19),sep="")
  descrip.file=paste("descriptive_statistics",date.heure.mod,".txt",sep="")
  cat("some descriptive statistics will progressively be added to the file ",descrip.file,"\n\n")
 }
  
# si les IBD ne sont pas fournies, on calcule le kinship a priori de toutes les paires de sujets a l'interieur de chaque famille
# et on met le tout dans le meme format que pour les IBD a posteriori.
if(is.null(ibd.loci)|is.null(ibd.dat.list[[1]]))
 {
  pere=ped.x.all$pere
  mere=ped.x.all$mere
  
  ibd.dat=NULL
  
  fam.u=unique(ped[,1])
  for(j in 1:length(fam.u))
   {
    indices=ped[,1]==fam.u[j]
	# traiter seulement les familles contenant plus d'un sujet
	if(sum(indices)>1)
	 {
      ped.tmp=ped[indices,]
	  sujet.tmp=ped.tmp[,2]
	  pere.tmp=pere[indices]
	  mere.tmp=mere[indices]

	  # pour que le calcul de kinship soit valide, il faut que la structure de famille soit complete:
	  # par consequent, on exige que chaque sujet satisfasse une des 2 conditions suivantes:
	  # 1- parent d'un autre sujet present dans la meme famille ou
	  # 2- enfant de 2 parents inclus dans la famille (les 2).
	  if(!all(sujet.tmp%in%pere.tmp|sujet.tmp%in%mere.tmp|(pere.tmp%in%sujet.tmp&mere.tmp%in%sujet.tmp))) 
	  stop(paste("When IBD is not provided, pedigree structures must be complete. Family",fam.u[j],"is not complete."))
	  
      matk<-2*kinship(sujet.tmp,pere.tmp,mere.tmp)
	
      pi.tmp=as.numeric(matk[lower.tri(matk)])
	  pi.all=pi.tmp
	  if(n.loc>1) for(i in 1:(n.loc-1)) pi.all=cbind(pi.all,pi.tmp)
	
      paires=outer(sujet.tmp,sujet.tmp,paste)
      kin.tmp=data.frame(matrix(unlist(strsplit(paires[lower.tri(paires)],split=" ")),ncol=2,byrow=T),pi.all)
  
      kin.diag=as.numeric(diag(matk))
  	  if(n.loc>1) for(i in 1:(n.loc-1)) kin.diag=cbind(kin.diag,as.numeric(diag(matk)))

	  kin.diag=data.frame(sujet.tmp,sujet.tmp,kin.diag)
	  colnames(kin.diag)=colnames(kin.tmp)
      kin.tmp=rbind(kin.diag,kin.tmp)
      ibd.dat=rbind(ibd.dat,cbind(rep(fam.u[j],nrow(kin.tmp)),kin.tmp))
	 }
   }
  colnames(ibd.dat)[1:3]=c("FAMILY","ID1","ID2")
 }   
  
# Boucle principale sur les marqueurs testes commence ici
for(s in 1:nrow(snp.names.mat)){

snp.names=snp.names.mat[s,]

if(length(snp.names)==1) cat(paste("analyzing SNP",snp.names,"\n"))
if(length(snp.names)==2) cat(paste("analyzing SNP pair",paste(snp.names,collapse=" ")),"\n")

if(contingency.file) cat(snp.names,"\n",file=descrip.file,append=(s!=1))

# si les IBD sont fournies, on extrait celles qu'il nous faut parmi celles-ci.
if(!is.null(ibd.loci)&!is.null(ibd.dat.list[[1]]))
 {
  ibd.loci.tmp=ibd.loci[s,]
  if(contingency.file)
   {
    cat(ibd.loci.tmp,"\n",file=descrip.file,append=TRUE)
    cat("\n",file=descrip.file,append=TRUE)
   }

  ##################### combiner les IBD des marqueurs (ou positions de marqueurs) de ibd.loci.tmp ###################################
  ibd.dat=ibd.dat.list[[1]]
  if(sum((ibd.dat$MARKER)==ibd.loci.tmp[1])==0) stop(paste(ibd.loci.tmp[1],"is absent from",ibdfilenames[1]))
  ibd.dat=ibd.dat[(ibd.dat$MARKER)==ibd.loci.tmp[1],]
  pi1=(ibd.dat$P1)/2+ibd.dat$P2
  ibd.dat=cbind(ibd.dat[,1:3],pi1)

  if(n.loc>1)
   {
    for(loc.num in 2:n.loc)
     {
      ibd.dat.tmp=ibd.dat.list[[loc.num]]
	  if(sum((ibd.dat.tmp$MARKER)==ibd.loci.tmp[loc.num])==0) stop(paste(ibd.loci.tmp[loc.num],"is absent from",ibdfilenames[loc.num]))

      ibd.dat.tmp=ibd.dat.tmp[(ibd.dat.tmp$MARKER)==ibd.loci.tmp[loc.num],]
      pi.tmp=(ibd.dat.tmp$P1)/2+ibd.dat.tmp$P2
      ibd.dat.tmp=cbind(ibd.dat.tmp[,1:3],pi.tmp)
      ibd.dat=merge(ibd.dat,ibd.dat.tmp,by=c("FAMILY","ID1","ID2"),all.x=TRUE,all.y=TRUE,sort=FALSE)
     }
   }
  ######################################################################################################################################
 }

# pour le moment on suppose que snp.names.mat a une ou 2 colonnes (il y a un ou 2 locus)
ifelse(ncol(snp.names.mat)==1,x <- as.matrix(x.all[,c.names==snp.names]),x <- x.all[,c(which(c.names==snp.names[1]),which(c.names==snp.names[2]))])

# extraire n.levels de design.constraint a l'aide d'un x bidon.
n.levels=attributes(design.constraint(array(1,c(2,ncol(x))),par.constrained=par.constrained,constraints=constraints))$n.levels

if(n.levels==4){
# enlever les sujets pour lesquels soit y1 est manquant, ou y2 est manquant, ou un des genotypes des differents locus est manquant.
ped.filtre=ped[ped$y1!=0&ped$y2!=0&!apply(is.na(x),1,any),]
fam.id=ped.filtre[,1]
subject.ids=ped.filtre[,2]
y1=ped.filtre$y1
y2=ped.filtre$y2
dims=c(sum(ped$y1!=0&ped$y2!=0&!apply(is.na(x),1,any)),ncol(x))
x=array(x[ped$y1!=0&ped$y2!=0&!apply(is.na(x),1,any),],dims)

# creation du vecteur de reponse polytomique a partir de y1 et y2 
y=rep(NA,length(y1))
y[y1==1&y2==1]=4
y[y1==1&y2==2]=2
y[y1==2&y2==1]=1
y[y1==2&y2==2]=3
y=factor(y,levels=1:n.levels)

# enlever les familles n'ayant pas plus d'une categorie representee et emettre un avertissement s'il y a lieu.
nb.cat.par.fam=tapply(y,fam.id,function(x) length(unique(x)))
if(any(nb.cat.par.fam<=1))
 {
  fams.enleve=as.character(names(nb.cat.par.fam)[nb.cat.par.fam<=1])
  ped.filtre=ped.filtre[!(as.character(fam.id)%in%fams.enleve),]
  dims=c(nrow(ped.filtre),ncol(x))
  x=array(x[!(as.character(fam.id)%in%fams.enleve),],dims)
  y=y[!(as.character(fam.id)%in%fams.enleve)]  
  fam.id=ped.filtre[,1]
  subject.ids=ped.filtre[,2]

  cat(paste("The families",paste(fams.enleve,collapse=", "),"\n","are uninformative and have been excluded.  For a family to be informative, its genotyped subjects must not all belong to the same category (phenotype combination)."),"\n\n")
 }

if(contingency.file)
 {
  cat("frequency table\n",file=descrip.file,append=TRUE)
  cat(c("Y1=2,Y2=1","Y1=1,Y2=2","Y1=2,Y2=2","Y1=1,Y2=1"),"\n",file=descrip.file,append=TRUE)
  cat(table(y),"\n",file=descrip.file,append=TRUE)
  cat("\n",file=descrip.file,append=TRUE)
 }
 
if(any(as.numeric(table(y))==0)){
 categ=which(table(y)==0)
 if(categ==1) categ="Y1=2, Y2=1"
 else if(categ==2) categ="Y1=1, Y2=2"
 else if(categ==3) categ="Y1=2, Y2=2"
 else if(categ==4) categ="Y1=1, Y2=1"
 stop(paste("There is no subject with level",categ,"in any of the families (each level must be represented in at least one family)"))
 }
 
# JC, 9 juillet 2013: s'il y a 2 locus, estimer les alphas par regression polytomique GEE de y sur X du locus lc.
if(!is.null(alpha)&is.null(lc)) stop("If alpha is not null, locus number on which to condition must be given")

if(ncol(x)==2&!is.null(lc)&is.null(alpha))
 {
  if (!(lc %in% 1:ncol(x))) stop("The index lc does not correspond to a valid locus index in",1:ncol(x))
  dat=data.frame(fam.id,subject.ids,y,as.numeric(x[,lc]))
  colnames(dat)=c("fam.id","sub","y","geno")
  un.loc=nomLORgee(y~geno,data=dat,id=fam.id,repeated=sub,LORstr = "independence")
  coefs.all=un.loc$coef
  alpha=coefs.all[c(2,4,6)]
 }
}

else if(n.levels==2){
cat("Analysis with only one dichotomous phenotype","\n")
cat("Phenotype used is",ped.x.all$y1.name,"\n")

# enlever les sujets pour lesquels soit le phenotype est manquant ou un des genotypes des differents locus est manquant.
ped.filtre=ped[ped$y1!=0&!apply(is.na(x),1,any),]
fam.id=ped.filtre[,1]
subject.ids=ped.filtre[,2]
y=factor(3-ped.filtre$y1,levels=1:n.levels)
dims=c(sum(ped$y1!=0&!apply(is.na(x),1,any)),ncol(x))
x=array(x[ped$y1!=0&!apply(is.na(x),1,any),],dims)

# enlever les familles n'ayant pas plus d'une categorie representee et emettre un avertissement s'il y a lieu.
nb.cat.par.fam=tapply(y,fam.id,function(x) length(unique(x)))
if(any(nb.cat.par.fam<=1))
 {
  fams.enleve=as.character(names(nb.cat.par.fam)[nb.cat.par.fam<=1])
  ped.filtre=ped.filtre[!(as.character(fam.id)%in%fams.enleve),]
  dims=c(nrow(ped.filtre),ncol(x))
  x=array(x[!(as.character(fam.id)%in%fams.enleve),],dims)
  y=y[!(as.character(fam.id)%in%fams.enleve)]
  fam.id=ped.filtre[,1]
  subject.ids=ped.filtre[,2]

  cat(paste("The families",paste(fams.enleve,collapse=", "),"are uninformative and have been excluded.  For a family to be informative, both categories must be represented."),"\n\n")
 }

if(contingency.file)
 {
  cat("frequency table\n",file=descrip.file,append=TRUE)
  cat(c("Y2=2","Y2=1"),"\n",file=descrip.file,append=TRUE)
  cat(table(y),"\n",file=descrip.file,append=TRUE)
  cat("\n",file=descrip.file,append=TRUE)
 }
 
 # JC, 11 juillet 2013: s'il y a 2 locus, estimer alpha par regression logistique de y sur X du locus lc.
if(!is.null(alpha)&is.null(lc)) stop("If alpha is not null, locus number on which to condition must be given")

if(ncol(x)==2&!is.null(lc)&is.null(alpha))
 {
  if (!(lc %in% 1:ncol(x))) stop("The index lc does not correspond to a valid locus index in",1:ncol(x))
  dat=data.frame(2-as.numeric(y),as.numeric(x[,lc]))
  colnames(dat)=c("y","geno")
  glm.fit=glm(y~geno,data=dat,family=binomial(link = "logit"))
  alpha=(summary(glm.fit)$coef)[2,1]
 }
}
else stop("attribute n.levels of the design.constraint function must be 2 or 4")

# creation des matrices individuelles pour chaque categorie k de 1 a K-1 a partir de la matrice des "0,0.5,1" a l'aide d'une fonction specifique
tmp <- design.constraint(x,par.constrained=par.constrained,constraints=constraints)

# creation de la matrice de design xp pour le calcul du score 
tmp2=design.polytomous(x=tmp$x.e,K=n.levels,x.loc.in=tmp$x.loc.e,par.constrained=par.constrained,constraints=constraints)
xp = tmp2$xp
# Liste du vecteur de locus impliques dans chaque parametre
xp.loc = tmp2$x.loc.out

# Creation de la matrice de design xl pour le calcul des covariances (On prend par.constrainted et constraints ressortis par design.constraint)
# rep.par est le vecteur du nombre de parametres dans chaque matrice x
rep.par <- unlist(lapply(tmp$x.e,ncol))
if (missing(constraints)) tmp2=design.polytomous(x=tmp$x.l,K=n.levels,x.loc.in=tmp$x.loc.l,rep.par=rep.par)
else tmp2=design.polytomous(x=tmp$x.l,K=n.levels,x.loc.in=tmp$x.loc.l,par.constrained=tmp$par.constrained,constraints=tmp$constraints,rep.par=rep.par)
xl = tmp2$xp
# a partir de version 15, xl.loc est caractere au lieu de numerique
xl.loc = unlist(tmp2$x.loc.out)
# Conversion des termes de produits en indices 
il <- strsplit(xl.loc,split="")
xibd.loc= sapply(il,converti.terme,n.loc=n.loc)

# ind.par nous donne les indices des locus pour la categorie a laquelle chaque terme appartient
ind.par = tmp2$ind.par
# ind.cat nous donne la categorie a laquelle appartient chaque parametre
# et ind.catl la categorie a laquelle appartient chaque terme
ind.catl = rep(1:(n.levels-1),lapply(unique(ind.par),length))
ind.cat = rep(1:(n.levels-1),rep.par)

res[[s]]=scores.covs(subject.ids,fam.id,y,n.levels,ibd.dat,n.loc,xp,xp.loc,xl,il,xibd.loc,ind.par,rep.par,ind.catl,ind.cat,contingency.file,descrip.file,calculpoids=pairweights,lc=lc,alpha.vec=alpha)
}

list(scores.covs.all.SNPs=res,snp.names.mat=snp.names.mat)
}
