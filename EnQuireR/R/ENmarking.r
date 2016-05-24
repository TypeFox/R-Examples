#MARQUAGE SEMANTIQUE

"ENmarking"=function(dataset,var.int,proba=0.05){                           #dataset:jeu de données, num.var:variables quali d'intérêt SANS la variable à marquer, var.int:variable quali à marquer, ntripl:nb de triplets que l'utilisateur souhaite garder

res=vector("list",nlevels(dataset[,var.int]))

#isolement de la variable d'intérêt
var.interet=dataset[,var.int]
dataset2=dataset[,-var.int]

#selection des numéros de variables quali et quanti
quali=c()
for (i in 1:ncol(dataset2)){
if (is.factor(dataset2[,i])==TRUE){
quali=c(quali,i)}}

total=1:ncol(dataset2)
quanti=setdiff(total,quali)

dataset2=cbind(dataset2[,quali],var.interet)
colnames(dataset2)[ncol(dataset2)]= colnames(dataset)[var.int]


for (m in 1:nlevels(dataset[,var.int])){                                      #boucle pour faire toutes les modalités de la variable d'intérêt

#PREMIER NIVEAU : VARIABLES
#Test du chi² entre la variable de classe et toutes les variables du jeu de données
chi.res=chisq.desc(dataset2,ncol(dataset2),c(1:(ncol(dataset2)-1)),print=FALSE) #tableau des distances du chi² entre la variable d'intérêt et les autres variables


#Détermination des variables dont au moins une modalité est significative pour chaque classe

#généralisation
keep_nom=matrix(0,nlevels(dataset2[,ncol(dataset2)]),(ncol(dataset2)-1))
keep_pc=matrix(0,nlevels(dataset2[,ncol(dataset2)]),(ncol(dataset2)-1))

for (k in 1:length(chi.res)){
	for (i in 1:nlevels(dataset2[,ncol(dataset2)])){
			if (min(chi.res[[k]][[3]][i,])<proba)
			keep_nom[i,k]=names(dataset2)[k]
			keep_pc[i,k]=min(chi.res[[k]][[3]][i,])
	}
	}

rownames(keep_nom)=rownames(keep_pc)=c(levels(dataset2[,ncol(dataset2)]))                                                 #obtention d'un tableau avec autant de lignes que de modalités de la variable à marquer avec les noms des variables pour lesquelles au moins une modalité est significative et des 0 sinon. + une ligne avec les tests du chi²


  class=rbind(keep_nom[m,],keep_pc[m,])
  delete=which(class[1,]==0)
if(length(delete)==0)
	class2=class
if(length(delete)!=0)
	class2=class[,-delete]

if (length(class2)==2){
	res1=t(as.matrix(class2))
}
if (length(class2)!=2){
	sortclass=class2[,order(as.numeric(as.character(class2[2,])))]       #tri par probabilités croissantes
	no.delete=which(as.numeric(as.character(sortclass[2,]))<=proba)
	resclass=sortclass[,no.delete]
	if (length(class2)<2){
    if (length(class2)==0){
	res1=paste("No significant variable for the category")
	res2=paste("No significant pair for the category")
	res3=paste("No significant triplet for the category")
  catdes_lev1=catdes_lev2=catdes_lev3=NULL
  res1=matrix(0,0,0)}
	  if (length(class2)==1){
	res1=paste("Only one significant variable for the category : ",class2)
	res2=paste("No significant pair for the category")
	res3=paste("No significant triplet for the category")}
  }
	else{
	res1=resclass[,1:min(10,ncol(resclass))]                             #premier niveau : on ne garde qu'un maximum de 10 variables caractérisantes
	res1=t(res1)}
}
if (nrow(res1)==1){                                                    #si une seule variable significative, pas de doublets ni triplets possibles
	res2=paste("No significant pair for the category")
	res3=paste("No significant triplet for the category")
	cat=cbind.data.frame(var.interet,dataset2[,res1[1,1]])
	colnames(cat)=c(colnames(dataset2)[ncol(dataset2)],paste(res1[1,1]))
	cat2=catdes2(cat,1,proba)
	res4=cat2[[2]][[m]]
	catdes_lev1=res4
	catdes_lev2=catdes_lev3=NULL

}
if (nrow(res1)>1){

	cat=cbind.data.frame(var.interet,dataset2[,res1[,1]])
	colnames(cat)=c(colnames(dataset2)[ncol(dataset2)],paste(res1[,1]))
	cat2=catdes2(cat,1,proba)
	res4=cat2[[2]][[m]]
	catdes_lev1=res4
	catdes_lev2=catdes_lev3=NULL
#DEUXIEME NIVEAU : DOUBLETS

tab=dataset2[,res1[,1]]                                                #on ne garde que les colonnes du jeu de données qui correspondent aux variables obtenues à l'étape précédente

pair=matrix(0,nrow(dataset2),1)                                        #création des variables doubles à partir des variables de res1
for (i in 1:nrow(dataset2)){
	pair[i]=paste(c(paste(tab[i,1]),paste(tab[i,2])),collapse="_")
}
colnames(pair)=paste(c(paste(colnames(tab)[1]),paste(colnames(tab)[2])),collapse="_")

if(nrow(res1)==2){                                                    #si seulement deux variables simples significatives, alors une seule variable doublet possible et aucun triplet
varint=as.matrix(dataset2[,ncol(dataset2)])
colnames(varint)=names(dataset2)[ncol(dataset2)]
concaten=as.data.frame(cbind(pair,varint))
chi.res=chisq.desc(concaten,ncol(concaten),c(1:(ncol(concaten)-1)),print=FALSE) #chi² entre la variable double et la variable d'intérêt

keep_nom=matrix(0,nlevels(dataset2[,ncol(dataset2)]),1)
keep_pc=matrix(0,nlevels(dataset2[,ncol(dataset2)]),1)
for (i in 1:nlevels(dataset2[,ncol(dataset2)])){
		if (min(chi.res[[1]][[3]][i,])<proba)                            #suppression des probabilités inférieures à 5%
		keep_nom[i,1]=colnames(pair)[1]
		keep_pc[i,1]=min(chi.res[[1]][[3]][i,])

}

rownames(keep_nom)=rownames(keep_pc)=c(levels(dataset2[,ncol(dataset2)]))

tabb=rbind(keep_nom[m,],keep_pc[m,])

doubl=as.matrix(tabb)
delete=which(doubl[2,]>=proba & doubl[2,]==0)                                #sélection des doublets significatifs
if(length(delete)!=0)
	doubl2=doubl[-delete,]
if(length(delete)==0)
	doubl2=doubl
res2=doubl2	                              					    #sélection d'un maximum de 10 doublets significatifs
res2=t(res2)
rownames(res2)=c()
res3=paste("No significant triplet for the category")

cat=cbind.data.frame(varint,dataset2[,res1[1,1]]:dataset2[,res1[2,1]])
colnames(cat)=c(colnames(varint),paste(c(res1[1,1],res1[2,1]),collapse=":"))
cat2=catdes2(cat,1,proba)
res4=cat2[[2]][[m]]
catdes_lev2=res4
catdes_lev3=NULL
}


if(nrow(res1)!=2){                                                          #si plusieurs variables simples significatives, création de tous les doublets possibles
for (i in 1:ncol(tab)){
	for (j in 3:ncol(tab)){
		pair2=matrix(0,nrow(dataset2),1)
		if (j>i){
			for (k in 1:nrow(tab)){
				pair2[k]=paste(c(paste(tab[k,i]),paste(tab[k,j])),collapse="_")
			}
			colnames(pair2)=paste(c(paste(colnames(tab)[i]),paste(colnames(tab)[j])),collapse="_")
			pair=cbind(pair, pair2)
		}
	}
}
varint=as.matrix(dataset2[,ncol(dataset2)])
colnames(varint)=names(dataset2)[ncol(dataset2)]
concaten=as.data.frame(cbind(pair,varint))
chi.res=chisq.desc(concaten,ncol(concaten),c(1:(ncol(concaten)-1)),print=FALSE) #chi² entre les doublets et la variable d'intérêt

#initialisation
keep_nom=matrix(0,nlevels(dataset2[,ncol(dataset2)]),ncol(pair))
keep_pc=matrix(0,nlevels(dataset2[,ncol(dataset2)]),ncol(pair))

#généralisation
for (k in 1:length(chi.res)){

	for (i in 1:nlevels(dataset2[,ncol(dataset2)])){
			if (min(chi.res[[k]][[3]][i,])<proba)
			keep_nom[i,k]=colnames(pair)[k]
      keep_pc[i,k]=min(chi.res[[k]][[3]][i,])
	}

rownames(keep_nom)=rownames(keep_pc)=c(levels(dataset2[,ncol(dataset2)]))                                             #obtention d'un tableau avec autant de lignes que de modalités de la variable à marquer avec les noms des variables pour lesquelles au moins une modalité est significative et des 0 sinon. + une ligne avec les tests du chi²
}

#Ensemble des doublets gardés pour une classe
class=rbind(keep_nom[m,],keep_pc[m,])                                     #récupération des lignes correspondant à la modalité choisie et aux probabilités
delete=which(class[1,]==0)                                           #suppression des 0
if(length(delete)==0)
	class2=class
if(length(delete)!=0)
	class2=class[,-delete]
sortclass=class2[,order(as.numeric(as.character(class2[2,])))]       #tri par probabilités croissantes
no.delete=which(as.numeric(as.character(sortclass[2,]))<=proba)       #suppression des doublets non significatifs
resclass=sortclass[,c(no.delete)]                                    #deuxième niveau

if(is.null(dim(resclass))==TRUE){
	if(length(no.delete)==0){                                          #si resclass est vide, res2=NULL
		res2=NULL}
	if(length(no.delete)!=0){						                               #si resclass n'est pas vide, existence de res2=maximum 10 doublets significatifs
		res2=resclass[,1:min(10,dim(resclass)[2])]
		res2=t(res2)}}
		
if(is.null(dim(resclass))==FALSE){
	if(dim(resclass)[1]==0){						                               #si 0 lignes, res2=NULL
		res2=NULL}
	if(dim(resclass)[1]!=0){
    	res2=resclass[,1:min(10,dim(resclass)[2])]                       #sinon, sélection d'un max de 10 doublets significatifs
	res2=t(res2)}}
	
if(is.null(res2)==TRUE){           									                 #si res2 n'existe pas, res3 non plus
	res2=paste("No significant pair for the category")
	res3=paste("No significant triplet for the category")
	catdes_lev2=catdes_lev3=NULL
}
if(is.null(res2)==FALSE){                                            #si res2 n'est pas nul, calcul de res3

	cat=cbind.data.frame(varint,pair[,res2[,1]])
  cat2=catdes2(cat,1,proba)
	res4=cat2[[2]][[m]]
	catdes_lev2=res4
	catdes_lev3=NULL

#TROISIEME NIVEAU : TRIPLETS

triplets=matrix(,nrow(dataset2),0)                                   #création des variables triples (à partir des variables de res1)

for (i in 1:ncol(tab)){
	for (j in 2:ncol(tab)){
		for (k in 3:ncol(tab)){
			triplets2=matrix(0,nrow(tab),1)
			if (k>j & j>i){
			for (l in 1:nrow(dataset2)){
				triplets2[l]=paste(c(paste(tab[l,i]),paste(tab[l,j]),paste(tab[l,k])),collapse="_")
			}
			colnames(triplets2)=paste(c(paste(colnames(tab)[i]),paste(colnames(tab)[j]),paste(colnames(tab)[k])),collapse="_")
			triplets=cbind(triplets, triplets2)
			}
		}
	}
}

varint=as.matrix(dataset2[,ncol(dataset2)])
colnames(varint)=names(dataset2)[ncol(dataset2)]
concaten=as.data.frame(cbind(triplets,varint))
chi.res=chisq.desc(concaten,ncol(concaten),c(1:(ncol(concaten)-1)),print=FALSE) #chi² entre les variables triples et la variable d'intérêt


#généralisation
keep_nom=matrix(0,nlevels(dataset2[,ncol(dataset2)]),length(chi.res))
keep_pc=matrix(0,nlevels(dataset2[,ncol(dataset2)]),length(chi.res))

for (k in 1:length(chi.res)){
	for (i in 1:nlevels(dataset2[,ncol(dataset2)])){
			if (min(chi.res[[k]][[3]][i,])<proba)
			keep_nom[i,k]=colnames(triplets)[k]
			keep_pc[i,k]=min(chi.res[[k]][[3]][i,])
	}
	}

rownames(keep_nom)=rownames(keep_pc)=c(levels(dataset2[,ncol(dataset2)]))                                                 #obtention d'un tableau avec autant de lignes que de modalités de la variable à marquer avec les noms des variables pour lesquelles au moins une modalité est significative et des 0 sinon. + une ligne avec les tests du chi²


  class=rbind(keep_nom[m,],keep_pc[m,])
  delete=which(class[1,]==0)
  if(length(delete)==0)
	 class2=class
  if(length(delete)!=0)
	 class2=class[,-delete]
sortclass=as.matrix(class2[,order(as.numeric(as.character(class2[2,])))])       #tri par probabilités croissantes
no.delete=which(as.numeric(as.character(sortclass[2,]))<=proba)       #suppression des doublets non significatifs
resclass=as.matrix(sortclass[,c(no.delete)])

if(is.null(dim(resclass))==TRUE){
	if(length(no.delete)==0){                                          #si resclass est vide, res2=NULL
		res3=NULL}
	if(length(no.delete)!=0){						                               #si resclass n'est pas vide, existence de res2=maximum 10 doublets significatifs
		res3=resclass[,1:min(10,dim(resclass)[2])]
		res3=t(res3)
    rownames(res3)=c()}}
		
if(is.null(dim(resclass))==FALSE){
	if(dim(resclass)[1]==0){						                               #si 0 lignes, res2=NULL
		res3=NULL}
	if(dim(resclass)[1]!=0){
    	res3=resclass[,1:min(10,dim(resclass)[2])]                       #sinon, sélection d'un max de 10 doublets significatifs
	res3=t(res3)}}
	
if(is.null(res3)==TRUE){           									                 #si res2 n'existe pas, res3 non plus
	res3=paste("No significant triplet for the category")}
	
if (is.null(res3)==FALSE){
  cat=matrix(,nrow(dataset2),nrow(res3)+1)
	cat=cbind.data.frame(varint,triplets[,res3[,1]])
	if(dim(triplets)[2]==1)
  colnames(cat)=c(colnames(varint),colnames(triplets))
  cat2=catdes2(cat,1,proba)
	res4=cat2[[2]][[m]]
	catdes_lev3=res4}

}}}

#si aucune variable liée
if (nrow(res1)<1)
res1=paste("No significant variable for the category")

res[[m]]=list()
if (length(res1)==1)                                                            #remplissage de res par les trois niveaux de marquage
res[[m]]$lev_1$marking=res1
if (length(res1)>1){
res[[m]]$lev_1$marking=as.matrix(res1[,1])
colnames(res[[m]]$lev_1$marking)=c("Variable(s)")}
res[[m]]$lev_1$catdes=catdes_lev1
if (length(res2)==1)                                                            #remplissage de res par les trois niveaux de marquage
res[[m]]$lev_2$marking=res2
if (length(res2)>1){
res[[m]]$lev_2$marking=as.matrix(res2[,1])
colnames(res[[m]]$lev_2$marking)=c("Pair(s)")}
res[[m]]$lev_2$catdes=catdes_lev2
if (length(res3)==1)                                                            #remplissage de res par les trois niveaux de marquage
res[[m]]$lev_3$marking=res3
if (length(res3)>1){
res[[m]]$lev_3$marking=as.matrix(res3[,1])
colnames(res[[m]]$lev_3$marking)=c("Triplet(s)")}
res[[m]]$lev_3$catdes=catdes_lev3
}

names(res)=levels(dataset[,var.int])
return(res)
}
