#ENMCA
"ENMCA"=function (dataset,ind.sup = NULL, quanti.sup = NULL, quali.sup = NULL, axes = c(1, 2), row.w = NULL,level.ventil=0,signif=FALSE,proba=0.05,report=FALSE,language="english"){              #num.var=numéros des variables quali, var.class : variable d'intérêt pour catdes
if(40-length(dev.list())<16){
print("There are too many devices open. Please close them and launch ENMCA again")
}else{
a=getwd()
if(report==TRUE){
    dir.create(paste(a,"/EnQuireR/",sep=""))
}
fn=function (dataset,ind.sup,quanti.sup,quali.sup,proba,signif,report){

acm=MCA2(dataset,ind.sup=ind.sup,quali.sup=quali.sup,quanti.sup=quanti.sup,level.ventil=level.ventil,graph=FALSE)                                                #ACM
#ACM sans les éléments supplémentaires
dataset.min=dataset

if (!is.null(ind.sup)){
dataset.min=dataset.min[-ind.sup,] }      #si individus supplémentaires

if (!is.null(quali.sup)){
if (!is.null(quanti.sup)){
 dataset.min=dataset.min[,-c(quali.sup,quanti.sup)]}  #si quali.sup et quanti.sup
else {
 dataset.min=dataset.min[,-quali.sup]}}     #si quali.sup uniquement
 else{
 if (!is.null(quanti.sup)){
 dataset.min=dataset.min[,-quanti.sup]}      #si quanti sup uniquement
}       


acm_ss_clust=MCA2(dataset.min,ind.sup=NULL,quali.sup=NULL,quanti.sup=NULL,graph=FALSE)                                                #ACM
if (signif==TRUE){
res.signif=p_inertia(dataset.min)} #ré-echantillonnage

tab=data.frame(acm$ind$coord)                                                     #récupération des coordonnées des individus
dist=dist(tab, method = "euclidean")                                              #matrice des distances
hclust=hclust(dist, method = "ward")                          #classification
dev.new()
tab$hclus.label=assignCluster(tab,tab,cuttree.graph(hclust))
group=tab$hclus.label                       #création de la variable de classe
if(report==TRUE){
    setwd(paste(a,"/EnQuireR",sep=""))
    savePlot(filename = "Clustplot",type = "pdf",device = dev.cur(),restoreConsole = TRUE)
    setwd(a)
}
if (!is.null(ind.sup)){             #ajout de la variable de classe au jeu de données
tab2=cbind(acm$call$X[-ind.sup,],group)}
else {
tab2=cbind(acm$call$X,group)}
                                              
res1=catdes(tab2,num.var=ncol(tab2),proba=proba)                                   #description des classes
if (!is.null(ind.sup)){ 
acm$call$X=cbind(acm$call$X[-ind.sup,],tab$hclus.label)
plot.MCA(acm,choix="ind",col.ind=as.numeric(acm$call$X[,ncol(acm$call$X)]),invisible=c("var","ind.sup","quali.sup"),col.ind.sup="none",cex=0.7)
legend("topleft", legend = levels(as.factor(acm$call$X[,ncol(acm$call$X)])), text.col = 1:nlevels(as.factor(acm$call$X[,ncol(acm$call$X)])), cex = 0.8)}

else{
acm$call$X=cbind(acm$call$X,tab$hclus.label)
plot.MCA(acm,choix="ind",habillage=ncol(acm$call$X),invisible=c("var","quali.sup"),cex=0.7)}    #plot avec habillage par classe. variables rendues invisibles                                          #ajout de la variable de classe au jeu de données
                   
if(report==TRUE){
    setwd(paste(a,"/EnQuireR",sep=""))
    savePlot(filename = "MCAplotindclasse",type = "pdf",device = dev.cur(),restoreConsole = TRUE)
    setwd(a)
}
plot.MCA(acm,choix="ind",habillage="none",invisible=c("ind","ind.sup"),cex=0.7)                   #plot avec habillage par classe. variables rendues invisibles
if(report==TRUE){
    setwd(paste(a,"/EnQuireR",sep=""))
    savePlot(filename = "MCAplotvar",type = "pdf",device = dev.cur(),restoreConsole = TRUE)
    setwd(a)
}
plot.MCA(acm,choix="ind",habillage="none",invisible=c("var","quali.sup"),cex=0.7,label="none")                   #plot variables rendues invisibles
if(report==TRUE){
    setwd(paste(a,"/EnQuireR",sep=""))
    savePlot(filename = "MCAplotind",type = "pdf",device = dev.cur(),restoreConsole = TRUE)
    setwd(a)
}

if (!is.null(quanti.sup)){
plot.MCA(acm,choix="quanti.sup")                   #plot variables quanti sup si elles existent
if(report==TRUE){
    setwd(paste(a,"/EnQuireR",sep=""))
    savePlot(filename = "MCAplotquantisup",type = "pdf",device = dev.cur(),restoreConsole = TRUE)
    setwd(a)
}}

matP=acm$ind$coord[,c(1,2)]
eig=acm$eig #valeurs propres de l'acm
v=tab2[,ncol(tab2)] #variable de classe
dev.new()
ENellipse(matP,cluster=v,eig=eig, resampling="cluster")
if(report==TRUE){
    setwd(paste(a,"/EnQuireR",sep=""))
    savePlot(filename = "Ellipse",type = "pdf",device = dev.cur(),restoreConsole = TRUE)
    setwd(a)
}


ENlisib(acm,0.25,0.25,axes=c(1,2))
if(report==TRUE){
    setwd(paste(a,"/EnQuireR",sep=""))
    savePlot(filename = "Lisib",type = "pdf",device = dev.cur(),restoreConsole = TRUE)
    setwd(a)
}

res=list()
res$MCA=acm
res$catdes=res1
res$data=tab2
res.acm=acm_ss_clust
if (signif==TRUE){
res$signif=res.signif}

dev.new()
ENdensity(res$MCA,res$data[,ncol(res$data)]) #trace les courbes de densité selon la variable de classe
if(report==TRUE){
    setwd(paste(a,"/EnQuireR",sep=""))
    savePlot(filename = "Density",type = "pdf",device = dev.cur(),restoreConsole = TRUE)
    setwd(a)
}

dev.new()
ENdensity(res$MCA) #trace les corbes de densité pour l'ensemble des individus
if(report==TRUE){
    setwd(paste(a,"/EnQuireR",sep=""))
    savePlot(filename = "Densitypop",type = "pdf",device = dev.cur(),restoreConsole = TRUE)
    setwd(a)
}

dev.new()
barplotgroupes(res$data[,ncol(res$data)]) #barplot pour la variables de classe
if(report==TRUE){
    setwd(paste(a,"/EnQuireR",sep=""))
    savePlot(filename = "barplotgr",type = "pdf",device = dev.cur(),restoreConsole = TRUE)
    setwd(a)
}

if(length(rownames(res$catdes$test.chi))>10){ #barplots pour les 10 variables les plus liées à la variable de classe
for(i in 1:10){
    dev.new()
    barplotclustervsvar(res$data[,rownames(res$catdes$test.chi)[i]],res$data[,ncol(res$data)],res$data)
    if(report==TRUE){
        setwd(paste(a,"/EnQuireR",sep=""))
        savePlot(filename = paste("barvarclus",i,sep=""),type = "pdf",device = dev.cur(),restoreConsole = TRUE)
        setwd(a)
    }
}
}else{
for(i in 1:length(rownames(res$catdes$test.chi))){
    dev.new()
    barplotclustervsvar(res$data[,rownames(res$catdes$test.chi)[i]],res$data[,ncol(res$data)],res$data)
    if(report==TRUE){
        setwd(paste(a,"/EnQuireR",sep=""))
        savePlot(filename = paste("barvarclus",i,sep=""),type = "pdf",device = dev.cur(),restoreConsole = TRUE)
        setwd(a)
    }
}
}

aa=res$catdes
bb=res.acm
cc=dataset
dd=res$data
num.var=ncol(dataset)-length(quali.sup)-length(quanti.sup)

assign("hclus",tab2,envir=.GlobalEnv)
assign("num.var",num.var,envir=.GlobalEnv)
assign("catdes",aa,envir=.GlobalEnv)
assign("res.acm",bb,envir=.GlobalEnv)
assign("dataset",cc,envir=.GlobalEnv)
assign("dataset2",dd,envir=.GlobalEnv)
#ACM avec les variables supplémentaires, la variable de groupes et les individus supplémentaires
assign("res.acm.tot",res$MCA,envir=.GlobalEnv)
#retourne le résultat du réechantillonnage sur les pourcentages d'inertie si demandé
if (signif==TRUE){
assign("res.signif",res.signif$resum,envir=.GlobalEnv)}
return(res)
#res=list()
#res$MCA=acm
#res.acm=acm
#bb=res.acm
#cc=dataset
#assign("res.acm",bb,envir=.GlobalEnv)
#assign("dataset",cc,envir=.GlobalEnv)
#return(res)
#return(dataset)
}
env=environment(fn)

Sweave_f=function(){
if(report==TRUE){
nomsfichiers=list.files(path = paste(.libPaths()[1],"/EnQuireR/Sweave/sty/",sep=""), pattern = NULL, all.files = FALSE,full.names = FALSE, recursive = FALSE,ignore.case = FALSE)
for(i in 1:length(nomsfichiers)){
file.copy(from=paste(.libPaths()[1],"/EnQuireR/Sweave/sty/",nomsfichiers[i],sep=""),to=paste(a,"/EnQuireR/",nomsfichiers[i],sep=""), overwrite = TRUE)
}
#file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/subfigure.sty",sep=""),paste(a,"/MCA/subfigure.sty",sep=""))
#file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/fancyvrb.sty",sep=""),paste(a,"/MCA/fancyvrb.sty",sep=""))
#file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/Sweave.sty",sep=""),paste(a,"/MCA/Sweave.sty",sep=""))
#file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/upquote.sty",sep=""),paste(a,"/MCA/upquote.sty",sep=""))
#file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/algorithmic.sty",sep=""),paste(a,"/MCA/algorithmic.sty",sep=""))
setwd(paste(a,"/EnQuireR",sep=""))
if (language=="english"){
Sweave(paste(.libPaths()[1],"/EnQuireR/Sweave/MCA/En/Multivariate_report.Rnw",sep=""), driver = RweaveLatex(),syntax = getOption("SweaveSyntax"))
Sweave(paste(.libPaths()[1],"/EnQuireR/Sweave/MCA/En/Multivariate_presentation.Rnw",sep=""), driver = RweaveLatex(),syntax = getOption("SweaveSyntax"))
tools::texi2dvi(paste(a,"/EnQuireR/Multivariate_report.tex",sep=""), pdf=TRUE)
tools::texi2dvi(paste(a,"/EnQuireR/Multivariate_presentation.tex",sep=""), pdf=TRUE)
}
if (language=="french"){
Sweave(paste(.libPaths()[1],"/EnQuireR/Sweave/MCA/Fr/Multivariate_report.Rnw",sep=""), driver = RweaveLatex(),syntax = getOption("SweaveSyntax"))
}
if (language=="german"){
Sweave(paste(.libPaths()[1],"/EnQuireR/Sweave/MCA/De/Multivariate_report.Rnw",sep=""), driver = RweaveLatex(),syntax = getOption("SweaveSyntax"))
}
tools::texi2dvi(paste(a,"/EnQuireR/Multivariate_report.tex",sep=""), pdf=TRUE)
setwd(a)
}
}
res=fn(dataset,ind.sup=ind.sup,quali.sup=quali.sup,quanti.sup=quanti.sup,proba=proba,signif=signif,report=report)
Sweave_f()
return(res)
}
}
