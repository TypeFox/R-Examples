fahst=function(don,group,alpha=0.05,graph=TRUE,axes=c(1,2),name.group=NULL,ncp=5,B=200,ncp.boot=2){

don=as.data.frame(don)
I=nrow(don)
J=length(group)


for (i in 1:ncol(don)){
don[,i]=as.factor(don[,i])}

if (is.null(name.group)) 
name.group <- paste("Sj", 1:length(group), sep = ".")

########################################Graphiques préliminaires
if (graph){

#Nombre de niveaux par sujet
group2=as.factor(group)
plot(group2,main="Number of levels per subject")

#Nombre de groupes au niveau 1 de chaque sujet
niv1=cumsum(group)-group+1
lev1=rep(NA,J)
for (i in 1:J){
lev1[i]=length(levels(don[,niv1[i]]))}
lev1b=as.factor(lev1)
dev.new()  
plot(lev1b,main="Number of groups formed from first levels")

#Nombre de groupes au dernier niveau de chaque sujet
nivd=cumsum(group)
levd=rep(NA,J)
for (i in 1:J){
levd[i]=length(levels(don[,nivd[i]]))}
levdb=as.factor(levd)
dev.new()
plot(levdb,main="Number of groups formed from last levels")

#Nombre d'objets par groupe au niveau 1 de chaque sujet
nbp1=strsplit(summary(don[,niv1],maxsum=max(lev1)),":")
agg1=rep(0,J*max(lev1))
for (i in 1:(J*max(lev1))){
agg1[i]=nbp1[[i]][2]}
agg1b=na.omit(agg1)
agg1b=as.factor(agg1b)
dev.new()
plot(agg1b,main="Number of objects per group for the first levels")

#Nombre d'objets par groupe au dernier niveau de chaque sujet
nbpd=strsplit(summary(don[,nivd],maxsum=max(levd)),":")
aggd=rep(0,J*max(levd))
for (i in 1:(J*max(levd))){
aggd[i]=nbpd[[i]][2]}
aggdb=na.omit(aggd)
aggdb=as.factor(aggdb)
dev.new()
plot(aggdb,main="Number of objects per group for the last levels")

}
########################################Fin graphiques préliminaires

#AFM
afm=MFA(don,group=group,type=rep("n",J),name.group=name.group,graph=F,ncp=ncp)


##################rapport de corrélation##################################
eta2 <- function(x, gpe) {
vartot <- function(x) {
res <- sum((x - mean(x))^2)
return(res)}
varinter <- function(x, gpe) {
moyennes <- tapply(x, gpe, mean)
effectifs <- tapply(x, gpe, length)
res <- (sum(effectifs * (moyennes - mean(x))^2))
return(res)}
res <- varinter(x, gpe)/vartot(x)
return(res)}
###########################################################################

#calcul pour chaque dimension du rapport de corrélation
coord.niv=matrix(0,ncol(don),ncp)
rownames(coord.niv)=colnames(don)
colnames(coord.niv)=colnames(afm$ind$coord)
for (i in 1:ncol(don)){
coord.niv[i,]=apply(afm$ind$coord,2,eta2,don[,i])}

########################################Graphiques de l'AFM
if (graph){
#Graph des individus
plot.MFA(afm,choix="ind",invisible="quali",axes=axes,new.plot=TRUE)
#Graph des mots
plot.MFA(afm,choix="ind",invisible="ind",axes=axes,new.plot=TRUE)
#Graph des individus et des mots
plot.MFA(afm,choix="ind",axes=axes,new.plot=TRUE)
#Graph des groupes au niveau des sujets
plot.MFA(afm,choix="group",axes=axes,new.plot=TRUE)
#Graph des niveaux
dev.new(width = 8, height = 8)
plot(coord.niv[,axes], xlab = paste("Dim ", axes[1], " (", signif(afm$eig[axes[1],2], 4), "%)", sep = ""), ylab = paste("Dim ", axes[2], " (", signif(afm$eig[axes[2],2], 4), "%)", sep = ""), xlim = c(0,1), ylim = c(0, 1), pch = 17, main = "Levels representation")
text(coord.niv[, axes[1]], y = coord.niv[, axes[2]], labels = rownames(coord.niv),pos = 3)
#Graph des niveaux et trajectoires
dev.new(width = 8, height = 8)
plot(coord.niv[,axes], xlab = paste("Dim ", axes[1], " (", signif(afm$eig[axes[1],2], 4), "%)", sep = ""), ylab = paste("Dim ", axes[2], " (", signif(afm$eig[axes[2],2], 4), "%)", sep = ""), xlim = c(0,1), ylim = c(0, 1), pch = 17, main = "Levels representation and trajectories")
text(coord.niv[, axes[1]], y = coord.niv[, axes[2]], labels = rownames(coord.niv),pos = 3)
subj=0 
for (j in 1:length(group)){
if (group[j]!=1){
for (i in 1:(group[j]-1)){
lines(x=coord.niv[(subj+i):(subj+i+1),axes[1]],y=coord.niv[(subj+i):(subj+i+1),axes[2]])}}
subj=subj+group[j]
}
}
########################################Fin graphiques de l'AFM

########################################Ellipses
if (graph){
boot (don,method="hsort",group=group,ncp=ncp.boot,nbsim=B, level.conf= 1-alpha)}
########################################Ellipses

afm_call = list(X = afm$call$X, col.w = afm$call$col.w, 
        row.w = afm$call$row.w, ncp = afm$call$ncp, 
        name.group = afm$call$name.group,group=group,mfa=afm)
var_afm = list(coord=afm$quali.var$coord,contrib=afm$quali.var$contrib,cos2=afm$quali.var$cos2,v.test=afm$quali.var$v.test,coord.lev=coord.niv)

res = list(eig = afm$eig, var = var_afm, ind = afm$ind, group = afm$group,
        call = afm_call)
class(res) <- c("fahst", "list ")

return(res)}