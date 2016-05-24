fasnt=function(don,first="nappe",B=100,axes=c(1,2),alpha=0.05,ncp=5,graph=TRUE,name.group=NULL,sep.word=" ",word.min=5,ncp.boot=2){

don=data.frame(don)
I=nrow(don)
J=ncol(don)/3

#Si pour chaque sujet d'abord la catégorisation  puis la nappe
don_input=don
if (first=="catego"){
don2=don
  for (i in 1:J){
  don2[,(3*(i-1)+1):(3*(i-1)+2)]=don[,(3*(i-1)+2):(3*i)]
  don2[,(3*i)]=don[,(3*(i-1)+1)]
  colnames(don2)[(3*(i-1)+1):(3*(i-1)+2)]=colnames(don)[(3*(i-1)+2):(3*i)]
  colnames(don2)[(3*i)]=colnames(don)[(3*(i-1)+1)]
  }
  don=don2
}

if(is.null(name.group)){
name.group=list(rep("aa",2*J),c(paste("J",1:J,sep="")))
name.group[[1]][seq(from=1,to=2*J,by=2)]=paste("J",1:J,"n",sep="")
name.group[[1]][seq(from=2,to=2*J,by=2)]=paste("J",1:J,"c",sep="") }

##############################FONCTIONS NECESSAIRES#############################

plot.HMFA.group.niv1 <- function(res.hmfa, coord = c(1, 2),title = NULL, cex = 1, sub.title = NULL) {
    lab.x <- paste("Dim ", coord[1], " (", signif(res.hmfa$eig[coord[1],
        2], 4), " %)", sep = "")
    lab.y <- paste("Dim ", coord[2], " (", signif(res.hmfa$eig[coord[2],
        2], 4), " %)", sep = "")

        if (is.null(title))
            title <- "Groups representation"
        else sub.title <- "Groups representation"

            coord.actif <- res.hmfa$group$coord[[1]][, coord]
             gp=rep(c("1","2"),(length(coord.actif[,1])/2))
                plot(coord.actif, xlab = lab.x, ylab = lab.y,
                  xlim = c(0, 1), ylim = c(0, 1), pch = 17, col = color[2],
                  cex = cex, main = title, cex.main = cex, asp = 1)
                points(coord.actif[, 1], coord.actif[, 2], col = gp,
                pch = 17, cex = cex)
                text(coord.actif[, 1], y = coord.actif[, 2],
                  labels = rownames(coord.actif), pos = 3, col = gp)

        title(sub = sub.title, cex.sub = cex, font.sub = 2, col.sub = "steelblue4",
            adj = 0, line = 3.8)
    }

 color = c("black", "red", "green3", "blue", "cyan", "magenta",
            "darkgray", "darkgoldenrod", "darkgreen", "violet",
            "turquoise", "orange", "lightpink", "lavender", "yellow",
            "lightgreen", "lightgrey", "lightblue", "darkkhaki",
            "darkmagenta", "darkolivegreen", "lightcyan", "darkorange",
            "darkorchid", "darkred", "darksalmon", "darkseagreen",
            "darkslateblue", "darkslategray", "darkslategrey",
            "darkturquoise", "darkviolet", "lightgray", "lightsalmon",
            "lightyellow", "maroon")

#####################################DEBUT ANALYSE##############################


#AFMH
hierar=list(rep(c(2,1),J),rep(2,J))
afmh=HMFA(don,H=hierar,type=rep(c("c","n"),J),name.group=name.group,graph=F,ncp=ncp)
if (graph){
lab.x <- paste("Dim ", axes[1], " (", signif(afmh$eig[axes[1],2], 4), " %)", sep = "")
lab.y <- paste("Dim ", axes[2], " (", signif(afmh$eig[axes[2],2], 4), " %)", sep = "")

#graph des individus
xmin <- min(afmh$ind$coord[, axes[1]])
xmax <- max(afmh$ind$coord[, axes[1]])
ymin <- min(afmh$ind$coord[, axes[2]])
ymax <- max(afmh$ind$coord[, axes[2]])
x <- c(xmin, xmax) * 1.1
y <- c(ymin, ymax) * 1.1
dev.new()
plot(0, 0, xlab = lab.x, ylab = lab.y, xlim = x,ylim = y, col = "white",main="Individuals factor map", asp = 1)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
points(afmh$ind$coord[, axes], pch = 15)
text(afmh$ind$coord[, axes[1]],afmh$ind$coord[, axes[2]], rownames(afmh$ind$coord),pos = 3, offset = 0.2, cex = 0.8)

#graph des individus et des variables quali
xmin <- min(afmh$ind$coord[, axes[1]],afmh$quali.var$coord[, axes[1]])
xmax <- max(afmh$ind$coord[, axes[1]],afmh$quali.var$coord[, axes[1]])
ymin <- min(afmh$ind$coord[, axes[2]],afmh$quali.var$coord[, axes[2]])
ymax <- max(afmh$ind$coord[, axes[2]],afmh$quali.var$coord[, axes[2]])
x <- c(xmin, xmax) * 1.1
y <- c(ymin, ymax) * 1.1
dev.new()
plot(0, 0, xlab = lab.x, ylab = lab.y, xlim = x,ylim = y, col = "white",main="Individuals factor map", asp = 1)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
points(afmh$ind$coord[, axes], pch = 15,col=2)
points(afmh$quali.var$coord[, axes], pch = 15)
text(afmh$ind$coord[, axes[1]],afmh$ind$coord[, axes[2]],col=2, rownames(afmh$ind$coord),pos = 3, offset = 0.2, cex = 0.8)
text(afmh$quali.var$coord[, axes[1]],afmh$quali.var$coord[, axes[2]], rownames(afmh$quali.var$coord),pos = 3, offset = 0.2, cex = 0.8)

#Graph des points partiels
inter <- afmh$partial[[2]][, axes, 1]
for (i in 2:J) inter <- rbind(inter, afmh$partial[[2]][,axes, i])
xmin <- min(afmh$ind$coord[, axes[1]],inter[, axes[1]])
xmax <- max(afmh$ind$coord[, axes[1]],inter[, axes[1]])
ymin <- min(afmh$ind$coord[, axes[2]],inter[, axes[2]])
ymax <- max(afmh$ind$coord[, axes[2]],inter[, axes[2]])
x <- c(xmin, xmax) * 1.1
y <- c(ymin, ymax) * 1.1
dev.new()
plot(0, 0, xlab = lab.x, ylab = lab.y, xlim = x,ylim = y, col = "white",main="Individuals factor map", asp = 1)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
points(afmh$ind$coord[, axes], pch = 15,col=1:I)
text(afmh$ind$coord[, axes[1]],afmh$ind$coord[, axes[2]],col=1:I, rownames(afmh$ind$coord),pos = 3, offset = 0.2, cex = 0.8)

for (j in 1:J) {
points(afmh$partial[[2]][, axes, j], col=rep(1:I,times=J),pch = 20, cex = 0.8)
text(afmh$partial[[2]][, axes, j], col=rep(1:I,times=J),labels=rownames(afmh$group$coord[[2]])[j],pos=3,cex=0.5)
for (i in 1:nrow(afmh$partial[[2]]))
  lines(c(afmh$ind$coord[i,axes[1]], afmh$partial[[2]][i, axes[1],j]), c(afmh$ind$coord[i, axes[2]], afmh$partial[[2]][i,axes[2], j]),col=i)
}

#Graph des variables quanti
dev.new()
plot(0, 0, xlab = lab.x, ylab = lab.y,xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1),main="Correlation circle", col = "white",asp = 1)
x.cercle <- seq(-1, 1, by = 0.01)
y.cercle <- sqrt(1 - x.cercle^2)
lines(x.cercle, y = y.cercle)
lines(x.cercle, y = -y.cercle)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)

coord.var <- afmh$quanti.var$cor[, axes]
for (v in 1:nrow(coord.var)) {
  arrows(0, 0, coord.var[v, 1], coord.var[v,2], length = 0.1, angle = 15, code = 2)
  if (abs(coord.var[v, 1]) > abs(coord.var[v,2])) {
    if (coord.var[v, 1] >= 0)
      pos <- 4
    else pos <- 2
  }
  else {
    if (coord.var[v, 2] >= 0)
      pos <- 3
    else pos <- 1
  }
  text(coord.var[v, 1], y = coord.var[v, 2],labels = rownames(coord.var)[v], pos = pos)
}

#graph des sujets
coord.actif <- afmh$group$coord[[2]][, axes]
dev.new()
plot(coord.actif, xlab = lab.x, ylab = lab.y,xlim = c(0, 1), ylim = c(0, 1), pch = 17, col = color[2], main = "Subjects representation", asp = 1)
points(coord.actif[, 1], coord.actif[, 2], col = color[2],pch = 17)
text(coord.actif[, 1], y = coord.actif[, 2],labels = rownames(coord.actif), pos = 3, col = color[2])

#graph des méthodes
coord.actif <- afmh$group$coord[[1]][, axes]
gp=rep(c("1","2"),(length(coord.actif[,1])/2))
dev.new()
plot(coord.actif, xlab = lab.x, ylab = lab.y,xlim = c(0, 1), ylim = c(0, 1), pch = 17, col = color[2],main = "Method's representation", asp = 1)
points(coord.actif[, 1], coord.actif[, 2], col = gp,pch = 17)
text(coord.actif[, 1], y = coord.actif[, 2],labels = rownames(coord.actif), pos = 3, col = gp)
}


########################################Ellipses

if (graph){
boot (don,method="sortnapping", level.conf = 1-alpha,nbsim=B,ncp=ncp.boot)}



#Indicateurs du napping
#Récupérer que les nappes
nappe_c=data.frame(matrix(NA,I,2*J))
for (i in 1:J){
nappe_c[,(2*(i-1)+1)]=don[,(3*(i-1)+1)]
colnames(nappe_c)[(2*(i-1)+1)]=colnames(don)[(3*(i-1)+1)]
nappe_c[,(2*i)]=don[,(3*(i-1)+2)]
colnames(nappe_c)[(2*i)]=colnames(don)[(3*(i-1)+2)]}

res.nappe=function(don){
I=nrow(don)
J=ncol(don)/2
X=Y=matrix(NA,I,J)
for (i in 1:J){
X[,i]=don[,(2*i-1)]
Y[,i]=don[,(2*i)]}

resultat=matrix(NA,5,3)
rownames(resultat)=c("St. dev. X","St. dev. Y","Range X","Range Y","Dimensionnality")
colnames(resultat)=c("min","median","max")

ecart.type.X=apply(X,2,sd)
ecart.type.Y=apply(Y,2,sd)
min.X=apply(X,2,min)
max.X=apply(X,2,max)
etendue.X=max.X-min.X
min.Y=apply(Y,2,min)
max.Y=apply(Y,2,max)
etendue.Y=max.Y-min.Y

res.pca=matrix(NA,2,J)
dimen=rep(NA,J)
for (i in 1:J){
acp=PCA(don[,(2*(i-1)+1):(2*i)],graph=F,scale.unit=F)
res.pca[,i]=acp$eig[,1]
dimen[i]=1+(res.pca[2,i]^2/res.pca[1,i]^2)}

resultat[1,1]=min(ecart.type.X)
resultat[1,2]=median(ecart.type.X)
resultat[1,3]=max(ecart.type.X)
resultat[2,1]=min(ecart.type.Y)
resultat[2,2]=median(ecart.type.Y)
resultat[2,3]=max(ecart.type.Y)
resultat[3,1]=min(etendue.X)
resultat[3,2]=median(etendue.X)
resultat[3,3]=max(etendue.X)
resultat[4,1]=min(etendue.Y)
resultat[4,2]=median(etendue.Y)
resultat[4,3]=max(etendue.Y)
resultat[5,1]=min(dimen)
resultat[5,2]=median(dimen)
resultat[5,3]=max(dimen)
return(resultat)}

res.nappe2=function(don){
I=nrow(don)
J=ncol(don)/2
X=Y=matrix(NA,I,J)
for (i in 1:J){
X[,i]=don[,(2*i-1)]
Y[,i]=don[,(2*i)]}

resultat=matrix(NA,3,2)
rownames(resultat)=c("PCA X","PCA Y","PCA F1")
colnames(resultat)=c("% inertia dim 1","% inertia plane 1-2")

acp.X=PCA(X,scale.unit=T,graph=F)
acp.Y=PCA(Y,scale.unit=T,graph=F)

res.pca=matrix(NA,I,J)
for (i in 1:J){
acp=PCA(don[,(2*(i-1)+1):(2*i)],graph=F,scale.unit=F)
res.pca[,i]=acp$ind$coord[,1]}

acp.F1=PCA(res.pca,scale.unit=T,graph=F)

resultat[1,1]=acp.X$eig[1,2]
resultat[1,2]=acp.X$eig[2,3]
resultat[2,1]=acp.Y$eig[1,2]
resultat[2,2]=acp.Y$eig[2,3]
resultat[3,1]=acp.F1$eig[1,2]
resultat[3,2]=acp.F1$eig[2,3]
return(list(resultat,res.pca))}

res1_nappe=res.nappe(nappe_c)
res2_nappe=res.nappe2(nappe_c)

#Indicateurs de la catégorisation
catego_c=data.frame(matrix(NA,I,J))
for (i in 1:J){
catego_c[,i]=don[,(3*i)]
colnames(catego_c)[i]=colnames(don)[3*i]}

ordre_prod=order(afmh$ind$coord[,1])
coocc=tab.disjonctif(catego_c)%*%t(tab.disjonctif(catego_c))
coocc_reord=coocc[ordre_prod,ordre_prod]
colnames(coocc_reord)=rownames(coocc_reord)=rownames(afmh$ind$coord)[ordre_prod]

#Nombre de produits par classe
lev=rep(NA,J)
for (i in 1:J){
lev[i]=length(levels(catego_c[,i]))}
nbp=strsplit(summary(catego_c,maxsum=max(lev)),":")
agg=rep(0,J*max(lev))
for (i in 1:(J*max(lev))){
agg[i]=nbp[[i]][2]}
agg2=na.omit(agg)
agg2=as.factor(agg2)
if (graph){
dev.new()
plot(agg2,main="Number of products per group",xlab=c("Number of products"),ylab=c("Frequency"))}

#Nombre de groupes
lev2=as.factor(lev)
if (graph){
dev.new()
plot(lev2,main="Number of groups provided during sorting task",xlab=c("Number of groups"),ylab=c("Frequency"))}


##Analyse textuelle
texte=matrix(NA,(I*J),3)
texte=data.frame(texte)
texte[,1]=rep(rownames(don),J)
texte[,2]=rep(colnames(catego_c),each=I)
for (i in 1:J){
texte[((I*(i-1))+1):(I*i),3]=paste(catego_c[,i])}

restext=textual(texte,3,1,sep.word=sep.word)

#Suppression des modalité g1, ..., g99 (attention tout est mis en minuscule avce textual)
mod.suppr=paste("g",1:99,sep="")
mod.suppr=intersect(colnames(restext$cont.table),mod.suppr)
if (length(mod.suppr)!=0){
num.mod.suppr=which(colnames(restext$cont.table)%in%mod.suppr)
restext$cont.table=restext$cont.table[,-num.mod.suppr]
num.mod.suppr2=which(rownames(restext$nb.words)%in%mod.suppr)
restext$nb.words=restext$nb.words[-num.mod.suppr2,] }


#Nombre de mots différents
nb_mot_diff=nrow(restext$nb.words)
cat("Number of different words : ",nb_mot_diff,"\n")


#Nombre de mots par classe
mots=rep(NA,sum(lev))
grp=0
for (i in 1:J){
mots[(grp+1):(grp+lev[i])]=levels(catego_c[,i])
grp=grp+lev[i]}
mots_split=strsplit(mots,split=sep.word)
nb_mots=rep(NA,length(mots_split))
for (i in 1:length(mots_split)){
if (mots_split[[i]][1] %in% paste("G",1:99,sep="")){
nb_mots[i]=0}
else {
nb_mots[i]=length(mots_split[[i]])}}
nb_mots2=as.factor(nb_mots)

if (graph){
dev.new()
plot(nb_mots2,main="Number of words per group")}

#Seuil minimum à mettre en paramètre...
freq_min=which(apply(restext$cont.table,2,sum)<=word.min)
if (length(freq_min)!=0){
restext$cont.table=restext$cont.table[,-freq_min]}


caract_prod=descfreq(restext$cont.table)
###########################Analyse textuelle



indicator=list(catego=coocc_reord,napping=list(res1_nappe,res2_nappe[[1]]))

call=list(hmfa=afmh,X=don_input)

ind=list(coord=afmh$ind$coord,cos2=afmh$ind$cos2,contrib=afmh$ind$contrib,partial=afmh$partial)

res = list(eig=afmh$eig,ind=ind,quali.var=afmh$quali.var,quanti.var=afmh$quanti.var,group=afmh$group,indicator=indicator,textual=caract_prod,call=call)

class(res) <- c("fasnt", "list ")

return(res)
}
