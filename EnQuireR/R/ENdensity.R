##graphe de densité pour les différents groupes de l'ACM.
##res est le tableau des coordonnées des individus sur les deux premiers axes de l'ACM.
##classe est la variable de classe.

"ENdensity"=function(res,classe=NULL,nlevels=20) 
{
res.ind=res$ind$coord[,1:2]
res2<-cbind(res.ind,classe)

a<-nlevels(classe)
col<-c("black", "red", "green3", "blue", "cyan", "magenta", 
            "darkgray", "darkgoldenrod", "darkgreen", "violet", 
            "turquoise", "orange", "lightpink", "lavender", "yellow", 
            "lightgreen", "lightgrey", "lightblue", "darkkhaki", 
            "darkmagenta", "darkolivegreen", "lightcyan", "darkorange", 
            "darkorchid", "darkred", "darksalmon", "darkseagreen", 
            "darkslateblue", "darkslategray", "darkslategrey", 
            "darkturquoise", "darkviolet", "lightgray", "lightsalmon", 
            "lightyellow", "maroon")
if (is.null(classe)) {
    kde2d_table <- kde2d(res2[,1],res2[,2],lims = c(min(res.ind[,1])*1.5,max(res.ind[,1])*1.5,min(res.ind[,2])*1.5,max(res.ind[,2])*1.5))
    contour(kde2d_table,nlevels=nlevels,col=col[4])
    #Ajout des labels des axes
    title(main="Density curbs",xlab=paste("Dim 1 (",signif(res$eig[1,2],digits=4),"%)"),ylab=paste("Dim 2 (",signif(res$eig[2,2],digits=4),"%)"))

            }

else {
table<-res.ind[(which(res2[,3]==1)),]       
        kde2d_table_groupe<-kde2d(table[,1],table[,2],lims = c(min(res.ind[,1])*1.5, max(res.ind[,1])*1.5, min(res.ind[,2])*1.5, max(res.ind[,2])*1.5))
        contour(kde2d_table_groupe,add=F,nlevels=nlevels,col=col[1])

for (i in 2:a){
        
        table<-res.ind[(which(res2[,3]==i)),]       
        kde2d_table_groupe<-kde2d(table[,1],table[,2],lims = c(min(res.ind[,1])*1.5, max(res.ind[,1])*1.5, min(res.ind[,2])*1.5, max(res.ind[,2])*1.5))
        contour(kde2d_table_groupe,add=T,nlevels=nlevels,col=col[i],xlim=range(c(-10,10)),ylim=range(c(-10,10)))
            }
        #Ajout des labels des axes
        title(main="Density curbs",xlab=paste("Dim 1 (",signif(res$eig[1,2],digits=4),"%)"),ylab=paste("Dim 2 (",signif(res$eig[2,2],digits=4),"%)"))

        }
}



