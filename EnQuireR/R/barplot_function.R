#Barplot simple

"barplot_function"=function(dataset,X,spl=FALSE,numr=NULL,numc=NULL,cex=1,colour=NULL){      #X=vecteur des variables choisies, numr=nb de lignes de la fenêtre graphique, numc=nb de colonnes de la fenête graphique, col=couleurs
mult=length(X)                                                             #mult=taille de X
col=rep(0)
count_mod<-rep(0)
a<-1/mult
for (m in 1:mult){
    if (spl){cont<-sort(table(dataset[,X[m]]))
    }
    else {cont<-table(dataset[,X[m]])}
         
        NR <- dim(cont) ## on veut le nombre de lignes      
        count_mod<-c(count_mod,NR)
                    ## permet de diviser la palette de couleurs en autant de couleurs qu'on a de colonnes.
            
                
                        max<-max(cont)
                        coli<-rep(0,NR)         ## pour chaque colonne de la matrice, on crée un répertoire qui va répertorier les différentes teintes dégradées de la couleur de la colonne i.
                            for (j in 1:NR){
                                if (cont[j]==max){
                                                coli[j]<-hsv(h=a*m,s=1,v=1,1,1)
                                                
                                            }
                        else {  
                                        coli[j]<-hsv(h=a*m,s=0.4+(cont[j]/max)*0.6,v=1,1,1) ## pour toutes les colonnes, on fait comme expliqué 13 lignes plus haut
                            }
                                        }
                col<-c(col,coli)

                }  
            col<-col[-1]
            count_mod<-count_mod[-1]
            summ<-cumsum(count_mod)
        
na=matrix(0,1,mult)                                                        #création d'une matrice nulle
if (is.null(numr)){                                                        #si numr et numc non précisée, 2 et 2 par défaut
    if (is.null(numc)){
    numr=numc=2}}
par(mfrow=c(numr,numc),yaxt="n")
tpolice <- par("cex") #pour le mtext
par(xpd=T,mar=par()$mar+c(0,1,0,0))
for (m in 1:mult) {
for (i in 1:((dim(dataset)[2])/(numr*numc))){                                        #création de nouvelles fenêtres graphiques si le nb de variables dépasse numr et numc
        if (m == ((numr*numc)*i)+1){ 
          x11()
        par(mfrow=c(numr,numc),yaxt="n")}
    }
k=0                                                                        #initialisation de k
    for (i in 1:length(dataset[,X[m]])){                                     #pour chaque cellule du tableau de données, si NA alors incrémentation de k
        if ((is.na(dataset[,X[m]])[[i]])==TRUE){
        k=k+1}
    }

na[m]=k/length(dataset[,X[m]])                                             #calcul du pourcentage de données manquantes (nb de NA sur nb total de cellules)
        if (spl==TRUE){                                                               #si spl=T, tri par taille de bâton
        coord=barplot(sort(table(dataset[,X[m]])),beside=TRUE,las=2,horiz=TRUE,main=names(dataset)[X[m]],col=col[(summ[m]-count_mod[m]+1:summ[m])])
text(x=2,y=coord,labels=names(sort(table(dataset[,X[m]]))),adj=0,cex=cex,col=colour)}
        if(spl==FALSE){
        coord=barplot(table(dataset[,X[m]]),beside=TRUE,las=2,horiz=TRUE,main=names(dataset)[X[m]],col=col[(summ[m]-count_mod[m]+1:summ[m])])  #sinon, tri par ordre alphabétique
text(x=2,y=coord,labels=levels(dataset[,X[m]]),adj=0,cex=cex,col=colour)}
mtext(paste(c("Percentage of missing values  =",round(na[m],2)*100,"%"),collapse=" "), side = 3, line=-0.1,cex=tpolice,adj=0)   #phrase : % de NA=...
}}

#ENbarplot(tea,c(names(tea[,-22])),spl=TRUE)
#ENbarplot(tea,c(names(tea[,-22])))
#ENbarplot(tea,c("price","sex","socio.professional.category","age_Q"))
#ENbarplot(tea,c("price","sex","socio.professional.category","age_Q"),spl=TRUE)
