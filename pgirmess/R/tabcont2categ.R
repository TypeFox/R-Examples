"tabcont2categ" <-
function(tab){
# Giraudoux 1.8.2004 convert a contingency table into
# a presence/absence table of categories

cnoms<-names(tab)
rnoms<-row.names(tab)
tab3<-NULL
    for (j in 1:length(tab)) {
        for (i in 1:length(tab[,1])) {
        if (tab[i,j]!=0){
            cate<-rep(cnoms[j],tab[i,j])
            rec<-rep(rnoms[i],tab[i,j])
            
            tab2<-cbind(rec,cate)
            tab3<-rbind(tab3,tab2)
            }
        }
    }
tab3<-as.data.frame(tab3)
tab3[,1]<-as.factor(tab3[,1]);tab3[,2]<-as.factor(tab3[,2])
return(tab3)
}

