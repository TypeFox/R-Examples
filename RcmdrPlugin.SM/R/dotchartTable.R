dotchartTable<-
function(tableau,...){
x<-as.numeric(tableau)
names(x)<-names(tableau)
dotchart(x,...)
}

