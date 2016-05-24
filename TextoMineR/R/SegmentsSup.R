SegmentsSup <-
function(DocVar,Table,Tab.SegR,SegSelec=NULL, cos2=NULL, 
ncp=5,num.agg=NULL, graph=TRUE){
ncp <- min(ncp, (nrow(Table) - 1), (ncol(Table) - 1))

if(!is.null(SegSelec))
Tab.SegR<-Tab.SegR[,SegSelec]
if(!is.null(DocVar)) 
base<-DocVar
if(!is.null(num.agg)){ 
  if(is.character(num.agg))
      num.agg<-which(colnames(base)%in%num.agg)
  if(is.numeric(num.agg))
     num.agg<-num.agg
	agg<-(base[,num.agg])
	dis.X<-tab.disjonctif(agg)
	Tagreg.LexSeg<-t(Tab.SegR)%*%dis.X
	Tagreg.LexSeg<-t(Tagreg.LexSeg) 
	Table.Seg<-cbind(Table,Tagreg.LexSeg)
  }else{
      Tagreg.LexSeg=NULL
      Table<-Table[apply(Table,1,sum)>0,]      
      Tab.SegR<-Tab.SegR[which(rownames(Tab.SegR)%in%rownames(Table)),]
      Table.Seg<-cbind(Table,Tab.SegR)
 }
	n<-ncol(Table)+1
	m<-ncol(Table.Seg)
	res.ca<-CA(Table.Seg[apply(Table.Seg,1,sum)>0,],col.sup=c(n:m),ncp, graph=F)
if(graph){
if(!is.null(cos2)){
Cosen<-paste("cos2", cos2, sep = " ")        
plot.CA(res.ca,selectCol=Cosen,unselect=1,invisible=c("col","row","row.sup"),title="Repeated Segments")
dev.new()
plot.CA(res.ca,selectCol=Cosen,unselect=1)
}else{
plot.CA(res.ca,invisible=c("col","row","row.sup"),title="Repeated Segments")
dev.new()
plot.CA(res.ca)
}
}
res=(list(Table.Seg=Table.Seg, TagSeg=Tagreg.LexSeg, res.ca=res.ca))
return(res)
}
