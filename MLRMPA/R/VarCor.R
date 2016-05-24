VarCor <-
function(tst,activity,descriptor,deleted_descriptor,var.lower,var.upper,xy.cor){
  activity_01<-scale_MinMax(activity)
  tr.expr<-as.data.frame(activity_01[-tst,1])
  colnames(tr.expr)<-c("expr")
  tst.expr<-as.data.frame(activity_01[tst,1])
  colnames(tst.expr)<-c("expr")
 
  
  descriptor_01<-scale_MinMax(descriptor)
  tr.dscrp<-descriptor_01[-tst,]
  tst.dscrp<-descriptor_01[tst,]
 
 
  HM_delete.dscrp<-deleted_descriptor
  index<-which(colnames(tr.dscrp) %in% HM_delete.dscrp[,1])
  tr.dscrp<-tr.dscrp[,-index]
  dscrp.var<-apply(tr.dscrp,2,var)
  y.var<-var(tr.expr)
  var.range_lower<- var.lower
  var.range_upper<- var.upper
  x_varselect <- tr.dscrp[,which(dscrp.var>=as.numeric(var.range_lower)&dscrp.var<=as.numeric(var.range_upper))]
  x_corselect <- x_varselect[,which(cor(x_varselect,as.numeric(tr.expr[,1]))>=xy.cor|cor(x_varselect,as.numeric(tr.expr[,1]))<=-xy.cor)]
  aa<-vector()
  aa[1]<-dim(x_varselect)[2]
  aa[2]<-dim(x_corselect)[2]
  tr.tst<-list(expr.tr=tr.expr,expr.tst=tst.expr,dscrp.all=tr.dscrp,dscrp.tr=x_corselect,dscrp.tst=tst.dscrp,VarCordim=aa)
  return
  tr.tst  
}
