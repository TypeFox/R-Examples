calculate.xtab<-function(v1,v2,varnames=NULL) {
 counts<-table(v1,v2)
 row.margin<-apply(counts,1,sum)
 col.margin<-apply(counts,2,sum)
 vl1<-attr(v1,"value.labels")
 if(!is.null(vl1) && length(vl1)) {
  vl1<-sort(vl1)
  if(length(rownames(counts)) < length(vl1)) 
   vl1<-vl1[vl1 %in% rownames(counts)]
  if(is.null(names(vl1))) names(row.margin)<-vl1
  else names(row.margin)<-names(vl1)
 }
 vl2<-attr(v2,"value.labels")
 if(!is.null(vl2) && length(vl2)) {
  vl2<-sort(vl2)
  if(length(colnames(counts)) < length(vl2))
   vl2 <- vl2[vl2 %in% colnames(counts)]
  if(is.null(names(vl2))) names(col.margin)<-vl2
  else names(col.margin)<-names(vl2)
 }
 if(is.null(varnames)) {
  rowname<-deparse(substitute(v1))
  colname<-deparse(substitute(v2))
 }
 else {
  rowname<-varnames[1]
  colname<-varnames[2]
 }
 xt<-list(counts=counts,row.margin=row.margin,col.margin=col.margin, 
  varnames=c(rowname,colname))
 attr(xt, "class") <- "xtab"
 return(xt)
}
