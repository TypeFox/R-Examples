iqr_cell <-
function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min, digits=3, add_n=FALSE){
out<- ''
if(length(cell_ids)>= n_min){
if(all(w==1)) w<- NULL
if(is.null(w)){
out<-  paste(
a.format_n.ade(quantile(x[cell_ids], probs=c(0.5),  na.rm=TRUE, type=8), digits=digits), ' (',
a.format_n.ade(quantile(x[cell_ids], probs=c(0.25), na.rm=TRUE, type=8), digits=digits), '/',
a.format_n.ade(quantile(x[cell_ids], probs=c(0.75), na.rm=TRUE, type=8), digits=digits), ')',
sep='')
if(add_n) out<-  paste(out, 'N:', length(cell_ids))
}
if(!is.null(w)){

out<-  paste(
a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.5)),  digits=digits), ' (',
a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.25)), digits=digits), '/',
a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.75)), digits=digits), ')',
sep='')
if(add_n) out<-  paste(out, 'N:', round(sum(w[cell_ids][which(!is.na(x[cell_ids]))], na.rm=TRUE)))
}
}
return(out)
}
