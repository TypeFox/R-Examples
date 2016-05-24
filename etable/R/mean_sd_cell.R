mean_sd_cell <-
function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min, digits=3, style=1, nsd=1){

out<- ''
if(length(cell_ids)>= n_min){
if(is.null(w))   w<- rep(1, length(x))
w[which(is.na(x))]<-NA

if(style==1){
out<-  paste(
a.format_n.ade(    wtd.mean(x[cell_ids], w[cell_ids], na.rm=TRUE),  digits=digits), ' (',
a.format_n.ade(sqrt(wtd.var(x[cell_ids], w[cell_ids], na.rm=TRUE)), digits=digits),')', sep='')
}
if(style==2){
m<- wtd.mean(x[cell_ids], w[cell_ids], na.rm=TRUE)
msd<-sqrt(wtd.var(x[cell_ids], w[cell_ids], na.rm=TRUE))
out<-  paste(a.format_n.ade(m, digits=digits), ' (', a.format_n.ade(m-msd*nsd, digits=digits)  ,'~',  a.format_n.ade(m+msd*nsd, digits=digits), ')', sep='')
}
if(style==3){
out<-  paste(
a.format_n.ade(    wtd.mean(x[cell_ids], w[cell_ids], na.rm=TRUE),  digits=digits), '\u00b1',
a.format_n.ade(sqrt(wtd.var(x[cell_ids], w[cell_ids], na.rm=TRUE)), digits=digits),  sep='')
}
}
return(out)
}
