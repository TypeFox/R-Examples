quantile_cell <-
function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min, digits=3, probs=0.5, plabels=FALSE){

out<- ''
if(length(cell_ids)>= n_min){
if(all(w==1)) w<- NULL
for(i in 1:length(probs)){

if(!is.null(w)){
out_part <-  paste(a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=probs[i]), digits=digits),sep='')
}
if(is.null(w)){
out_part<-  paste(a.format_n.ade(quantile(x[cell_ids], probs=probs[i], na.rm=TRUE, type=8), digits=digits),sep='')
}
if(i==1 & !plabels) out<- out_part
if(i==1 & plabels)  out<- paste(((probs[i])*100),'% (',out_part, ')', sep='')
if(i>1 & !plabels)  out<- paste(out, ' | ', out_part, sep='')
if(i>1 & plabels)   out<- paste(out, ' | ',((probs[i])*100),'% (', out_part, ')', sep='')
}
}

return(out)
}
