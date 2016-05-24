n_cell <-
function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min, digits=0, digits2=1, type='n'){
if(is.null(w))   w<- rep(1, length(x))
w[which(is.na(x))]<-NA

if(type=='n')    out<- a.round.ade(sum(w[cell_ids], na.rm=TRUE), digits)
if(type=='pct')  out<- paste(a.round.ade((sum(w[cell_ids], na.rm=TRUE)/sum(w, na.rm=TRUE))*100, digits2), '%', sep='')
if(type=='pctn') out<- paste(a.round.ade((sum(w[cell_ids], na.rm=TRUE)/sum(w, na.rm=TRUE))*100, digits2), '% (', a.round.ade(sum(w[cell_ids], na.rm=TRUE),digits), ')', sep='')
if(type=='rowpct') out<- paste(a.round.ade((sum(w[cell_ids], na.rm=TRUE)/sum(w[row_ids], na.rm=TRUE))*100, digits2), '%', sep='')
if(type=='colpct') out<- paste(a.round.ade((sum(w[cell_ids], na.rm=TRUE)/sum(w[col_ids], na.rm=TRUE))*100, digits2), '%', sep='')
if(type=='rowpctn') out<- paste(a.round.ade((sum(w[cell_ids], na.rm=TRUE)/sum(w[row_ids], na.rm=TRUE))*100, digits2), '% (', a.round.ade(sum(w[cell_ids], na.rm=TRUE),digits), ')', sep='')
if(type=='colpctn') out<- paste(a.round.ade((sum(w[cell_ids], na.rm=TRUE)/sum(w[col_ids], na.rm=TRUE))*100, digits2), '% (', a.round.ade(sum(w[cell_ids], na.rm=TRUE),digits), ')', sep='')
if(type=='all'){
out<- paste(
a.round.ade((sum(w[cell_ids], na.rm=TRUE)/sum(w, na.rm=TRUE))*100, digits2), '% - ',
a.round.ade((sum(w[cell_ids], na.rm=TRUE)/sum(w[row_ids], na.rm=TRUE))*100, digits2), '% row - ',
a.round.ade((sum(w[cell_ids], na.rm=TRUE)/sum(w[col_ids], na.rm=TRUE))*100, digits2), '% col', sep='')
}

return(out)
}
