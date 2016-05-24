miss_cell <-
function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min, pct=FALSE, digits=0, prefix='', suffix=''){
n<- sum(is.na(x[cell_ids]))
N<- length(cell_ids)
if(pct)  out<- paste(prefix, a.round.ade(n, digits), ' (' ,a.round.ade((n/N)*100, 1), '%)', suffix, sep='')
if(!pct) out<- paste(prefix, a.round.ade(n, digits),              suffix, sep='')
return(out)
}
