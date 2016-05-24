mode_cell <-
function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min, digits=3){
if(is.null(w))   w<- rep(1, length(x))
w[which(is.na(x))]<-NA

if(is.factor(x))    v<- factor(x[cell_ids])
if(is.numeric(x))   v<- x[cell_ids]
if(is.character(x)) v<- x[cell_ids]
ul <- sort(unique(v))
tab <- wtd.table(v, w[cell_ids], type='table')
on <- which(tab==max(tab, na.rm=TRUE))
value<- ul[on]
N<- tab[on]
if(is.numeric(x))  out<- paste(a.format_n.ade(value,digits), ' (N:', round(N) , ')',  sep='', collapse='| ')
if(!is.numeric(x)) out<- paste(value, ' (N:', round(N) , ')',  sep='', collapse='|')

return(out)
}
