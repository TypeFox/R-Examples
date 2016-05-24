eventpct_cell <-
function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min, digits=1, digits2=0, event=2, type=1){
if(is.null(w))   w<- rep(1, length(x))
out<-''
x<- as.factor(x)
w[which(is.na(x))]<- NA
N<- sum(w[cell_ids], na.rm=TRUE)
xin<- x[cell_ids]
win<- w[cell_ids]
events<- sum(win[which(xin==levels(xin)[event])], na.rm=TRUE)
pct<- paste(a.round.ade((events/N)*100, digits), '%', sep='')
n<-a.round.ade(events, digits2)
if(type==1) out<- paste(pct, ' (', n, ')', sep='')
if(type==2) out<- paste(n, ' (', pct, ')', sep='')
if(type==3) out<- paste(pct , sep='')
if(type==4) out<- paste(n , sep='')
if(type==5) out<- paste(pct, ' (', n, '/', a.round.ade(N, digits2) ,')', sep='')
return(out)
}
