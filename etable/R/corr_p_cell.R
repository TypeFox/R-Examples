corr_p_cell <-
function(x, y, z, w, cell_ids, row_ids, col_ids, vnames,  vars, n_min, digits=3){

if(is.null(w))   w<- rep(1, length(x))
out<- ''
if(length(cell_ids)>= n_min){
xin<-x[cell_ids]
yin<-y[cell_ids]
win<-w[cell_ids]
inout<- which(is.na(xin) | is.na(yin) | is.na(win))
if(length(inout)>0){
xin<-xin[-inout]
yin<-yin[-inout]
win<-win[-inout]
}
corrvalue<- cov.wt(cbind(xin, yin), wt = win, cor = TRUE, center = TRUE, method = c("unbiased", "ML"))$cor[1,2]
out<- a.round.ade(corrvalue, digits=digits)
}
return(out)
}
