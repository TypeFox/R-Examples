combi_cell <-
function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min, digits=3, style=1){
out<- NULL
if(!is.null(y)){
x<- y
y<- NULL
}
if(is.numeric(w)) if(length(w)==0) w<-NULL


#####################################
if(style==1){
if(length(unique(x)[which(!is.na(unique(x)))])==1){
out<- paste('N:', n_cell(x=x, y=y, z=z, w=w, cell_ids=cell_ids, row_ids=row_ids, col_ids=col_ids,vnames=vnames, n_min=n_min))
}
if(length(unique(x)[which(!is.na(unique(x)))])==2){
out<- eventpct_cell(x=x, y=y, z=z, w=w, cell_ids=cell_ids, row_ids=row_ids, col_ids=col_ids,
vnames=vnames, n_min=n_min, digits=1, digits2=0, event=2, type=1)
}
if(length(unique(x)[which(!is.na(unique(x)))])>2 & !is.factor(x)){
out<- iqr_cell(x=x, y=y, z=z, w=w, cell_ids=cell_ids, row_ids=row_ids, col_ids=col_ids,
vnames=vnames, n_min=n_min, digits=digits)
}
if(length(unique(x)[which(!is.na(unique(x)))])>2 & is.factor(x)){
out<- mode_cell(x=x, y=y, z=z, w=w, cell_ids=cell_ids, row_ids=row_ids, col_ids=col_ids,vnames=vnames, n_min=n_min)
}
}
#####################################
if(style==2){
if(length(unique(x)[which(!is.na(unique(x)))])==1){
out<- paste('N:', n_cell(x=x, y=y, z=z, w=w, cell_ids=cell_ids, row_ids=row_ids, col_ids=col_ids,vnames=vnames, n_min=n_min))
}
if(length(unique(x)[which(!is.na(unique(x)))])==2){
out<- eventpct_cell(x=x, y=y, z=z, w=w, cell_ids=cell_ids, row_ids=row_ids, col_ids=col_ids,
vnames=vnames, n_min=n_min, digits=1, digits2=0, event=2, type=1)
}
if(length(unique(x)[which(!is.na(unique(x)))])>2 & !is.factor(x)){
out<- mean_sd_cell(x=x, y=y, z=z, w=w, cell_ids=cell_ids, row_ids=row_ids, col_ids=col_ids,
vnames=vnames, n_min=n_min, digits=digits, style=2, nsd=1.96)
}
if(length(unique(x)[which(!is.na(unique(x)))])>2 & is.factor(x)){
out<- mode_cell(x=x, y=y, z=z, w=w, cell_ids=cell_ids, row_ids=row_ids, col_ids=col_ids,vnames=vnames, n_min=n_min)
}
}
#####################################

return(out)
}
