return1colmatrix <- function(x,grp.name="",sortx=FALSE) {
#Converts vectors, data frames, lists, arrays etc. to 1-column matrix with optional column name for use in orddom functions
varname<-match.call()$x
if (is.vector(x)) { x<-t(as.matrix(x)) } 
if (is.list(x)==TRUE&&is.data.frame(x)==FALSE) {
 x<-as.data.frame(x[[1]]) }
if (is.data.frame(x)) { 
  if(nrow(x)==1&&ncol(x)>1) {x<-t(as.matrix(x))} else
 {if(nrow(x)>1&&ncol(x)==1) {x<-as.matrix(x)} else
 {if(nrow(x)>1&&ncol(x)>1)  {name_col<-names(x)
  x<-as.matrix(x[,1])
  colnames(x) <- name_col[1]
 }}}}
if (is.array(x)==TRUE&&is.matrix(x)==FALSE) {
  att<-attributes(x)
  if(length(att$dim)>=2){
  x<-matrix(x,att$dim[1],att$dim[2],dimnames=list(att$dimnames[[1]],att$dimnames[[2]]))}
  else{
  x<-matrix(x,att$dim[1],dimnames=list(att$dimnames[[1]]))}}
if (is.matrix(x)) {
 if(nrow(x)==1&&ncol(x)>1) {x<-t(as.matrix(x))} else
 {if(nrow(x)>1&&ncol(x)==1){x<-as.matrix(x)}}
 if(nrow(x)>1&&ncol(x)>1)  {x<-x[,1,drop=FALSE]}
if (is.null(colnames(x))) { name_col<-c(varname) } else { name_col<-c(colnames(x)[1]) }
if(sortx==TRUE)
  {x <- t(matrix(as.numeric(x[,1][order(x[,1])]),1))} else { x<-t(matrix(as.numeric(x[,1]),1)) }
if(grp.name!=""){name_col<-grp.name}  
colnames(x) <- name_col[1]
return (x)} else{ return (NA) }
}

