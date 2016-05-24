maxi.false <-
function(vector){
if((sum(is.na(vector))<length(vector))&(sum(sum(is.na(vector)),sum(vector,na.rm=TRUE))<length(vector))){
vector[1:which(!vector)[length(which(!vector))]]=FALSE
vector
}
vector
}
