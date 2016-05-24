update_designinlist <- function(designsin,groupsize,ni,xt,x,a,it,iGroup){

design$xt = xt
design$a = a
design$x = x
design$groupsize = groupsize
design$ni = ni

if((it==-1) ){
    design$num = length(designsin)+1
} else {
    design$num = it
}
design$iGroup = iGroup
    
designsin=matrix(c(designsin, design),nrow=1,byrow=T)

return( designsin ) 
}
