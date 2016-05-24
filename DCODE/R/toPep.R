toPep <-
function(vec, frame=0){
  
    res=apply(t(apply(vec,1,s2c)),1,translate,frame=frame)
  
    if(!(is.null(dim(res)))){
        res=apply(res,2,c2s)
    }
    return(res)
}
