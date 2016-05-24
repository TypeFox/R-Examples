v.y <-
function(v){
v.star=rep(NA,length(v)-1)
v=v[v>0]
v=logit(v[-1]/v[-length(v)])
v.star[1:length(v)]=v
return(v.star)
}

