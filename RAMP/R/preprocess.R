preprocess <-
function(trx, tex, size=ncol(trx)){
minx = min(min(trx),min(tex))
if(minx<1){
trx = trx+1-minx
tex = tex+1-minx
}
#trx[trx<1] = 1
#tex[tex<1] = 1

varx = apply(trx,2,var)
ind = sort(varx,index=TRUE, decreasing=TRUE)

trx = trx[,ind$ix[1:size]]
tex = tex[,ind$ix[1:size]]

trx = log(trx)/log(10)
tex = log(tex)/log(10)
trx = scale(trx)
tex = scale(tex,attr(trx,'scaled:center'), attr(trx,'scaled:scale'))

return(list(trx=trx,tex=tex))	
}
