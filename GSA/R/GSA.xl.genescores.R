GSA.xl.genescores=function(GSA.listsets.obj, genesets,  GSA.obj,  genenames){

o1=as.numeric(GSA.listsets.obj$pos[,1])
pos=vector("list",length(o1))
ii=0
for(i in o1){
ii=ii+1
pos[[ii]]=GSA.genescores(i,genesets,  GSA.obj,  genenames)
}
o2=as.numeric(GSA.listsets.obj$neg[,1])
neg=vector("list",length(o2))
ii=0
for(i in o2){
ii=ii+1
neg[[ii]]=GSA.genescores(i,genesets,  GSA.obj,  genenames, negfirst=TRUE)
}

return(list(posi=o1, negi=o2,pos=pos,neg=neg))
}

