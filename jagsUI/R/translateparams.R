
translate.params <- function(x,params.sub){


params = colnames(x$samples[[1]])

params.simple.sub = unique(sapply(strsplit(params.sub, "\\["), "[", 1))
params.simple <- unique(sapply(strsplit(params, "\\["), "[", 1))
n = length(params.simple.sub)

if(sum(params.simple.sub%in%params.simple)!=n){stop('One or more specified parameters are not in model output./n')}


params.sub.1 <- sapply(strsplit(params.sub, "\\]"), "[", 1)
params.2 <- sapply(strsplit(params.sub.1, "\\["), "[", 2)
expand <- sapply(strsplit(params, "\\["), "[", 1)

dim = get.dim(params)

gen.samp.mat <- function(x){
  out = x
  for(i in 1:length(x)){
    if(!is.na(x[[i]][1])){
      if(length(x[[i]])>1){
        out[[i]] = array(params[expand==names(x)[i]],dim=x[[i]])
    }
    if(length(x[[i]])==1){
      out[[i]] = params[expand==names(x)[i]]
    }
    
  } else {out[[i]] = NA}
}
return(out)
}

mats = gen.samp.mat(dim)

mats.sub = mats[params.simple.sub]

index=1
params.new = character()
for (i in 1:length(params.sub)){

  if(!is.na(mats.sub[i])||!is.na(params.2[i])){
    if(params.sub[i]==params.simple.sub[i]){
      st = paste('mats.sub$',params.simple.sub[i],"[]",sep="")
    } else {
      st = paste('mats.sub$',params.simple.sub[i],"[",params.2[i],']',sep="")
    }
    ind = eval(parse(text=st))
    params.new[index:(index+length(ind)-1)] = ind
    index = index+length(ind)
  } else {
    params.new[index]=params.sub[i]
    index=index+1
  }
}
return(params.new)
}
