norm.data <-
function(rd_sub,...)

{ Lst<-list(...)

  if (!is.null(Lst$u))

  {u<-Lst$u}else

  { if(is.vector(rd_sub))

    { u<-seq(1,length(rd_sub))}else

    { u<-seq(1,dim(rd_sub)[2])}}



if (is.vector(rd_sub))

{rd_sub_norm=(rd_sub[u]-min(rd_sub[u]))/(max(rd_sub[u])-min(rd_sub[u]))}

 else

{ k_max=apply(rd_sub[,u],1,max)

  k_min=apply(rd_sub[,u],1,min)

  rd_sub_norm=(rd_sub-k_min)/(k_max-k_min)}

  return(rd_sub_norm)}
