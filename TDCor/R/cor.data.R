cor.data <-
function(rd0,rtf0)

{  if (is.vector(rtf0))

{rtf0=t(as.matrix(rtf0))

 rownames(rtf0)=c("X")}



if (is.vector(rd0))

{rd0=t(as.matrix(rd0))

 rownames(rd0)=c("X")}



nw_cor=matrix(0,dim(rd0)[1],dim(rtf0)[1])

for (j in 1:dim(nw_cor)[2])

{ nw_cor[,j]<-round(apply(rd0,1,cor,as.numeric(rtf0[j,]),use="complete.obs"),3)}

nw_cor=data.frame(nw_cor)

names(nw_cor)=rownames(rtf0)

rownames(nw_cor)=rownames(rd0)

return(as.matrix(nw_cor))}
