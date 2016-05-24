DocVarTable <-
function(base,VarSel)
{
     if(is.character(VarSel))
      VarSel<-which(colnames(base)%in%VarSel)
      if(is.numeric(VarSel))
	  VarSel<-VarSel
      DocVar<-as.data.frame(base[,VarSel])
      rownames(DocVar)<-rownames(base)
      if(length(VarSel)==1)
      colnames(DocVar)<-colnames(base[VarSel])
return(DocVar)
 }
