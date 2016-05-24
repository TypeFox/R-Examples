fitdi <-
function(z,T,xint=NA,itype=c("step","pulse"),delta=NA,delta0=NA,...){

  itype=match.arg(itype)
  ###Check whether delta exists in the model##
  if (exists("delta",mode="numeric")==TRUE){
    Mod.Out<-GetIntMod(delta=delta,z=z,T=T,xint=xint,itype=itype,...)
    output<-list("delta"=delta,"Int.Mod"=Mod.Out$Model,"xi"=Mod.Out$xi)  
  }

   else 
     {
       ans <- optimize(f = GetLikeDI, interval = c(0.1, 2), maximum = TRUE,z=z,T=T,xint=xint,itype=itype,...)
       delta1<-ans$maximum
       Mod.Out<-GetIntMod(delta=delta1,z=z,T=T,xint=xint,itype=itype,...)
       output<-list("delta"=delta1,"Int.Mod"=Mod.Out$Model,"xi"=Mod.Out$xi)  
     }
  if (exists("delta0",mode="numeric")==TRUE)
    {
    pval<-GetPV(delta0=delta0,z=z,T=T,xint=xint,itype=itype,...)
    output<-list("delta"=output$delta,"Int.Mod"=output$Int.Mod,"xi"=output$xi,"pvalue"=pval)
  }
  class(output)<-"dynia"
  return(output)
  }


