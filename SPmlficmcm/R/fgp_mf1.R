fgp_mf1 <-
function(data,d,gmname,gcname,outc,HW=TRUE){
                   # data : les donnees
                   # gmname : le nom de la variable du genotype de la mere 
                   # gcname : le nom de la variable du genotype de l enfant
                   # d est di=Ni/ni d[1] temoin d[2] cas
                   # outc : est la variable outcome
                   ll<-data[,outc];dof<-data[,c(gmname,gcname)]
                   vd<-ifelse(ll==1,d[2],d[1])
                   if(HW==TRUE){fctp<-function(p){return(sum(vd*logDgmc_HW1(dof,p)))}
                                }else{fctp<-NULL} 
                   
                   return(fctp)
                   }
