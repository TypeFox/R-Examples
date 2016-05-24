fctcd <-
function(data1,gcname,yname){
                               # gcname : est le genotype de l enfant
                               datnmv<-data1[is.na(data1[gcname])!=TRUE,] # donnees non manquantes
                               datmv<-data1[is.na(data1[gcname])==TRUE,] # donnees manquante
                               vecg<-c(0,1,2)
                               vecy<-0:1
                               if(dim(datmv)[1]!=0){
                                                   # construction du genotype de l enfant sur les donnees manquantes 
                                                     gc<-rep(vecg,2*dim(datmv)[1])
                                                     yy<-rep(rep(vecy,rep(length(vecg),2)),dim(datmv)[1])
                                                     datdm<-datmv[rep(1:nrow(datmv),rep(2*length(vecg),nrow(datmv))),]
                                                     datdm$id = rep(1:nrow(datmv),rep(2*length(vecg),nrow(datmv)))
                                                     datdm$vdcop<-as.numeric(datdm[yname]==yy)
                                                     datdm[gcname]<-gc; # table des donnes complete
                                                     datdm[yname]<-yy
                                                   }else{datdm<-NULL}
                               # contruction du nouveau genotype de l enfant sur les donnees complet
                               datnmv0<-datnmv
                               gc<-rep(vecg,2*dim(datnmv0)[1])
                               yy<-rep(rep(vecy,rep(length(vecg),2)),dim(datnmv0)[1])
                               datdmcp<-datnmv0[rep(1:nrow(datnmv0),rep(2*length(vecg),nrow(datnmv0))),]
                               datdmcp$vdcop<-as.numeric(datdmcp[gcname]==gc&datdmcp[yname]==yy)
                               datdmcp[gcname]<-gc
                               datdmcp[yname]<-yy
                               return(list(datnmv=datnmv,datdm=datdm,datmv=datmv,datdmcp=datdmcp))
                               }
