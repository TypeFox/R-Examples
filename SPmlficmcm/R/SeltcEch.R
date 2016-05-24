SeltcEch <-
function(outc,n1,n0,id,datf){
                                  # n0 : le nombre de controles et n1 le nombre des cas
                                  # n1 : le nombre des cas
                                  # id : la variable identifiant
                                  id1<-sample(t(datf[datf[outc]==1,][id]),n1,replace=TRUE)
                                  id2<-sample(t(datf[datf[outc]==0,][id]),n0,replace=TRUE)
                                  ntab1<-NULL
                                  for(u in id1){ntab1<-rbind(ntab1,datf[datf[id]==u,])}
                                  ntab2<-NULL
                                  for(u in id2){ntab2<-rbind(ntab2,datf[datf[id]==u,])}
                                  datR<-rbind(ntab1,ntab2);obs<-c(1:dim(datR)[1])
                                  datR$obs<-NULL;datR1<-data.frame(obs,datR)
                                  return(datR1)
                                  }
