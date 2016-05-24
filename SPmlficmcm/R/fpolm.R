fpolm <-
function(dets,vec,vars,nom,garder=NULL){
                                       # dets : la table 
                                       # vec : vecteur de variable a fixe,
                                       # vars :  variables a sommer
                                       # garder : variables a garder tel quel
                                       # nom : nom de la nouvelle var
                                       n<-length(vec)
                                       ux=0
                                       for(i in 1:n){xs=max(dets[vec[i]])
                                                     if(xs!=0){ux<-ux+t(t(dets[,vec[i]]/xs))*10^{n+1-i}
                                                              }else{ux<-ux}
                                                     }
                                       ux<-round(ux,digits=4)
                                       uxf<-duplicated(ux)
                                       dats1<-data.frame(dets[,c(vec,garder)],ux,uxf)
                                       names(dats1)<-c(vec,garder,"ux","uxf")
                                       ndats1<-dats1[dats1["uxf"]==FALSE,]
                                       ndets1<-ndats1[,c(vec,garder,"ux")]                          
                                       #nvt<-dets1[duplicated(dets1)==FALSE,] 
                                       names(ndets1)<-c(vec,garder,"ux")
                                       vs<-tapply(dets[,vars[1]],ux,sum)
                                       id<-as.numeric(names(vs))
                                       for (j in 2:length(vars))              
                                         vs<-cbind(vs,tapply(dets[,vars[j]],ux,sum))
                                       rownames(vs)<-NULL
                                       tab<-data.frame(id,vs)
                                       names(tab)<-c("ux",nom)
                                       tab2<-merge(ndets1,tab,by=c("ux","ux"))
                                       tab2$ux<-NULL
                                       return(tab2)
                                      }
