#' Remove mismatched entities
#' 
#' Remove those entities that have been associated to more than one adduct, retaining only the most probable.
#' @param CommonEntities (Matrix generated from the FindCommon function).
#' @examples
#' \dontrun{
#' CommonEntitiesImproved<-RemoveMismatch(CommonEntities)
#' }
#' @export
#' @return CommonEntitiesImproved The matrix without mismatched entities.
RemoveMismatch<-function(CommonEntities) {colnames(CommonEntities)[colnames(CommonEntities)=="Adduct+"] <- "Adduct."
                                          colnames(CommonEntities)[colnames(CommonEntities)=="Adduct-"] <- "Adduct..1"
                                          colnames(CommonEntities)[colnames(CommonEntities)=="RT+"] <- "RT."
                                          colnames(CommonEntities)[colnames(CommonEntities)=="RT-"] <- "RT..1"
                                          colnames(CommonEntities)[colnames(CommonEntities)=="Mass+"] <- "Mass."
                                          colnames(CommonEntities)[colnames(CommonEntities)=="Mass-"] <- "Mass..1"
                                          colnames(CommonEntities)[colnames(CommonEntities)=="CpdID+"] <- "CpdID."
                                          colnames(CommonEntities)[colnames(CommonEntities)=="CpdID-"] <- "CpdID..1"
                                          colnames(CommonEntities)[colnames(CommonEntities)=="Mean+"] <- "Mean."
                                          colnames(CommonEntities)[colnames(CommonEntities)=="Mean-"] <- "Mean..1"
                                          colnames(CommonEntities)[colnames(CommonEntities)=="N+"] <- "N."
                                          colnames(CommonEntities)[colnames(CommonEntities)=="N-"] <- "N..1"
                                          
                                          CpdID.=NULL
                                          Adduct.=NULL
                                          CpdID..1=NULL
                                          Adduct..1=NULL
                                          colp <- c('CpdID.','Adduct.')
                                          coln <- c('CpdID..1','Adduct..1')
                                          CommonEntities$IDp<-apply(CommonEntities[,colp],1,paste,collapse="-")
                                          CommonEntities$IDn<-apply(CommonEntities[,coln],1,paste,collapse="-")
                                          
                                          nPOS<-ddply(CommonEntities,.(CpdID., Adduct.),summarise,nPos=sum(!is.na(CpdID.))); #Lista de Compuestos y aductos con su n, es decir, si un compuesto sigue saliendo repetido es porque tiene varios tipos de aductos, esos son los casos que hay que estudiar
                                          nPOSmustbedeletedX<-ddply(nPOS,"CpdID.",summarise,nPos=sum(!is.na(CpdID.)));
                                          nPOSmustbedeleted<-nPOSmustbedeletedX[nPOSmustbedeletedX$nPos>1,]; #Saca la lista de los compuestos que ademas de estar repetidos tienen distintos tipos de aductos
                                          if(dim(nPOSmustbedeleted)[2]!=0){
                                          Cpdtodel1<-nPOSmustbedeleted[,"CpdID."];
                                          Cpdtodel2<-nPOS[nPOS$CpdID. %in% Cpdtodel1[],];
                                          Cpdtodel4<-Cpdtodel2[0,]
                                          lengthCpdtodel1<-length(Cpdtodel1); #Contiene el numero de compuestos que estaban repetidos y cuyos aductos se encuentran en Cpdtodel2
                                          dd=1;
                                          while(dd<(lengthCpdtodel1+1)){
                                            Cpdtodel3 <- Cpdtodel2[Cpdtodel2$CpdID.==Cpdtodel1[dd],];
                                            if(any(Cpdtodel3[,"nPos"]!=max(Cpdtodel3$nPos)))
                                            {Cpdtodel3 [Cpdtodel3==max(Cpdtodel3$nPos)] <- NA;
                                             Cpdtodel3 <- na.omit(Cpdtodel3);
                                             Cpdtodel4 <- rbind(Cpdtodel3,Cpdtodel4);
                                             dd <- dd+1;
                                            }
                                            
                                            else {
                                              if("M+H" %in% Cpdtodel3[,"Adduct."]){
                                                Cpdtodel3[Cpdtodel3=="M+H"] <- NA;
                                                Cpdtodel3 <- na.omit(Cpdtodel3);
                                                Cpdtodel4 <- rbind(Cpdtodel3,Cpdtodel4);
                                                dd <- dd+1;
                                              }
                                              else {Cpdtodel3[Cpdtodel3=="M+Na"] <- NA;
                                                    Cpdtodel3 <- na.omit(Cpdtodel3);
                                                    Cpdtodel4 <- rbind(Cpdtodel3,Cpdtodel4);
                                                    dd <- dd+1;
                                              }
                                            }
                                          }
                                          Cpdtodel4$IDp<-apply(Cpdtodel4[,colp],1,paste,collapse="-")
                                          CommonEntitiesPreimproved <- CommonEntities[-which(CommonEntities$IDp %in% Cpdtodel4$IDp),]}
                                          else{CommonEntitiesPreimproved<-CommonEntities}
                                          
                                          nNEG<-ddply(CommonEntities,.(CpdID..1, Adduct..1),summarise,nNeg=sum(!is.na(CpdID..1))); #Lista de Compuestos y aductos con su n, es decir, si un compuesto sigue saliendo repetido es porque tiene varios tipos de aductos, esos son los casos que hay que estudiar
                                          nNEGmustbedeletedX<-ddply(nNEG,"CpdID..1",summarise,nNeg=sum(!is.na(CpdID..1)));
                                          nNEGmustbedeleted<-nNEGmustbedeletedX[nNEGmustbedeletedX$nNeg>1]; #Saca la lista de los compuestos que ademas de estar repetidos tienen distintos tipos de aductos
                                          if(dim(nNEGmustbedeleted)[2]!=0){
                                          Cpdtodel1neg<-nNEGmustbedeleted[,"CpdID..1"];
                                          Cpdtodel2neg<-nNEG[nNEG$CpdID..1 %in% Cpdtodel1neg[],];
                                          Cpdtodel4neg<-Cpdtodel2neg[0,]
                                          lengthCpdtodel1neg<-length(Cpdtodel1neg); #Contiene el numero de compuestos que estaban repetidos y cuyos aductos se encuentran en Cpdtodel2
                                          ff=1;
                                          while(ff<(lengthCpdtodel1neg+1)){
                                            Cpdtodel3neg <- Cpdtodel2neg[Cpdtodel2neg$CpdID..1==Cpdtodel1neg[ff]];
                                            if(any(Cpdtodel3neg[,"nNeg"]!=max(Cpdtodel3neg$nNeg)))
                                            {Cpdtodel3neg [Cpdtodel3neg==max(Cpdtodel3neg$nNeg)] <- NA;
                                             Cpdtodel3neg <- na.omit(Cpdtodel3neg);
                                             Cpdtodel4neg <- rbind(Cpdtodel3neg,Cpdtodel4neg);
                                             ff <- ff+1;
                                            }
                                            
                                            else {
                                              if("M-H" %in% Cpdtodel3neg[,"Adduct..1"]){
                                                Cpdtodel3neg[Cpdtodel3neg=="M-H"] <- NA;
                                                Cpdtodel3neg <- na.omit(Cpdtodel3neg);
                                                Cpdtodel4neg <- rbind(Cpdtodel3neg,Cpdtodel4neg);
                                                ff <- ff+1;
                                              }
                                              else {Cpdtodel3neg[Cpdtodel3neg=="M+FA-H"] <- NA;
                                                    Cpdtodel3neg <- na.omit(Cpdtodel3neg);
                                                    Cpdtodel4neg <- rbind(Cpdtodel3neg,Cpdtodel4neg);
                                                    ff <- ff+1;
                                              }
                                            }
                                          }
                                                                                    
                                          Cpdtodel4neg$IDn<-apply(Cpdtodel4neg[,coln],1,paste,collapse="-")
                                          
                                          CommonEntitiesImproved <- CommonEntitiesPreimproved[-which(CommonEntitiesPreimproved$IDn %in% Cpdtodel4neg$IDn),]}
                                          else{CommonEntitiesImproved<-CommonEntitiesPreimproved}
                                          
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="Adduct."] <- "Adduct+"
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="Adduct..1"] <- "Adduct-"
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="RT."] <- "RT+"
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="RT..1"] <- "RT-"
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="Mass."] <- "Mass+"
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="Mass..1"] <- "Mass-"
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="CpdID."] <- "CpdID+"
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="CpdID..1"] <- "CpdID-"
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="Mean."] <- "Mean+"
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="Mean..1"] <- "Mean-"
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="N."] <- "N+"
                                          colnames(CommonEntitiesImproved)[colnames(CommonEntitiesImproved)=="N..1"] <- "N-"
                                          
                                          CommonEntitiesImproved$RTdiff <- CommonEntitiesImproved[,'RT+']-CommonEntitiesImproved[,'RT-']
                                          
                                          write.table(CommonEntitiesImproved[,-(14:15)],file="CommonEntitiesImproved.csv",sep=",",row.names=TRUE,na="")
                                          CommonEntitiesImproved }
