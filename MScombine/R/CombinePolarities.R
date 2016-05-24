#' Combine positive and negative matrices
#' 
#' Take positive and negative matrices and combine them by deleting redundat entities.
#' @param POSITIVE A matrix with positive entities information (Compound Name, Mass, RT, and multiple columns with the area of the compound in samples)
#' @param NEGATIVE A matrix with positive entities information (Compound Name, Mass, RT, and multiple columns with the area of the compound in samples)
#' @param CommonEntitiesFiltered The data set generated with the FilterbyRT function.
#' @examples
#' \dontrun{
#' CombinePolarities(POSITIVE,NEGATIVE,CommonEntitiesFiltered)
#' }
#' @export
CombinePolarities<-function(POSITIVE,NEGATIVE,CommonEntitiesFiltered) {e=2;
                                                                       f=2;
                                                                       IDtoDel=c("CpdIDtoDelete");
                                                                       k=dim(CommonEntitiesFiltered);
                                                                       nentities=k[1];
                                                                       IDtoDelNeg = rbind(IDtoDel,c(0));
                                                                       IDtoDelPos = rbind(IDtoDel,c(0));
                                                                       d=1;
                                                                       while(d<(nentities+1)){
                                                                         if(CommonEntitiesFiltered[d,9]>CommonEntitiesFiltered[d,10]){
                                                                           IDtoDelNeg[e] = CommonEntitiesFiltered[d,2];
                                                                           e=e+1;
                                                                           d=d+1;
                                                                         } 
                                                                         else {
                                                                           IDtoDelPos[f] = CommonEntitiesFiltered[d,1];
                                                                           f=f+1;
                                                                           d=d+1;
                                                                         }
                                                                       }
                                                                       
                                                                       IDtoDelNegDef<-unique(IDtoDelNeg)
                                                                       IDtoDelPosDef<-unique(IDtoDelPos)
                                                                       
                                                                       POSITIVEDef<-POSITIVE[-which(POSITIVE[,1] %in% IDtoDelPosDef),]
                                                                       NEGATIVEDef<-NEGATIVE[-which(NEGATIVE[,1] %in% IDtoDelNegDef),]
                                                                       
                                                                       write.table(POSITIVEDef,file="Positive-Negative.csv",sep=",",row.names=FALSE,na="")
                                                                       write.table(NEGATIVEDef,file="Positive-Negative.csv",sep=",",row.names=FALSE,col.names=FALSE,na="",append=TRUE)}