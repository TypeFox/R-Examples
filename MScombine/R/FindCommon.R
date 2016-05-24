#' Find entities presented in both polarities
#' 
#' Takes matrices from positive and negative ionization mode and find entities in common.
#' @param POSITIVE A matrix with positive entities information (Compound Name, Mass, RT, and multiple columns with the area of the compound in samples).
#' @param NEGATIVE A matrix with positive entities information (Compound Name, Mass, RT, and multiple columns with the area of the compound in samples).
#' @param ADDUCTS A matrix with positive adducts, negative adducts and their difference in mass.
#' @param RTtolerance The tolerance of retention time when comparing both polarities. It should be in the same units as the RT in POSITIVE and NEGATIVE matrices.
#' @param Masstolerance The tolerance in Da when considering the adducts that can be present in both matrices.
#' @examples
#' \dontrun{
#' CommonEntities<-FindCommon(POSITIVE,NEGATIVE,ADDUCTS,Masstolerance=0.02,RTtorelance=0.5)
#' }
#' @export 
FindCommon<-function(POSITIVE,NEGATIVE,ADDUCTS,Masstolerance,RTtolerance) {ctr1 = 2;
                                                                           mh = c("CpdID+","CpdID-","Adduct+","Adduct-","Mass+","Mass-","RT+","RT-","Mean+","Mean-","N+","N-","Correlation");
                                                                           x = dim(POSITIVE);
                                                                           npos = x[1];
                                                                           POSITIVE_subset<-POSITIVE[,-1];
                                                                           POSITIVE_subsetX<-POSITIVE_subset[,-1];
                                                                           POSITIVE_subset2<-POSITIVE_subsetX[,-1];
                                                                           mpos<-(dim(POSITIVE_subset2))[2];
                                                                           POS <- as.matrix(POSITIVE_subset2);
                                                                           y = dim(NEGATIVE);
                                                                           nneg = y[1];
                                                                           NEGATIVE_subset<-NEGATIVE[,-1];
                                                                           NEGATIVE_subsetX<-NEGATIVE_subset[,-1];
                                                                           NEGATIVE_subset2<-NEGATIVE_subsetX[,-1];
                                                                           mneg<-(dim(NEGATIVE_subset2))[2];
                                                                           NEG <- as.matrix(NEGATIVE_subset2);
                                                                           z = dim(ADDUCTS);
                                                                           posnegaddu = z[1];
                                                                           cont=Masstolerance;
                                                                           a=1;
                                                                           while(a<(posnegaddu+1)){
                                                                             b=3;
                                                                             valor=ADDUCTS[a,b];
                                                                             i=1;
                                                                             while(i<(npos+1)){ 
                                                                               j=1; 
                                                                               while(j<(nneg+1)){ 
                                                                                 RTd = abs(POSITIVE[i,3]-NEGATIVE[j,3]);
                                                                                 if(RTd <= RTtolerance){
                                                                                   diffmass = POSITIVE[i,2]-NEGATIVE[j,2];
                                                                                   diffmass2=diffmass-valor;  		
                                                                                   if(abs(diffmass2)<cont){
                                                                                     NEGj=as.vector(NEG[j,]);
                                                                                     POSi=as.vector(POS[i,]);
                                                                                     mh=rbind(mh,c(0));
                                                                                     mh[ctr1,1] = POSITIVE[i,1];
                                                                                     mh[ctr1,2] = NEGATIVE[j,1];
                                                                                     mh[ctr1,3] = ADDUCTS[a,1]; 
                                                                                     mh[ctr1,4] = ADDUCTS[a,2];
                                                                                     mh[ctr1,5] = POSITIVE[i,2]; 
                                                                                     mh[ctr1,6] = NEGATIVE[j,2]; 
                                                                                     mh[ctr1,7] = POSITIVE[i,3]; 
                                                                                     mh[ctr1,8] = NEGATIVE[j,3];
                                                                                     mh[ctr1,9] = round(apply(POSITIVE_subset2[i,c(1:mpos)],1,mean,dim=1,na.rm=TRUE),digits=6)
                                                                                     mh[ctr1,10] = round(apply(NEGATIVE_subset2[j,c(1:mneg)],1,mean,dim=1,na.rm=TRUE),digits=6)
                                                                                     mh[ctr1,11] = sum(!is.na(POSi))
                                                                                     mh[ctr1,12] = sum(!is.na(NEGj))
                                                                                     mh[ctr1,13] = cor(POSi, NEGj,use="na.or.complete")
                                                                                     ctr1=ctr1+1;
                                                                                     j=j+1;
                                                                                   }else j=j+1;
                                                                                   
                                                                                 }else j=j+1;
                                                                               }
                                                                               i = i+1;
                                                                             }
                                                                             a=a+1;
                                                                           }
                                                                           
                                                                           write.table(mh,file="CommonEntities.csv",sep=",",row.names=FALSE,col.names=FALSE);
                                                                           CommonEntities<-read.table("CommonEntities.csv",sep=",",header=TRUE,as.is=c(1,2))
                                                                           colnames(CommonEntities) <- c("CpdID+","CpdID-","Adduct+","Adduct-","Mass+","Mass-","RT+","RT-","Mean+","Mean-","N+","N-","Correlation")
                                                                           CommonEntities}