X1X2BFtable <-
function(n1,n2,m=5,w1,w2,pH0=0.5,e1w=0.5,gfun){
               pX1X2<-matrix(0,nrow=n1+1,ncol=n2+1)
               dimnames(pX1X2)<-list(x1=0:n1,x2=0:n2)
               postNull<-matrix(0,nrow=n1+1,ncol=n2+1)
               dimnames(postNull)<-list(x1=0:n1,x2=0:n2)
               BFX1X2<-matrix(0,nrow=n1+1,ncol=n2+1)
               dimnames(BFX1X2)<-list(x1=0:n1,x2=0:n2)
               LogBFX1X2<-matrix(0,nrow=n1+1,ncol=n2+1)
               dimnames(LogBFX1X2)<-list(x1=0:n1,x2=0:n2)

                for(i in 1:(n1+1)){
                   for(j in 1:(n2+1)){
                        x1=i-1
                        x2=j-1
                        
                        
                        pX1X2[i,j]<-mx(x1,x2,n1,n2,w1=w1,w2=w2,m=m)
                        postNull[i,j]<-postNullFun(x1,x2,n1,n2,m=m,w1=w1,w2=w2,gfun=gfun)
                        
                        BFX1X2[i,j]<-((1-postNull[i,j])*pH0)/(postNull[i,j]*(1-pH0))
                        
                        if(1/BFX1X2[i,j]==Inf){LogBFX1X2[i,j]<- -Inf}else{LogBFX1X2[i,j]<-log(BFX1X2[i,j])}
                        #cat(c('i=',i,'j=',j,'BF=',BFX1X2[i,j],'numBF',numBFX1X2[i,j]),fill=T)
                                     }}

                #calculating cut-off the value
                pTHETA0x1x2<-postNull*pX1X2
                pTHETA1x1x2<-(1-postNull)*pX1X2
                    
                    logbf.1d<-as.vector(LogBFX1X2)
                    logbf.1d<-sort(logbf.1d)
            
                    n.logbf.1d<-length(logbf.1d)
                    L0_HB<- -Inf
                    L0_BE<- -Inf
                    
                    minWTE_HB<-Inf
                    minWTE_BE<-Inf

                    i=0
                    WTE_HB=0
                    WTE_BE=0
                    
                    for(i in 1:n.logbf.1d){


                                            L0<-logbf.1d[i]
                                            IndicLBFgtL0<-(LogBFX1X2>L0)
                                            IndicLBFleL0<-(LogBFX1X2<=L0)

                                            WTE_HB<-e1w*sum(IndicLBFgtL0*pTHETA0x1x2)/pH0+(1-e1w)*(1-
                                                  sum(IndicLBFgtL0*pTHETA1x1x2)/(1-pH0))

                                            if(WTE_HB<minWTE_HB){minWTE_HB<-WTE_HB; L0_HB<-L0}
                                            
                                            if(abs(L0)==Inf){WTE_BE=Inf}else{WTE_BE<-e1w*sum(IndicLBFgtL0*pTHETA0x1x2)/sum(IndicLBFgtL0*pX1X2)+
                                            (1-e1w)*sum(IndicLBFleL0*pTHETA1x1x2)/sum(IndicLBFleL0*pX1X2)}

                                            if(WTE_BE<minWTE_BE){minWTE_BE<-WTE_BE; L0_BE<-L0}

                                            }
                
                
                list(pX1X2=pX1X2,postNull=postNull,BFX1X2=BFX1X2,LogBFX1X2=LogBFX1X2,
                     L0_HB=L0_HB,minWTE_HB=minWTE_HB,L0_BE=L0_BE,minWTE_BE=minWTE_BE)}

