identifypathway<-function(componentList,PSS,pathType="KEGG",method="MPINet",
           weightnum=6,backgroundcid=Getenvir("getBackground"),annlim=1,bglim=6,order="pvalue",decreasing=FALSE){
          
 
   Veg<-Getenvir("getnodeseq")
   pathType<-unique(pathType)
  if((length(pathType)==1)&(length(pathType[which(pathType=="KEGG")])>0))
     {
        pathList<-Getenvir("getpathList('KEGG')")
             
        if(method=="MPINet")
         {
           componentList<-unique(intersect(componentList,backgroundcid))
           message(c("your input componentList have ",length(componentList)," components in background"))
           compoundListinnet<-intersect(Veg,componentList)
           message(c("your input componentList have ",length(compoundListinnet)," components in network"))
           
           annList<-list()
           N<-length(unique(backgroundcid)) ###number of background compounds
           
            if(length(pathList[,1])>0)
	         {
               for(i in 1:length(pathList[,1]))
                  {
                     ann<-list(pathwayId=character(),pathwayName="not known",annComponentList=character(),annComponentNumber=0,
                               annBgComponentList=character(),annBgNumber=0,componentNumber=0,bgNumber=0,pvalue=1,fdr=1,
                               InWeight=0,weight=1,anncompinNetworkNum=0,anncompinNetworkList=character(),riskcompinNetworkNum=0,riskcompinNetworkList=character())
                      
                      ann$pathwayId<-pathList[i,1]  ###pathwayId
                      ann$pathwayName<-pathList[i,2]###pathwayName
                      ann$annBgComponentList<-pathList[i,3]###annBgComponentList
                      ann$bgNumber<-N                      ### bgNumber
                      ann$componentNumber<-length(componentList) ###componentNumber
                      
                      tpathmetabolite<-strsplit(pathList[i,"annBgComponentList"],';')[[1]]####all metabolites in pathway i
                      m1<-length(unique(unlist(tpathmetabolite)))
                      inter<-unique(intersect(as.character(componentList),unlist(tpathmetabolite)))
                      x<-length(unique(inter))
                      
                      ann$annBgNumber<-m1                   ####annBgNumber
                      ann$annComponentNumber<-x             ####annComponentNumber
                      ann$annComponentList<-inter           ####annComponentList
                      
                  #######################################calcualte the anncompinNetworkNum,anncompinNetworkList,riskcompinNetworkNum,riskcompinNetworkList of pathway i ##########################
                      inter_2<-intersect(unique(inter),Veg)######## differential metabolites annotate in pathway i and network
                      riskinnet<-intersect(componentList,Veg)#######differential metabolites in network
                      
                      ann$anncompinNetworkNum<-length(inter_2) ###anncompinNetworkNum
                      ann$anncompinNetworkList<-inter_2        ###anncompinNetworkList
                      ann$riskcompinNetworkNum<-length(riskinnet) ###riskcompinNetworkNum
                      ann$riskcompinNetworkList<-riskinnet        ###riskcompinNetworkList
                      
                 #######################calculate pathway Inweight################
                      
                                       strinner<-0
                                       
                                       metaboliteinpath<-as.character(unlist(tpathmetabolite))
                                       metapathinnet<-intersect(Veg,unlist(tpathmetabolite))
                                       
                                       if(length(metapathinnet)>0)
                                           {
                                                pathw<-sum(PSS[metapathinnet,"pss"])
                                                strinner<-pathw/length(metaboliteinpath)
                                       
                                               
                                            }
                                       
                                       
                                       
                                       
                                       
                   
                     
                     ann$InWeight<-strinner ###InWeight
                     
                     
                      
                      
                     annList[[i]]<-ann
                  }#####for i in 1:length(pathList[,1])
                  
                  
                  #######################calculate pathway weight################
                  annList<-annList[sapply(annList,function(x) x$annComponentNumber>(annlim-1))]
                  
         if(length(annList)>0)
          {
                annList<-annList[sapply(annList,function(x) x$annBgNumber>(bglim-1))]
             
            if(length(annList)>0)
            {
             InWeightList<-c()
             for(j in (1:length(annList)))
               {
                   InWeightList<-c(InWeightList,annList[[j]]$InWeight)
               
               }
               meanweight<-mean(InWeightList)
              
               for(j in (1:length(annList)))
               {
                    
                         w<-annList[[j]]$InWeight/meanweight
                         
                           annList[[j]]$weight<-w^weightnum     
                          
                 }
                  
              ############################### calculate Pvalue###########################
                       
                    
                    for(i in (1:length(annList)))
                       {
             
                          tv<-as.numeric(annList[[i]]["weight"])
                          
                          if(tv>0)
                              {
                               
                                pvalue<- pWNCHypergeo((as.numeric(annList[[i]]["annComponentNumber"])-1),as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["bgNumber"])-as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["componentNumber"]),tv,precision=1E-100,lower.tail=FALSE)
                                annList[[i]]["pvalue"]<-pvalue
                               
                              }
                          
                        }####### for(i in (1:length(annList)))
            
               }###########if(length(annList)>0)
            
            
            
            
            }###############if(length(annList)>0)
           
           }####if(length(pathList[,1])>0)
        
        
        
        
        }#####if method=="MPINet"   
        
   if(method=="Hyper")
         {
            componentList<-unique(intersect(componentList,backgroundcid))
            message(c("your input componentList have ",length(componentList)," components in background"))
            annList<-list()
            N<-length(unique(backgroundcid)) ### number of background compounds
            
	       if(length(pathList[,1])>0)
	         {
               for(i in 1:length(pathList[,1]))
                  {
                     ann<-list(pathwayId=character(),pathwayName="not known",annComponentList=character(),annComponentNumber=0,
                               annBgComponentList=character(),annBgNumber=0,componentNumber=0,bgNumber=0,pvalue=1,fdr=1)
                      
                      ann$pathwayId<-pathList[i,1]  ###pathwayId
                      ann$pathwayName<-pathList[i,2]###pathwayName
                      ann$annBgComponentList<-pathList[i,3]###annBgComponentList
                      ann$bgNumber<-N                      ### bgNumber
                      ann$componentNumber<-length(componentList) ###componentNumber
                      
                      tpathmetabolite<-strsplit(pathList[i,"annBgComponentList"],';')[[1]]####all metabolites in pathway i
                      m1<-length(unique(unlist(tpathmetabolite)))
                      inter<-unique(intersect(as.character(componentList),unlist(tpathmetabolite)))
                      x<-length(unique(inter))
                      
                      ann$annBgNumber<-m1                   ####annBgNumber
                      ann$annComponentNumber<-x             ####annComponentNumber
                      ann$annComponentList<-inter           ####annComponentList
                

                     annList[[i]]<-ann
                  }#####for i in 1:length(pathList[,1])
                  
                  
                  
              ###############################calculate Pvalue###########################
                       annList<-annList[sapply(annList,function(x) x$annComponentNumber>(annlim-1))]
          if(length(annList)>0)
            {
                      
                       annList<-annList[sapply(annList,function(x) x$annBgNumber>(bglim-1))]
                if(length(annList)>0)
                  {    
                    for(i in (1:length(annList)))
                       {
             
                          
                               
                                pvalue<- pWNCHypergeo((as.numeric(annList[[i]]["annComponentNumber"])-1),as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["bgNumber"])-as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["componentNumber"]),1,precision=1E-100,lower.tail=FALSE)
                                annList[[i]]["pvalue"]<-pvalue
                               
                            
                          
                        }####### for(i in (1:length(annList)))
            
                   }########if(length(annList)>0)
            
            }########if(length(annList)>0)
            
            }####if(length(pathList[,1])>0)
         
         
         }######if method=="Hyper"
          
             
             
             
             
             }############# if (pathtype=="KEGG")
           
       
       
      
           
           
           
           
           
           
 else
        {
         
         if((length(pathType)==1)&(length(pathType[which(pathType=="consensusPath")])>0))
           {
         
              pathList<-Getenvir("consensusPath")
            }  
         
         else
         {
            pathListall<-Getenvir("consensusPath")
            pathList<-pathListall[which(as.character(pathListall[,2])%in%pathType),]
          }
         if(method=="MPINet")
          {
           componentList<-unique(intersect(componentList,backgroundcid))
           message(c("your input componentList have ",length(componentList)," components in background"))
           compoundListinnet<-intersect(Veg,componentList)
           message(c("your input componentList have ",length(compoundListinnet)," components in network"))
           
           annList<-list()
           N<-length(unique(backgroundcid)) ###number of background compounds
           
            if(length(pathList[,1])>0)
	         {
               for(i in 1:length(pathList[,1]))
                  {
                     ann<-list(pathwayName="not known",pathsource=character(),annComponentList=character(),annComponentNumber=0,
                               annBgComponentList=character(),annBgNumber=0,componentNumber=0,bgNumber=0,pvalue=1,fdr=1,
                               InWeight=0,weight=1,anncompinNetworkNum=0,anncompinNetworkList=character(),riskcompinNetworkNum=0,riskcompinNetworkList=character())
                      
                      ann$pathsource<-pathList[i,2]  ###pathwayId
                      ann$pathwayName<-pathList[i,1]###pathwayName
                      ann$annBgComponentList<-pathList[i,3]###annBgComponentList
                      ann$bgNumber<-N                      ### bgNumber
                      ann$componentNumber<-length(componentList) ###componentNumber
                      
                      tpathmetabolite<-strsplit(pathList[i,"metabolites"],';')[[1]]####all metabolites in pathway i
                      m1<-length(unique(unlist(tpathmetabolite)))
                      inter<-unique(intersect(as.character(componentList),unlist(tpathmetabolite)))
                      x<-length(unique(inter))
                      
                      ann$annBgNumber<-m1                   ####annBgNumber
                      ann$annComponentNumber<-x             ####annComponentNumber
                      ann$annComponentList<-inter           ####annComponentList
                      
                  #######################################calculate anncompinNetworkNum,anncompinNetworkList,riskcompinNetworkNum,riskcompinNetworkList of pathway i##########################
                      inter_2<-intersect(unique(inter),Veg)########differential metabolites annotate in pathway i and network
                      riskinnet<-intersect(componentList,Veg)#######differential metabolites in network
                      
                      ann$anncompinNetworkNum<-length(inter_2) ###anncompinNetworkNum
                      ann$anncompinNetworkList<-inter_2        ###anncompinNetworkList
                      ann$riskcompinNetworkNum<-length(riskinnet) ###riskcompinNetworkNum
                      ann$riskcompinNetworkList<-riskinnet        ###riskcompinNetworkList
                      
                 #######################calculate pathway InWeigtht################
                      
                                       strinner<-0
                                       
                                       metaboliteinpath<-as.character(unlist(tpathmetabolite))
                                       metapathinnet<-intersect(Veg,unlist(tpathmetabolite))
                                       
                                       if(length(metapathinnet)>0)
                                           {
                                                pathw<-sum(PSS[metapathinnet,"pss"])
                                                strinner<-pathw/length(metaboliteinpath)
                                       
                                               
                                            }
                                       
                                       
                                       
                                       
                                       
                   
                     
                     ann$InWeight<-strinner ###InWeight
                     
                     
                      
                      
                     annList[[i]]<-ann
                  }#####for i in 1:length(pathList[,1])
                  
                  annList<-annList[sapply(annList,function(x) x$annComponentNumber>(annlim-1))]
                  
        if(length(annList)>0)
          {
                  annList<-annList[sapply(annList,function(x) x$annBgNumber>(bglim-1))]
             
          if(length(annList)>0)
            {
             InWeightList<-c()
             for(j in (1:length(annList)))
               {
                   InWeightList<-c(InWeightList,annList[[j]]$InWeight)
               
               }
               meanweight<-mean(InWeightList)
              
               for(j in (1:length(annList)))
               {
                    
                         w<-annList[[j]]$InWeight/meanweight
                         
                           annList[[j]]$weight<-w^weightnum     
                          
                           
                           
                          
                           
                         
               
               }
                  
              ###############################calculate Pvalue###########################
                       
                    
                    for(i in (1:length(annList)))
                       {
             
                          tv<-as.numeric(annList[[i]]["weight"])
                          
                          if(tv>0)
                              {
                               
                                pvalue<- pWNCHypergeo((as.numeric(annList[[i]]["annComponentNumber"])-1),as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["bgNumber"])-as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["componentNumber"]),tv,precision=1E-100,lower.tail=FALSE)
                                annList[[i]]["pvalue"]<-pvalue
                               
                              }
                          
                        }####### for(i in (1:length(annList)))
            
             }########if(length(annList)>0)
      }########if(length(annList)>0)      
            
            }####if(length(pathList[,1])>0)
        }#####if method=="MPINet"   
        
   if(method=="Hyper")
         {
            componentList<-unique(intersect(componentList,backgroundcid))
            message(c("your input componentList have ",length(componentList)," components in background"))
            annList<-list()
            N<-length(unique(backgroundcid)) ###number of background compounds
            
	       if(length(pathList[,1])>0)
	         {
               for(i in 1:length(pathList[,1]))
                  {
                     ann<-list(pathwayName="not known",pathsource=character(),annComponentList=character(),annComponentNumber=0,
                               annBgComponentList=character(),annBgNumber=0,componentNumber=0,bgNumber=0,pvalue=1,fdr=1)
                      
                    ann$pathsource<-pathList[i,2]  ###pathwayId
                      ann$pathwayName<-pathList[i,1]###pathwayName
                      ann$annBgComponentList<-pathList[i,3]###annBgComponentList
                      ann$bgNumber<-N                      ### bgNumber
                      ann$componentNumber<-length(componentList) ###componentNumber
                      
                      tpathmetabolite<-strsplit(pathList[i,"metabolites"],';')[[1]]#### all metabolites in pathway i
                      m1<-length(unique(unlist(tpathmetabolite)))
                      inter<-unique(intersect(as.character(componentList),unlist(tpathmetabolite)))
                      x<-length(unique(inter))
                      
                      ann$annBgNumber<-m1                   ####annBgNumber
                      ann$annComponentNumber<-x             ####annComponentNumber
                      ann$annComponentList<-inter           ####annComponentList
                

                     annList[[i]]<-ann
                  }#####for i in 1:length(pathList[,1])
                  
                  
                  
              ###############################calculate Pvalue###########################
                       annList<-annList[sapply(annList,function(x) x$annComponentNumber>(annlim-1))]
               if(length(annList)>0)
                 {        
                       annList<-annList[sapply(annList,function(x) x$annBgNumber>(bglim-1))]
                  if(length(annList)>0)
                    {
                      for(i in (1:length(annList)))
                       {
             
                          
                               
                                pvalue<- pWNCHypergeo((as.numeric(annList[[i]]["annComponentNumber"])-1),as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["bgNumber"])-as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["componentNumber"]),1,precision=1E-100,lower.tail=FALSE)
                                annList[[i]]["pvalue"]<-pvalue
                               
                            
                          
                        }####### for(i in (1:length(annList)))
            
                  }###########if(length(annList)>0)
             }####if(length(annList)>0)
            }####if(length(pathList[,1])>0)
         
         
         }######if method=="Hyper"
               
               
   }###########else
       
          
###############################calculate fdr#####################################
         if(length(annList)>0)
             {
	                 p_value<-sapply(annList,function(x) return(x$pvalue))
                       #fdrtool.List<-fdrtool(p_value,statistic="pvalue",plot=FALSE,verbose=FALSE)
         	 
                       #print(fdrtool.List$qval)
                       #for(i in seq(annList)){
                       #   annList[[i]]$qvalue<-fdrtool.List$qval[i]
		                   #	annList[[i]]$lfdr<-fdrtool.List$lfdr[i]
                       #}
		              fdr.List<-fdr.est(p_value)
		             for(i in seq(annList))
		                {
		                  annList[[i]]$fdr<-fdr.List[i]
		                }
                      #names(annList)<-sapply(graphList,function(x) x$number)
                 
                 annList<-annList[order(sapply(annList,function(x) x[[order]]),decreasing=decreasing)]   
	           }
           
return(annList)
           
           
       }######function identifypathway
       





#####################################################################
fdr.est<-function(p)
{
    m <- length(ind <- which(!is.na(p)))
    fdr <- rep(NA, length(p))
    stat <- cbind(1:length(p), p, fdr)
    stat[ind, 3] <- unlist(lapply(stat[ind, 2], function(x) {
        c <- length(which(stat[ind, 2] <= x))
        m * x/c
    }))
    stat[ind, ] <- stat[ind, ][order(stat[, 2], decreasing = TRUE), 
        ]
    stat[ind, 3] <- cummin(stat[ind, 3])
    fdr <- stat[order(stat[, 1]), 3]
    fdr[which(fdr > 1)] <- 1
    return(fdr)
}








#####################################################################
printGraph<-function(ann,detail=FALSE,method="MPINet",pathType="KEGG")
   {
	    pathType<-unique(pathType)
	if((length(pathType)==1)&(length(pathType[which(pathType=="KEGG")])>0)){
	    if(method=="MPINet")
	     {
	      if(detail==FALSE)
	       {
	            pathwayId<-sapply(ann,function(x) x$pathwayId)
              pathwayName<-sapply(ann,function(x) x$pathwayName)
              annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
              annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
              pvalue<-sapply(ann,function(x) x$pvalue)
              fdr<-sapply(ann,function(x) x$fdr)
              #InWeight<-sapply(ann,function(x) x$InWeight)
              #OutWeight<-sapply(ann,function(x) x$OutWeight)
              #AllWeight<-sapply(ann,function(x) x$AllWeight)
              #anncompinNetworkNum<-sapply(ann,function(x) x$anncompinNetworkNum)
              #riskcompinNetworkNum<-sapply(ann,function(x) x$riskcompinNetworkNum)
              annComponentinNetRatio<-sapply(ann,function(x) paste(x$anncompinNetworkNum,x$riskcompinNetworkNum,sep="/"))
              weight<-sapply(ann,function(x) x$weight)
     
      
                  #qvalue<-sapply(ann,function(x) x$qvalue)
	                #lfdr<-sapply(ann,function(x) x$lfdr)
	                # fdr<-sapply(ann,function(x) x$fdr)
                  #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annComponentRatio,
                  #                       annBgRatio,pvalue,qvalue,lfdr))
              ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annComponentRatio=annComponentRatio,
	            annBgRatio=annBgRatio,weight=weight,pvalue=pvalue,fdr=fdr,stringsAsFactors=FALSE,
	            annComponentinNetRatio=annComponentinNetRatio)							 
	        }
	     else
	       {	 
             pathwayId<-sapply(ann,function(x) x$pathwayId)	  
	           pathwayName<-sapply(ann,function(x) x$pathwayName)
	           annComponentList<-sapply(ann, function(x){ paste(x$annComponentList,collapse=";") })
             annBgComponentList<-sapply(ann, function(x) x$annBgComponentList)#######
	           annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
             annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
             pvalue<-sapply(ann,function(x) x$pvalue)
      
             #InWeight<-sapply(ann,function(x) x$InWeight)
             #OutWeight<-sapply(ann,function(x) x$OutWeight)
             #AllWeight<-sapply(ann,function(x) x$AllWeight)
            # anncompinNetworkNum<-sapply(ann,function(x) x$anncompinNetworkNum)
             anncompinNetworkList<-sapply(ann, function(x){ paste(x$anncompinNetworkList,collapse=";") })
             #riskcompinNetworkNum<-sapply(ann,function(x) x$riskcompinNetworkNum)
             annComponentinNetRatio<-sapply(ann,function(x) paste(x$anncompinNetworkNum,x$riskcompinNetworkNum,sep="/"))
             riskcompinNetworkList<-sapply(ann, function(x){ paste(x$riskcompinNetworkList,collapse=";") })
             weight<-sapply(ann,function(x) x$weight)
                 #qvalue<-sapply(ann,function(x) x$qvalue)
	               #lfdr<-sapply(ann,function(x) x$lfdr)
	           fdr<-sapply(ann,function(x) x$fdr)
                 #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annComponentRatio,
                 #                       annBgRatio,pvalue,qvalue,lfdr,annComponentList,annBgComponentList))
             ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annComponentRatio=annComponentRatio,
	           annBgRatio=annBgRatio,weight=weight,pvalue=pvalue,fdr=fdr,annComponentList=annComponentList,
	           annBgComponentList=annBgComponentList,
	           annComponentinNetRatio=annComponentinNetRatio,anncompinNetworkList=anncompinNetworkList,riskcompinNetworkList=riskcompinNetworkList,stringsAsFactors=FALSE)								 
	      }
      return(ann.data.frame)
      }############if method="MPINet"
    if(method=="Hyper")
      {
        if(detail==FALSE)
	       {
	            pathwayId<-sapply(ann,function(x) x$pathwayId)
              pathwayName<-sapply(ann,function(x) x$pathwayName)
              annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
              annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
              pvalue<-sapply(ann,function(x) x$pvalue)
              fdr<-sapply(ann,function(x) x$fdr)
              #InWeight<-sapply(ann,function(x) x$InWeight)
              #OutWeight<-sapply(ann,function(x) x$OutWeight)
              #AllWeight<-sapply(ann,function(x) x$AllWeight)
              #anncompinNetworkNum<-sapply(ann,function(x) x$anncompinNetworkNum)
              #riskcompinNetworkNum<-sapply(ann,function(x) x$riskcompinNetworkNum)
             # annComponentinNetRatio<-sapply(ann,function(x) paste(x$anncompinNetworkNum,x$riskcompinNetworkNum,sep="/"))
             # Tw<-sapply(ann,function(x) x$Tw)
     
      
                  #qvalue<-sapply(ann,function(x) x$qvalue)
	                #lfdr<-sapply(ann,function(x) x$lfdr)
	                
                  #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annComponentRatio,
                  #                       annBgRatio,pvalue,qvalue,lfdr))
              ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annComponentRatio=annComponentRatio,
	            annBgRatio=annBgRatio,pvalue=pvalue,fdr=fdr,stringsAsFactors=FALSE)							 
	        }
	     else
	       {	 
             pathwayId<-sapply(ann,function(x) x$pathwayId)	  
	           pathwayName<-sapply(ann,function(x) x$pathwayName)
	           annComponentList<-sapply(ann, function(x){ paste(x$annComponentList,collapse=";") })
             annBgComponentList<-sapply(ann, function(x) x$annBgComponentList)#######
	           annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
             annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
             pvalue<-sapply(ann,function(x) x$pvalue)
      
             #InWeight<-sapply(ann,function(x) x$InWeight)
             #OutWeight<-sapply(ann,function(x) x$OutWeight)
             #AllWeight<-sapply(ann,function(x) x$AllWeight)
            # anncompinNetworkNum<-sapply(ann,function(x) x$anncompinNetworkNum)
             #anncompinNetworkList<-sapply(ann, function(x){ paste(x$anncompinNetworkList,collapse=";") })
             #riskcompinNetworkNum<-sapply(ann,function(x) x$riskcompinNetworkNum)
             #annComponentinNetRatio<-sapply(ann,function(x) paste(x$anncompinNetworkNum,x$riskcompinNetworkNum,sep="/"))
             #riskcompinNetworkList<-sapply(ann, function(x){ paste(x$riskcompinNetworkList,collapse=";") })
             #Tw<-sapply(ann,function(x) x$Tw)
                 #qvalue<-sapply(ann,function(x) x$qvalue)
	               #lfdr<-sapply(ann,function(x) x$lfdr)
	           fdr<-sapply(ann,function(x) x$fdr)
                 #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annComponentRatio,
                 #                       annBgRatio,pvalue,qvalue,lfdr,annComponentList,annBgComponentList))
             ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annComponentRatio=annComponentRatio,
	           annBgRatio=annBgRatio,pvalue=pvalue,fdr=fdr,annComponentList=annComponentList,
	           annBgComponentList=annBgComponentList,stringsAsFactors=FALSE)								 
	      }
      return(ann.data.frame)
      }  ########if(method=="Hyper")
 }######if (pathType=="KEGG")
 
	    
else
	  {
	    if(method=="MPINet")
	     {
	      if(detail==FALSE)
	       {
	            pathsource<-sapply(ann,function(x) x$pathsource)
              pathwayName<-sapply(ann,function(x) x$pathwayName)
              annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
              annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
              pvalue<-sapply(ann,function(x) x$pvalue)
      
              #InWeight<-sapply(ann,function(x) x$InWeight)
              #OutWeight<-sapply(ann,function(x) x$OutWeight)
              #AllWeight<-sapply(ann,function(x) x$AllWeight)
              #anncompinNetworkNum<-sapply(ann,function(x) x$anncompinNetworkNum)
              #riskcompinNetworkNum<-sapply(ann,function(x) x$riskcompinNetworkNum)
              annComponentinNetRatio<-sapply(ann,function(x) paste(x$anncompinNetworkNum,x$riskcompinNetworkNum,sep="/"))
              weight<-sapply(ann,function(x) x$weight)
              fdr<-sapply(ann,function(x) x$fdr)
      
                  #qvalue<-sapply(ann,function(x) x$qvalue)
	                #lfdr<-sapply(ann,function(x) x$lfdr)
	                # fdr<-sapply(ann,function(x) x$fdr)
                  #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annComponentRatio,
                  #                       annBgRatio,pvalue,qvalue,lfdr))
              ann.data.frame<-data.frame(pathwayName=pathwayName,pathsource=pathsource,annComponentRatio=annComponentRatio,
	            annBgRatio=annBgRatio,weight=weight,pvalue=pvalue,fdr=fdr,stringsAsFactors=FALSE,
	            annComponentinNetRatio=annComponentinNetRatio)							 
	        }
	     else
	       {	 
             pathsource<-sapply(ann,function(x) x$pathsource)	  
	           pathwayName<-sapply(ann,function(x) x$pathwayName)
	           annComponentList<-sapply(ann, function(x){ paste(x$annComponentList,collapse=";") })
             annBgComponentList<-sapply(ann, function(x) x$annBgComponentList)#######
	           annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
             annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
             pvalue<-sapply(ann,function(x) x$pvalue)
      
             #InWeight<-sapply(ann,function(x) x$InWeight)
             #OutWeight<-sapply(ann,function(x) x$OutWeight)
             #AllWeight<-sapply(ann,function(x) x$AllWeight)
            # anncompinNetworkNum<-sapply(ann,function(x) x$anncompinNetworkNum)
             anncompinNetworkList<-sapply(ann, function(x){ paste(x$anncompinNetworkList,collapse=";") })
             #riskcompinNetworkNum<-sapply(ann,function(x) x$riskcompinNetworkNum)
             annComponentinNetRatio<-sapply(ann,function(x) paste(x$anncompinNetworkNum,x$riskcompinNetworkNum,sep="/"))
             riskcompinNetworkList<-sapply(ann, function(x){ paste(x$riskcompinNetworkList,collapse=";") })
             weight<-sapply(ann,function(x) x$weight)
                 #qvalue<-sapply(ann,function(x) x$qvalue)
	               #lfdr<-sapply(ann,function(x) x$lfdr)
	           fdr<-sapply(ann,function(x) x$fdr)
                 #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annComponentRatio,
                 #                       annBgRatio,pvalue,qvalue,lfdr,annComponentList,annBgComponentList))
             ann.data.frame<-data.frame(pathwayName=pathwayName,pathsource=pathsource,annComponentRatio=annComponentRatio,
	           annBgRatio=annBgRatio,weight=weight,pvalue=pvalue,fdr=fdr,annComponentList=annComponentList,
	           annBgComponentList=annBgComponentList,annComponentinNetRatio=annComponentinNetRatio,
	           anncompinNetworkList=anncompinNetworkList,riskcompinNetworkList=riskcompinNetworkList,stringsAsFactors=FALSE)								 
	      }
      return(ann.data.frame)
      }############if method="MPINet"
    if(method=="Hyper")
      {
          if(detail==FALSE)
	       {
	            pathsource<-sapply(ann,function(x) x$pathsource)
              pathwayName<-sapply(ann,function(x) x$pathwayName)
              annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
              annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
              pvalue<-sapply(ann,function(x) x$pvalue)
              fdr<-sapply(ann,function(x) x$fdr)
             
              ann.data.frame<-data.frame(pathwayName=pathwayName,pathsource=pathsource,annComponentRatio=annComponentRatio,
	            annBgRatio=annBgRatio,pvalue=pvalue,fdr=fdr,stringsAsFactors=FALSE)							 
	        }
	     else
	       {	 
             pathsource<-sapply(ann,function(x) x$pathsource)	  
	           pathwayName<-sapply(ann,function(x) x$pathwayName)
	           annComponentList<-sapply(ann, function(x){ paste(x$annComponentList,collapse=";") })
             annBgComponentList<-sapply(ann, function(x) x$annBgComponentList)#######
	           annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
             annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
             pvalue<-sapply(ann,function(x) x$pvalue)
      
             #InWeight<-sapply(ann,function(x) x$InWeight)
             #OutWeight<-sapply(ann,function(x) x$OutWeight)
             #AllWeight<-sapply(ann,function(x) x$AllWeight)
            # anncompinNetworkNum<-sapply(ann,function(x) x$anncompinNetworkNum)
             #anncompinNetworkList<-sapply(ann, function(x){ paste(x$anncompinNetworkList,collapse=";") })
             #riskcompinNetworkNum<-sapply(ann,function(x) x$riskcompinNetworkNum)
             #annComponentinNetRatio<-sapply(ann,function(x) paste(x$anncompinNetworkNum,x$riskcompinNetworkNum,sep="/"))
             #riskcompinNetworkList<-sapply(ann, function(x){ paste(x$riskcompinNetworkList,collapse=";") })
             #Tw<-sapply(ann,function(x) x$Tw)
                 #qvalue<-sapply(ann,function(x) x$qvalue)
	               #lfdr<-sapply(ann,function(x) x$lfdr)
	           fdr<-sapply(ann,function(x) x$fdr)
                 #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annComponentRatio,
                 #                       annBgRatio,pvalue,qvalue,lfdr,annComponentList,annBgComponentList))
             ann.data.frame<-data.frame(pathwayName=pathwayName,pathsource=pathsource,annComponentRatio=annComponentRatio,
	           annBgRatio=annBgRatio,pvalue=pvalue,fdr=fdr,annComponentList=annComponentList,
	           annBgComponentList=annBgComponentList,stringsAsFactors=FALSE)								 
	      }
      return(ann.data.frame)
      }  ########if(method=="Hyper")
	      
	      
	      }###########else
	   
}##########printGraph

