

################################################################
#______________    Function  ConsensusPartition    _____________
#               Generates the consensus partition
################################################################
# Computation of consensus of partitions based on Adjusted Rand
# Source : Krieger and Green J. of Classification, 16:63-89, 1999
################################################################

ConsensusPartition<-function(Part,ngroups=0,type="cutree",optim=FALSE,maxiter=100,plotDendrogram=FALSE,verbose=FALSE){
  
  if (!class(Part)=="SortingPartition"){
    return("The argument is not an object of class SortingPartition")
  } else  {
    if (type=="fusion"){
      cat("\nFusion algorithm. May be time consuming.\n\n")
    }
    nstim<-Part@nstimuli
    nsubjects<-Part@nsubjects
    Parti<-Part@Partition
    Labels<-Part@LabStim
    


    MatPart<-getPartition(Part)
    npart<-ncol(MatPart)
    
    if (type=="cutree"){
      
      ListDissimil<-Dissimil(Part)
      MatDissim<-apply(simplify2array(ListDissimil),c(1,2),'sum')
      #Hierarchical clustering of stimuli
      hres<-hclust(as.dist(MatDissim), method = "ward.D2")
      
      if (plotDendrogram==TRUE){
        plot(hres,labels=Labels,hang=-1)
      }
      
    }#end if cutree

    if (type=="fusion"){
      
      # fusion algorithm  
      if (ngroups!=0){
        ngrclass<-nstim-ngroups+1
      } else {
        ngrclass<-nstim-1
      } 
      
      tabsubjopt<-matrix(0,nstim,ngrclass-1)
      tabcritopt<-matrix(0,1,ngrclass-1)
      
      #initialisation with nstim groups
      Cons<-1:nstim
      #fusion of two groups
      for (iter in 1:(ngrclass-1)){
        #from nclasses to nclasses-1
        nclasses<-nstim-iter+1
        maxcrit<-0
        for (gr1 in 1:(nclasses-1)){
          for (gr2 in ((gr1+1):(nclasses))){
            #test for fusion of gr1 and gr2
            prov<-Cons
            prov[prov==gr2]<-gr1
            prov[prov>gr2]=prov[prov>gr2]-1
            
            # computation of criterium
            crit<-0
            for (p in 1:npart){
              crit<-crit+RandAdjusted(MatPart[,p],prov)
            }
            crit<-crit/npart
            
            if (crit>maxcrit){
              maxcrit<-crit
              maxgr1<-gr1
              maxgr2<-gr2
            }
          }
        }
        #fusion of maxgr1 and maxgr2
        Cons[Cons==maxgr2]<-maxgr1
        Cons[Cons>maxgr2]<-Cons[Cons>maxgr2]-1
        tabsubjopt[ ,ngrclass-iter]<-Cons
        tabcritopt[ngrclass-iter]<-maxcrit
        
      }# end iter
      names(tabcritopt)<-(nstim-ngrclass+1):(nstim-1)
      colnames(tabsubjopt)<-(nstim-ngrclass+1):(nstim-1)
    }#end if fusion    
    
    if (type=="medoid"){    
    #searching for the closest partition to the others
      maxcrit<-0
      pmax<-0
      for (p in 1:npart){
        crit<-0
        for (p1 in 1:npart){
          if (p!=p1){  
            crit<-crit+RandAdjusted(MatPart[,p],MatPart[,p1])
          }  
        }
        crit<-crit/(npart-1)
        if (crit>maxcrit){
          maxcrit<-crit
          pmax<-p
        }
      }
      consensus<-MatPart[,pmax]
      if (verbose==TRUE){
        cat("Medoid partition\n")
        cat("The consensus is the partition of subject ",pmax, " whose name is ",Part@LabSubj[pmax],".\n\n",sep="")
      }
    }  
      
      
    #searching for optimum  
      
    if (type!="medoid"){    
    if (ngroups!=0){
      if (type=="cutree"){
        consensus<-cutree(hres, k = ngroups)
        names(consensus)<-Labels
        # computation of criterium
        maxcrit<-0
        for (p in 1:npart){
          maxcrit<-maxcrit+RandAdjusted(MatPart[,p],consensus)
        }   
        maxcrit<-maxcrit/npart
      } else {  #type="fusion"
        consensus<-tabsubjopt[,1]
        names(consensus)<-Labels
        maxcrit<-tabcritopt[1]
      }
      
      if (verbose==TRUE){
        cat("Criterion: ",maxcrit,"\n",sep="")
        cat("\nConsensus:\n")
        print(consensus)
      }
      
    } else {
      if (type=="cutree"){
        
        ngrclass<-nstim-1
        tabsubjopt<-matrix(0,nstim,nstim-2)
        tabcritopt<-matrix(0,1,nstim-2)
        colnames(tabsubjopt)<-2:(nstim-1)
        names(tabcritopt)<-2:(nstim-1)
        for (groups in (2:(nstim-1))){
          Cons<-cutree(hres, k = groups)
          
          # computation of criterium
          crit<-0
          for (p in 1:npart){
            crit<-crit+RandAdjusted(MatPart[,p],Cons)
          }
          crit<-crit/npart
          tabsubjopt[,groups-1]<-Cons
          tabcritopt[1,groups-1]<-crit
        }
      }
      
      groupmax<-names(which.max(tabcritopt))
      maxcrit<-max(tabcritopt)
      consensus<-tabsubjopt[,colnames(tabsubjopt)==groupmax]
      names(consensus)<-Labels  
      
      if (verbose==TRUE){
        cat("Table of index as a function of the number of groups:\n")
        print(tabcritopt)
        cat("\nTable of consensus as a function of the number of groups:\n")
        print(tabsubjopt)
        cat("\nOptimal consensus with ",groupmax, " groups and criterion ",maxcrit," :\n\n",sep="")
        print(consensus)
      }
      
    }
    }
    
    if (type!="medoid" & optim==TRUE){
      converge<-FALSE
      iter<-0
      obj<-0
      nclass<-sum(unique(consensus)!=0)
      while  ( iter<maxiter & converge==FALSE ) {
        change<-FALSE
        #permutation of stimuli
        permut<-sample(1:nstim)
        Singleton<-FALSE
        for ( i in permut){
          prov<-consensus
          if (sum(prov==prov[permut[i]])>1){    # if i is not a singleton
            maxcrit<-0
            cmax<-1
            for ( c in (1:nclass)){
              prov[permut[i]]<-c                # stimulus i is now belonging to the class c
              crit<-0                   # computation of criterium
              for (p in (1:npart)){
                crit<-crit+RandAdjusted(MatPart[,p],prov)
              }
              crit<-crit/npart
              if (crit>maxcrit){         # better criterium is found
                maxcrit<-crit
                cmax<-c
              }
            } # end c
            if (consensus[permut[i]]!=cmax){ 
              consensus[permut[i]]<-cmax   # i is now belonging to class cmax
              change<-TRUE
            }
          } # end if i singleton
          else {
            Singleton<-TRUE
          }  
        } # end i in permut  
        converge<- abs(maxcrit-obj)<1e-06       # test for convergence : criterum is not improved
        if (converge==TRUE & Singleton==TRUE){
          #in case of singleton trying to exchange 2 stimuli
          for ( i in permut){
            if (sum(consensus==consensus[permut[i]])==1){    # if i is a singleton
              for (j in permut){
                #test for exchanging i and j
                prov<-consensus
                gri<-prov[permut[i]]
                prov[permut[i]]<-prov[permut[j]]
                prov[permut[j]]<-gri
                
                crit<-0                   # computation of criterium
                for (p in (1:npart)){
                  crit<-crit+RandAdjusted(MatPart[,p],prov)
                }
                crit<-crit/npart
                if (crit>maxcrit){
                  change=TRUE
                  maxcrit<-crit
                  candi<-i
                  candj<-j
                }
              }
            }#singleton
          }#end i
          if (change==TRUE){
            Converge<-FALSE
            #exchange i and j
            grj<-consensus[permut[candj]]
            consensus[permut[candj]]<-consensus[permut[candi]]
            consensus[permut[candi]]<-grj
          }
        }# end if Singleton==TRUE
        obj<-maxcrit
        iter<-iter+1
      }   # end while
      
      cat("Final consensus with ",length(unique(consensus)), " groups and criterion ",maxcrit," :\n",sep="")    
      print(consensus)
      
    }# end if optim==TRUE
    
    
    
    return(list(Consensus=consensus,Crit=maxcrit))
    
  }
} 
