rsmlToDART <- function(rsml.path, final.date, connect){
  
  # Create the plant object
  rsml <- xmlToList(xmlParse(rsml.path))
  plants <- rsml$scene
  
  # Create LIE and RAC files for each root system 
  
  n <- 0 # Initialise number of LIE/RAC files (as 1 RSML file can contain more than 1 first-order root)
  lie.all<-list()
  rac.all<-list()
  tps.all<-list()
  
  for (r0 in plants){# For each plant
    
    for (r1 in r0){ # For each first-order root
      
      r <- 0 # Initialise the number of roots consisting a root system
      
      if (class(r1)=="list"){# For each first order root, create LIE and RAC files
        
        n<-n+1
        currentMother<-0
        lie<-data.frame(Num=0,Date=0, Bran=0, Apic=0, Prec=0, Suiv=0, X=0, Y=0, dist=0)
        rac<-data.frame(Root=0, Mother=0, Ord=0, DBase=0, DApp=0, Lengths=0)
        tps<-data.frame(Num=1, Date=final.date)
        ns<-r1$geometry$polyline
        length1<-length(ns)
        r<-r+1
        
        #c(Num, Date, Bran, Apic, Prec, Suiv)
        lie[1:length1,1:6]<-c(c(1:length1),rep(1, length1),rep(0, length1), rep(0, length1), 0:(length1-1), 2:(length1+1))
        
        #c(X,Y)
        lie[1:length1,7]<-sapply(ns, xnodes)
        lie[1:length1,8]<-sapply(ns, ynodes)
        
        #c(dist)
        lie[1:length1, 9]<-c(0, cumsum(sqrt((diff(lie$X[1:length1]))^2+(diff(lie$Y[1:length1]))^2)))
        
        #For the first point of the first order root
        lie[1,c(2:3)]<-c(0,1)
        start1<-as.numeric(lie$Num[1])
        
        #For the last point of the first order root
        lie[length1,c(4,6)]<-c(1,0)
        stop1<-as.numeric(lie$Num[length1])
        
        # Fill RAC file for the first order root
        #c(Root, Mother, Ord, DBase, DApp, Length)
        cumulDist<-sum(sqrt((diff(lie$X[1:length1]))^2+(diff(lie$Y[1:length1]))^2))
        rac[r,1:6]<-c(0,-1,1,0,0,cumulDist)
        
        #----------------------------------------------------------------------------------------------------
        
        # If there are lateral roots
        
        if ("root" %in% names(r1)){
          
          for (r2 in r1){# For each 2-order root
            
            if ("geometry" %in% names(r2)){
              
              r<-r+1
              ns <- r2$geometry$polyline
              length2<-length(ns)
              lie.lines<-nrow(lie)
              
              #c(Num, Date, Bran, Apic, Prec, Suiv)
              lie[(lie.lines+1):(lie.lines+length2),1:6]<-c((lie.lines+1):(lie.lines+length2), rep(1, length2), rep(0, length2), rep(0, length2), lie.lines:(lie.lines+length2-1), (lie.lines+2):(lie.lines+length2+1))
              
              #c(X,Y)
              lie[(lie.lines+1):(lie.lines+length2),7]<-sapply(ns, xnodes)
              lie[(lie.lines+1):(lie.lines+length2),8]<-sapply(ns, ynodes)
              
              #c(dist)
              lie[(lie.lines+1):(lie.lines+length2), 9]<-c(0, cumsum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length2)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length2)]))^2)))
              
              # Search the closest point on the mother root (calculate DBase)
              dist1<-sqrt((lie$X[start1:stop1]-lie$X[lie.lines+1])^2+(lie$Y[start1:stop1]-lie$Y[lie.lines+1])^2)
              index<-match(min(dist1), dist1)
              dbase<-lie$dist[start1+index-1]
              prec<-lie$Num[start1+index-1]
              
              # Change Prec and Bran values for the first point of a lateral root
              if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1} else {lie[lie.lines+1,5]<-prec}
              lie[lie.lines+1,3]<-1
              start2<-as.numeric(lie$Num[lie.lines+1])
              
              # Change Suiv and Apic values for the last point of a lateral root
              lie[lie.lines+length2,c(4,6)]<-c(1, 0)
              stop2<-as.numeric(lie$Num[lie.lines+length2])
              
              # Fill RAC file for the 2-order root
              
              if (connect==TRUE){cumulDist<-sqrt((lie$X[prec]-lie$X[lie.lines+1])^2+(lie$Y[prec]-lie$Y[lie.lines+1])^2) + sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length2)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length2)]))^2))}
              else {cumulDist<-sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length2)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length2)]))^2))}
              rac[r, 1:6]<-c(max(rac$Root)+1, currentMother, 2, dbase, final.date/2, cumulDist)
              
              currentMother2<-rac$Root[r]
              
              #----------------------------------------------------------------------------------------------------            
              
              if ("root" %in% names(r2)){
                
                for (r3 in r2){# For each 3-order root
                  
                  if ("geometry" %in% names(r3)){
                    
                    r<-r+1
                    ns <- r3$geometry$polyline
                    length3<-length(ns)
                    lie.lines<-nrow(lie)
                    
                    #c(Num, Date, Bran, Apic, Prec, Suiv)
                    lie[(lie.lines+1):(lie.lines+length3), 1:6]<-c((lie.lines+1):(lie.lines+length3), rep(1, length3), rep(0, length3), rep(0, length3), lie.lines:(lie.lines+length3-1), (lie.lines+2):(lie.lines+length3+1))
                    
                    #c(X,Y)
                    lie[(lie.lines+1):(lie.lines+length3),7]<-sapply(ns, xnodes)
                    lie[(lie.lines+1):(lie.lines+length3),8]<-sapply(ns, ynodes)
                    
                    #c(dist)
                    lie[(lie.lines+1):(lie.lines+length3), 9]<-c(0, cumsum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length3)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length3)]))^2)))
                    
                    # Search the closest point on the mother root (calculate DBase)
                    dist1<-sqrt((lie$X[start2:stop2]-lie$X[lie.lines+1])^2+(lie$Y[start2:stop2]-lie$Y[lie.lines+1])^2)
                    index<-match(min(dist1), dist1)
                    dbase<-lie$dist[start2+index-1]
                    prec<-lie$Num[start2+index-1]
                    
                    # Change Prec and Bran values for the first point of a lateral root
                    if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1} else {lie[lie.lines+1,5]<-prec}
                    lie[lie.lines+1,3]<-1
                    start3<-as.numeric(lie$Num[lie.lines+1])
                    
                    # Change Suiv and Apic values for the last point of a lateral root
                    lie[lie.lines+length3,c(4,6)]<-c(1, 0)
                    stop3<-as.numeric(lie$Num[lie.lines+length3])
                    
                    # Fill RAC file for the 3-order root
                    
                    if (connect==TRUE){cumulDist<-sqrt((lie$X[prec]-lie$X[lie.lines+1])^2+(lie$Y[prec]-lie$Y[lie.lines+1])^2) + sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length3)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length3)]))^2))}
                    else {cumulDist<-sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length3)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length3)]))^2))}
                    rac[r, 1:6]<-c(max(rac$Root)+1, currentMother2, 3, dbase, final.date/2, cumulDist)
                    
                    currentMother3<-rac$Root[r]
                    
                    #----------------------------------------------------------------------------------------------------
                    
                    if ("root" %in% names(r3)){
                      
                      for (r4 in r3){# For each 4-order root
                        
                        if ("geometry" %in% names(r4)){
                          
                          r<-r+1
                          ns <- r4$geometry$polyline
                          length4<-length(ns)
                          lie.lines<-nrow(lie)
                          
                          #c(Num, Date, Bran, Apic, Prec, Suiv)
                          lie[(lie.lines+1):(lie.lines+length4), 1:6]<-c((lie.lines+1):(lie.lines+length4), rep(1, length4), rep(0, length4), rep(0, length4), lie.lines:(lie.lines+length4-1), (lie.lines+2):(lie.lines+length4+1))
                          
                          #c(X,Y)
                          lie[(lie.lines+1):(lie.lines+length4),7]<-sapply(ns, xnodes)
                          lie[(lie.lines+1):(lie.lines+length4),8]<-sapply(ns, ynodes)
                          
                          #c(dist)
                          lie[(lie.lines+1):(lie.lines+length4), 9]<-c(0, cumsum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length4)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length4)]))^2)))
                          
                          # Search the closest point on the mother root (calculate DBase)
                          dist1<-sqrt((lie$X[start3:stop3]-lie$X[lie.lines+1])^2+(lie$Y[start3:stop3]-lie$Y[lie.lines+1])^2)
                          index<-match(min(dist1), dist1)
                          dbase<-lie$dist[start3+index-1]
                          prec<-lie$Num[start3+index-1]
                          
                          # Change Prec and Bran values for the first point of a lateral root
                          if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1} else {lie[lie.lines+1,5]<-prec}
                          lie[lie.lines+1,3]<-1
                          start4<-as.numeric(lie$Num[lie.lines+1])
                          
                          # Change Suiv and Apic values for the last point of a lateral root
                          lie[lie.lines+length4,c(4,6)]<-c(1, 0)
                          stop4<-as.numeric(lie$Num[lie.lines+length4])
                          
                          # Fill RAC file for the 4-order root
                          
                          if (connect==TRUE){cumulDist<-sqrt((lie$X[prec]-lie$X[lie.lines+1])^2+(lie$Y[prec]-lie$Y[lie.lines+1])^2) + sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length4)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length4)]))^2))}
                          else {cumulDist<-sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length4)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length4)]))^2))}
                          rac[r, 1:6]<-c(max(rac$Root)+1, currentMother3, 4, dbase, final.date/2, cumulDist)
                          
                          currentMother4<-rac$Root[r]
                          
                          #----------------------------------------------------------------------------------------------------
                          
                          if ("root" %in% names(r4)){
                            
                            for (r5 in r4){# For each 5-order root
                              
                              if ("geometry" %in% names(r5)){
                                
                                r<-r+1
                                ns <- r5$geometry$polyline
                                length5<-length(ns)
                                lie.lines<-nrow(lie)
                                
                                #c(Num, Date, Bran, Apic, Prec, Suiv)
                                lie[(lie.lines+1):(lie.lines+length5), 1:6]<-c((lie.lines+1):(lie.lines+length5), rep(1, length5), rep(0, length5), rep(0, length5), lie.lines:(lie.lines+length5-1), (lie.lines+2):(lie.lines+length5+1))
                                
                                #c(X,Y)
                                lie[(lie.lines+1):(lie.lines+length5),7]<-sapply(ns, xnodes)
                                lie[(lie.lines+1):(lie.lines+length5),8]<-sapply(ns, ynodes)
                                
                                #c(dist)
                                lie[(lie.lines+1):(lie.lines+length5), 9]<-c(0, cumsum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length5)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length5)]))^2)))
                                
                                # Search the closest point on the mother root (calculate DBase)
                                dist1<-sqrt((lie$X[start4:stop4]-lie$X[lie.lines+1])^2+(lie$Y[start4:stop4]-lie$Y[lie.lines+1])^2)
                                index<-match(min(dist1), dist1)
                                dbase<-lie$dist[start4+index-1]
                                prec<-lie$Num[start4+index-1]
                                
                                # Change Prec and Bran values for the first point of a lateral root
                                if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1} else {lie[lie.lines+1,5]<-prec}
                                lie[lie.lines+1,3]<-1
                                
                                # Change Suiv and Apic values for the last point of a lateral root
                                lie[lie.lines+length5,c(4,6)]<-c(1, 0)
                                
                                # Fill RAC file for the 5-order root
                                if (connect==TRUE){cumulDist<-sqrt((lie$X[prec]-lie$X[lie.lines+1])^2+(lie$Y[prec]-lie$Y[lie.lines+1])^2) + sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length5)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length5)]))^2))}
                                else {cumulDist<-sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length5)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length5)]))^2))}
                                rac[r, 1:6]<-c(max(rac$Root)+1, currentMother4, 4, dbase, final.date/2, cumulDist)
                                
                              } 
                            }
                          } 
                        } 
                      }
                    } 
                  } 
                }
              } 
            } 
          }
        }
      }
      lie$Bran[lie$Bran==0]<-"false"
      lie$Bran[lie$Bran==1]<-"true"
      lie$Apic[lie$Apic==0]<-"false"
      lie$Apic[lie$Apic==1]<-"true"
      lie.all[[n]]<-lie
      rac.all[[n]]<-rac
      tps.all[[n]]<-tps
    }
  }
  result<-list(lie=lie.all, rac=rac.all, tps=tps.all)
  return(result)
}
