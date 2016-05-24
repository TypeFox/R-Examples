obtain.ix0<-function(x, y, s1, s2, family, nsis, iter, varISIS, perm, q, greedy, greedy.size, iterind){
if(iter == FALSE){
   margcoef = margcoef(x, y, family=family, null.model=FALSE, iterind=iterind)
   rankcoef = sort(margcoef, decreasing=TRUE, index.return=TRUE)      
   ix0 = rankcoef$ix[1:nsis]
}
else{  
   if(varISIS=="vanilla"){
      margcoef = margcoef(x, y, family=family, null.model=FALSE, iterind=iterind)
      rankcoef = sort(margcoef, decreasing=TRUE, index.return=TRUE)            
      if(perm == FALSE) 
         ix0 = rankcoef$ix[1:floor((2/3)*nsis)]                        
      else{
         repeat{
            randcoef = margcoef(x, y, family=family, null.model=TRUE, iterind=iterind)
            if(length(which(margcoef >= quantile(randcoef,q))) > 0)
               break                     
         }             
         if(greedy == FALSE){          
            if(length(which(margcoef >= quantile(randcoef,q))) >= 2){              
               length1 = length(which(margcoef >= quantile(randcoef,q)))       
               above.thresh = rankcoef$ix[1:length1]     
               ix0 = rankcoef$ix[1:floor((2/3)*nsis)]                 
               ix0 = sort(intersect(ix0,above.thresh))                     
            }else 
               ix0 = rankcoef$ix[1:2]   
         }else{
            if(greedy.size == 1)
               ix0 = rankcoef$ix[1:2]
            else
               ix0 = rankcoef$ix[1:greedy.size]
         }                        
      }
   }else{    
       margcoef1 = margcoef(x[s1,], y[s1], family=family, null.model=FALSE, iterind=iterind)
       margcoef2 = margcoef(x[s2,], y[s2], family=family, null.model=FALSE, iterind=iterind)
       rankcoef1 = sort(margcoef1, decreasing=TRUE, index.return=TRUE)
       rankcoef2 = sort(margcoef2, decreasing=TRUE, index.return=TRUE)        
       if(perm == FALSE){
          if(varISIS=="aggr"){
             ix01 = rankcoef1$ix[1:floor((2/3)*nsis)]
             ix02 = rankcoef2$ix[1:floor((2/3)*nsis)]
             ix0 = sort(intersect(ix01,ix02))
             if(length(ix0) <= 1) ix0 = int.size.k(rankcoef1$ix,rankcoef2$ix,2)
          }
          if(varISIS=="cons"){
             iensure = intensure(floor((2/3)*nsis),l1=rankcoef1$ix,l2=rankcoef2$ix,k=floor((2/3)*nsis))       
             ix01 = rankcoef1$ix[1:iensure]
             ix02 = rankcoef2$ix[1:iensure]
             ix0 = sort(intersect(ix01,ix02))  
          }        
       }else{
           repeat{
              randcoef1 = margcoef(x[s1,], y[s1], family=family, null.model=TRUE, iterind=iterind)
              randcoef2 = margcoef(x[s2,], y[s2], family=family, null.model=TRUE, iterind=iterind)
              if(length(which(margcoef1 >= quantile(randcoef1,q))) > 0 && length(which(margcoef2 >= quantile(randcoef2,q))) > 0)
                 break                     
           }                  
           if(greedy == FALSE){                     
              length1 = length(which(margcoef1 >= quantile(randcoef1,q)))    
              length2 = length(which(margcoef2 >= quantile(randcoef2,q)))   
              above.thresh.1 = rankcoef1$ix[1:length1]     
              above.thresh.2 = rankcoef2$ix[1:length2]
              ix01 = rankcoef1$ix[1:floor((2/3)*nsis)]                 
              ix02 = rankcoef2$ix[1:floor((2/3)*nsis)]  
              ix01 = sort(intersect(ix01,above.thresh.1))
              ix02 = sort(intersect(ix02,above.thresh.2))
              ix0 = sort(intersect(ix01,ix02))
              if(length(ix0) <= 1) ix0 = int.size.k(rankcoef1$ix,rankcoef2$ix,2) 
           }else{ 
              if(greedy.size == 1)
                 ix0 = int.size.k(rankcoef1$ix,rankcoef2$ix,2) 
              else         
                 ix0 = int.size.k(rankcoef1$ix,rankcoef2$ix,greedy.size)                          
           }
       }
   }          
}
return(ix0)
}

obtain.newix<-function(x, y, ix1, candind, s1, s2, family, pleft, varISIS, perm, q, greedy, greedy.size, iterind){  
if(varISIS=="vanilla"){
   margcoef = margcoef(x, y, ix1, family=family, null.model=FALSE, iterind=iterind)
   rankcoef = sort(margcoef, decreasing=TRUE, index.return=TRUE) 
   if(perm == FALSE){
      if(pleft>0)
         newix = candind[rankcoef$ix[1:pleft]]
      else
         newix = NULL
   }else{
       randcoef = margcoef(x, y, ix1, family=family, null.model=TRUE, iterind=iterind)
       if(length(which(margcoef >= quantile(randcoef,q))) > 0){
          if(greedy == FALSE){
             length1 = length(which(margcoef >= quantile(randcoef,q))) 
             above.thresh = candind[rankcoef$ix[1:length1]]           
             newix = candind[rankcoef$ix[1:pleft]]                        
             newix = sort(intersect(newix,above.thresh))                
          }else
             newix = candind[rankcoef$ix[1:greedy.size]]
       }      
       else
          newix = NULL
   }
}else{        
   margcoef1 = margcoef(x[s1,], y[s1], ix1, family=family, null.model=FALSE, iterind=iterind)
   margcoef2 = margcoef(x[s2,], y[s2], ix1, family=family, null.model=FALSE, iterind=iterind)
   rankcoef1 = sort(margcoef1, decreasing=TRUE, index.return=TRUE)
   rankcoef2 = sort(margcoef2, decreasing=TRUE, index.return=TRUE)     
   if(perm == FALSE){
      if(pleft>0){
         if(varISIS=="aggr"){
            newix1 = candind[rankcoef1$ix[1:pleft]]
            newix2 = candind[rankcoef2$ix[1:pleft]]
            newix = sort(intersect(newix1,newix2))
         }
         if(varISIS=="cons"){
            iensure = intensure(pleft,l1=rankcoef1$ix,l2=rankcoef2$ix,k=pleft)
            newix1 = candind[rankcoef1$ix[1:iensure]]
            newix2 = candind[rankcoef2$ix[1:iensure]]
            newix = sort(intersect(newix1,newix2))         
         }
      }else
         newix = NULL
   }else{
       randcoef1 = margcoef(x[s1,], y[s1], ix1, family=family, null.model=TRUE, iterind=iterind)
       randcoef2 = margcoef(x[s2,], y[s2], ix1, family=family, null.model=TRUE, iterind=iterind)     
       if(length(which(margcoef1 >= quantile(randcoef1,q))) > 0 && length(which(margcoef2 >= quantile(randcoef2,q))) > 0){
          if(greedy == FALSE){
             length1 = length(which(margcoef1 >= quantile(randcoef1,q)))
             length2 = length(which(margcoef2 >= quantile(randcoef2,q)))   
             above.thresh.1 = candind[rankcoef1$ix[1:length1]]
             above.thresh.2 = candind[rankcoef2$ix[1:length2]]                                
             newix1 = candind[rankcoef1$ix[1:pleft]]
             newix2 = candind[rankcoef2$ix[1:pleft]]      
             newix1 = sort(intersect(newix1,above.thresh.1))
             newix2 = sort(intersect(newix2,above.thresh.2))
             newix = sort(intersect(newix1,newix2))                            
          }else{
             length1 = length(which(margcoef1 >= quantile(randcoef1,q)))
             length2 = length(which(margcoef2 >= quantile(randcoef2,q)))
             newix1 = candind[rankcoef1$ix[1:length1]]
             newix2 = candind[rankcoef2$ix[1:length2]]    
             iensure = intensure(greedy.size,l1=newix1,l2=newix2,k=greedy.size)
             newix = sort(intersect(newix1[1:iensure],newix2[1:iensure]))                             
          }   
       }      
       else
         newix = NULL
   } 
}
return(newix)
}
