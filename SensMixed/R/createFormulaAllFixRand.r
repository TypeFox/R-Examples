

#create formula with all random and fixed effects and their interactions
createFormulaAllFixRand  <-  function(structure, data, response, fixed,
                                      random, corr)
{
   
   f1 <- function(x){
     is.factor(data[, which(colnames(data) == x)])
   }
  
   formula1 <- NULL
   is.cov.present <- FALSE
   covs <- fixed$Product[which(unlist(lapply(fixed$Product,f1))!=TRUE)]
   if(length(covs) > 0)
     is.cov.present <- TRUE
   formula.covs <- paste(covs, collapse="+")
     
     
   
   for (j in 1:length(c(fixed$Product, fixed$Consumer)))
   {
         #create all combinations of j'th order of fixed effects
         if(!is.list(structure) && (structure==1 & j>1))
           break
         else if(!is.list(structure) && (structure==2 & j>2))
           break
         else if(is.list(structure) && (structure$product_structure==1 & j>1))
           break
         else if(is.list(structure) && (structure$product_structure==2 & j>2))
           break
         inter.fix <- combn(c(fixed$Product, fixed$Consumer), j, simplify = FALSE)
         
         #put in formula all fixed factor combinations of j'th order
         for(k in 1:length(inter.fix))
         {
            eff.fix <- paste(inter.fix[[k]],collapse=":")
            if(checkComb(data,inter.fix[[k]]))
              next()
            if(is.null(formula1))
               formula1 <- paste(response, "~", eff.fix, sep="")
            else 
               formula1 <- paste(formula1,"+", eff.fix, sep="")
            
            
            # not to include factors associated with Consumer into random part
            if(length(which((inter.fix[[k]] %in% fixed$Consumer)==TRUE))>0)
              next
            #create all combinations of i'th order of random effects (interactions with fixed)
            for(i in 1:nrandTerms(random, structure))
            {
               #for sensmixed procedure
               if(is.list(random))
               {
                 if(structure$error_structure=="3-WAY")
                   inter.rand <- combn(unlist(random),i,simplify = FALSE)
                 else
                   inter.rand <- combn(random$individual,i,simplify = FALSE)
               }
               else inter.rand <- combn(random,i,simplify = FALSE)
               #put in formula all random (interactions with fixed) factor combinations of i'th order
               for(l in 1:length(inter.rand))
               { 
                                  
                 if(checkComb(data,c(inter.rand[[l]],inter.fix[[k]])))
                   next()
                 # if the consumer attributes are present in a fixed part
                 eff.rand <- paste(inter.rand[[l]],collapse=":")
                 ind.cov <- which(unlist(lapply(inter.fix[[k]],f1))!=TRUE)
                 if(length(ind.cov)==0)
                 {
                   # if the attributes associated with Consumer are present then
                   # eliminate this attribute from eff.fix
                   #eff.fix.rand <- if(!is.null(fixed$Consumer)) paste(inter.fix[[k]][!inter.fix[[k]] %in% fixed$Consumer], collapse=":") else eff.fix
                   eff.fix.rand  <-  eff.fix
                   #if(eff.fix.rand=="")
                   # next
                   if(is.cov.present)
                   {
                     if(corr)
                       formula1 <- paste(formula1,"+", "(1+",formula.covs,"|",
                                         eff.fix.rand,":",eff.rand,")",sep="")
                     else
                       formula1 <- paste(formula1,"+", 
                                         paste(lapply(covs, function(x) 
                                           paste("(0+",x,"|",eff.fix.rand,":",eff.rand,")",
                                                 sep="")),
                                           collapse="+"))
                     
                   }                      
                   if(!is.cov.present || !corr)
                      formula1 <- paste(formula1,"+", 
                                        "(1|",eff.fix.rand,":",eff.rand,")",
                                        sep="")
                 }
                    
               }
             }
          }
     }
     
   
     #create all combinations of i'th order of random effects 
     for(i in 1:nrandTerms(random, structure))
     {
       #for sensmixed procedure
       if(is.list(random))
       {
         if(structure$error_structure=="3-WAY")
           inter.rand <- combn(unlist(random),i,simplify = FALSE)
         else
           inter.rand <- combn(random$individual,i,simplify = FALSE)
       }
       else inter.rand <- combn(random,i,simplify = FALSE)
        #put in formula all random factor combinations of i'th order
        for(l in 1:length(inter.rand))
        { 
              
          if(checkComb(data,inter.rand[[l]]))
             next()
          
           eff.rand <- paste(inter.rand[[l]],collapse=":")
          
          if(is.cov.present)
          {
            if(corr)
              formula1 <- paste(formula1,"+", 
                                "(1+",formula.covs,"|",eff.rand,")",sep="")
            else
              formula1 <- paste(formula1,"+", paste(lapply(covs, function(x) 
                paste("(0+",x,"|",eff.rand,")",sep="")),collapse="+"))
           }            
          if(!is.cov.present || !corr)
            formula1 <- paste(formula1,"+", "(1|",eff.rand,")",sep="")
               
        }
     }     
   ## added code for 2.0-7:
   ## 2-WAY No-Rep + Replication effect + Assessor:Replication
   if(is.list(random) && structure$error_structure == "2-WAY"){
     formula1 <- paste(formula1,"+", "(1|",random$replication,")","+", 
                       "(1|",random$individual, ":", random$replication, ")",
                       sep="")
   }
 return(as.formula(formula1))
}


nrandTerms <- function(random, structure)
{
  if (is.list(random) && structure$error_structure!="3-WAY") 
    return(length(random$individual))
   else return(length(random))
}
