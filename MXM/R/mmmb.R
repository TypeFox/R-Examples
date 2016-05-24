mmmb = function(target, dataset , max_k = 3 , threshold = 0.05 , test = "testIndFisher" , user_test = NULL, robust = FALSE, ncores = 1, hold = FALSE) {

   durat = proc.time()  
   
   mmpcobject = MMPC(target, dataset, max_k = 3, threshold = threshold, test = test , user_test = user_test, hash = TRUE, robust = robust, ncores = 1, backward = FALSE)
   stat_hash = mmpcobject@hashObject$stat_hash
   pvalue_hash  = mmpcobject@hashObject$pvalue_hash
   varsToIterate = mmpcobject@selectedVars;
   met = 1:length(varsToIterate)
   
   lista = list()
   if ( length(varsToIterate) > 0 ) {
     for( i in 1:length(met) ) {
        tar <- dataset[ , varsToIterate[i] ];
        datas <- cbind( target, dataset[, -varsToIterate[i] ] )
        res = MMPC(tar, datas, max_k = 3, threshold = threshold, test = test , user_test = user_test, hash = FALSE, robust = robust, ncores = 1, backward = FALSE) 

        lista[[ i ]] = res@selectedVars
        
        if (hold == FALSE ) {
          if ( 1 %in% res@selectedVars == FALSE ) {
            met[i] = 0;
          }
        }
        
     }
   }
  
   pct = varsToIterate[met]
   lista = unlist( lista[met] )
   lista = unique(lista)
   lista = setdiff(lista, c(1, pct) )  ## remove the target and the PC set from the candidate MB set
   ina = 1:length(lista)
   
   if ( length(ina)>0 ) {
   
     ci_test = test = mmpcobject@test
     av_tests = c("testIndFisher", "testIndSpearman", "gSquare", NULL);
     
     #cat(test)
     
     if( length(test) == 1 ) #avoid vectors, matrices etc
     {
       test = match.arg(test , av_tests ,TRUE);
       #convert to closure type
       if(test == "testIndFisher")
       {
         #an einai posostiaio target
         if ( min(target) > 0 & max(target) < 1 ) {
           target = log( target/(1 - target) ) ## logistic normal 
         }

         test = testIndFisher;
       }
       else if(test == "testIndSpearman")
       {
         #an einai posostiaio target
         if ( min(target) > 0 & max(target) < 1 ) {
           target = log( target / (1 - target) ) ## logistic normal 
         }
         
         target = rank(target)
         dataset = apply(dataset, 2, rank)  
         test = testIndSpearman;  ## Spearman is Pearson on the ranks of the data

       }
       else if(test == "gSquare")
       {
         test = gSquare;
       }
       #more tests here
     }else{
       stop('invalid test option');
     }

     ps = mmpcobject@hashObject$pvalue_hash 
     a = log(threshold)
     p = as.matrix( as.list(ps) )

     ## Phase II of MMMB
     ## search for s such that Ind(X, T| s)
     ## For every potential child y if Dep(X, T| s, y), then X belongs to MB
    
     for ( l in ina ) {

       xIndex = lista[l]
       inds = grep( paste(xIndex, ""), rownames(p) )
       keys = rownames(p)[inds]
       keys = strsplit(keys," ") 

       cs = list()
       if (length(keys) > 0) {
         for ( i in 1:length(keys) )
         {
           vim = keys[[i]]
           #print(key)
           if ( vim[1] == as.character(xIndex) ) {
             if ( p[ inds[i], 1 ] > a) {
               cs[[ i ]] = as.numeric(vim[-1])
             }
           }
         }
       } else{ 
         cs[[ i ]] = 0
       }

       if ( length(cs) > 0 ) {
         ## need to speed it up here a bit, first the least associated cs and the most associated parents and children
         trials = numeric( length(cs) )
         for ( m in 1:length(pct) ) {
           ela = test(target, dataset, xIndex = xIndex, csIndex = as.vector( c( pct[m], trials[[ 1 ]] ) ), dataInfo=NULL, univariateModels=NULL, hash=FALSE, stat_hash=NULL, pvalue_hash=NULL, robust = robust)$pvalue
           k = 1
           while ( ela > 0.05 & k < length(trials) ) {
             k = k + 1
             ela = test(target, dataset, xIndex = xIndex, csIndex = as.vector( c( pct[m], trials[[ k ]] ) ), dataInfo=NULL, univariateModels=NULL, hash=FALSE, stat_hash=NULL, pvalue_hash=NULL, robust = robust)$pvalue
           }  
         }
         
       } else {
         ## need to speed it up here a bit, first the least associated cs and the most associated parents and children
         for ( m in 1:length(pct) ) {
           ela = test(target, dataset, xIndex = xIndex, csIndex = pct[m], dataInfo=NULL, univariateModels=NULL, hash=FALSE, stat_hash=NULL, pvalue_hash=NULL, robust = robust)$pvalue
           k = 1
           while ( ela > 0.05 ) {
             k = k + 1
             ela = test(target, dataset, xIndex = xIndex, csIndex = pct[m], dataInfo=NULL, univariateModels=NULL, hash=FALSE, stat_hash=NULL, pvalue_hash=NULL, robust = robust)$pvalue
           }  
         }
       }
       if ( ela > a ) 
       ina[l] = 0
     }

  }   
  
  runtime = proc.time() - durat   
  list( mb = mmpcobject@selectedVars[met], ci_test = ci_test, runtime = runtime )
}