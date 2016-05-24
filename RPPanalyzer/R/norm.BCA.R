`norm.BCA` <-
function(x,proteinc="BCA",writetable=F,vals="logged"){
         
         dat <- log2(x[[1]])
         
         cf <- log2(x[[4]][,proteinc])
         
         texas <- apply(dat,2,function(x){
                                    x-cf
                                    })
                                    
         ## change data scaling back to native values
                if(vals=="native"){
                   normalizer.median <- median(cf)
                   ## add median value of normalizer values
                   texas <- texas + normalizer.median
                   ## transform to native data
                   texas <- 2^(texas)
                }
                
         data <- list(  expression=texas,
                        dummy=texas,
                        arraydescription=x[[3]],
                        sampledescription=x[[4]])
                        
          if(writetable){
            write.Data(data)
            }
            
         return(data)
}

