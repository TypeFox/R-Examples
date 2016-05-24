`norm.static.II` <-
function (x,normalizer=c("housekeeping"),writetable=F,vals="logged") {
         
         dat <- log2(x[[1]])
         
         prot.cols <- which(x[[3]]["target",]==normalizer)
         
               ifelse ( length(prot.cols) > 1,
                     norm <- apply(dat[,prot.cols],1,mean),
                     norm <- dat[,prot.cols])
               
                temp <- apply(dat,2,function(x){
                              x-norm
                              })
                              
                ## optional: change data scaling back to native values
                if(vals=="native"){
                   normalizer.median <- median(norm)
                   ## add median value of normalizer values
                   temp <- temp+normalizer.median
                   ## transform to native data
                   temp <- 2^(temp)
                }
                
                      data <- list(expression=temp,
                                    dummy=temp,
                                   arraydescription=x[[3]],
                                   sampledescription=x[[4]])
                                   
                      if(writetable){
                        write.Data(data)
                        }
                        
                      return(data)
}

