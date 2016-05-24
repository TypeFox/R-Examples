`norm.texas` <-
function (x,writetable=F,vals="logged") {
	
   dat <- log2(x[[1]])
   
   spot.run <- levels(as.factor(x[[3]]["spotting_run",]))
   norm.dat <- c(NULL)
   arrayx <- c(NULL)
   
for (i in spot.run){
            run.cols <- which((x[[3]]["spotting_run",])==i)
            temp.dat <- dat[,run.cols]
            
            arrayx <- cbind(arrayx,x[[3]][,run.cols])
            
                cf <- apply(temp.dat,1,median)
                texas <- apply(temp.dat,2,function(x){
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
                                 
               norm.dat <- cbind(norm.dat,texas)
                                  
}
          data <- list(expression=norm.dat,
                        dummy=norm.dat,
                        arraydescription=x[[3]],
                        sampledescription=x[[4]])
                        
          if(writetable){
            write.Data(data)
            }
            
         return(data)
}

