`norm.static.I` <-
function (x,writetable=F,vals="logged") {
         
   dat <- log2(x[[1]])
         
         
   pads <-levels(as.factor(x[[3]]["pad",]))
   
   spot.run <- levels(as.factor(x[[3]]["spotting_run",]))
   
   norm.dat <- c(NULL)
   
   arrayx <- c(NULL)
       for (j in spot.run){
         for (i in pads){
         pad.cols <- which((x[[3]]["pad",])==i
                     &(x[[3]]["spotting_run",]==j))
                     
         temp <- dat[,pad.cols]
         arrayx <- cbind(arrayx,x[[3]][,pad.cols])
         prot.col <- which (x[[3]]["target",]=="protein" 
                              & x[[3]]["pad",]==i
                              & x[[3]]["spotting_run",]==j)
                              
               normalizer <- dat[,prot.col]

               
                temp.i <- apply(temp,2,function(x){
                              x-normalizer
                              })
                              
                ## change data scaling back to native values
                if(vals=="native"){
                   normalizer.median <- median(dat[,prot.col])
                   ## add median value of normalizer values
                   temp.i <- temp.i+normalizer.median
                   ## transform to native data
                   temp.i <- 2^(temp.i)
                }
                   
                      norm.dat <- cbind(norm.dat,temp.i)
                      }
            }
            
         data <- list(expression=norm.dat,
                     dummy=norm.dat,
                     arraydescription=arrayx,
                     sampledescription=x[[4]])
         
         if(writetable){
            write.Data(data)
            }
         
         return(data)
}

