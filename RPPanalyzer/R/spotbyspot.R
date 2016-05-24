
spotbyspot <- function(x,normalizer="housekeeping",vals="native"){

   ## remove blank and protein stain values
   tempx <- remove.arrays(x,param="target",arrays2rm=c("protein","blank"))
   
   ## log expression data
   tempx[[1]] <- log2(tempx[[1]])
   tempx[[2]] <- log2(tempx[[2]])
   ## housekeeping object
   hk.cols <- tempx[[3]]["target",]==normalizer
   
   hk <- vector("list",4)
   hk[[1]] <- tempx[[1]][,hk.cols]
   hk[[2]] <- tempx[[2]][,hk.cols]
   hk[[3]] <- tempx[[3]][,hk.cols]
   hk[[4]] <- tempx[[4]]
   

   ## target object
   tgt <- remove.arrays(tempx,param="target",arrays2rm=normalizer)
   
   ## order housekeep and target objects by column
   o.hk <- order(hk[[3]]["pad",],hk[[3]]["slide",],hk[[3]]["spotting_run",])
   hk[[1]] <- hk[[1]][,o.hk]
   hk[[2]] <- hk[[2]][,o.hk]
   hk[[3]] <- hk[[3]][,o.hk]
   
   
   o.tgt <- order(tgt[[3]]["pad",],tgt[[3]]["slide",],tgt[[3]]["spotting_run",])
   tgt[[1]] <- tgt[[1]][,o.tgt]
   tgt[[2]] <- tgt[[2]][,o.tgt]
   tgt[[3]] <- tgt[[3]][,o.tgt]
   
   ## normalize expression values
   tgt[[1]] <-  tgt[[1]]-hk[[1]]

   
   ## optional: change data scaling back to native values
                if(vals=="native"){
                   Fmedian <- median(hk[[1]])

                   ## add median value of normalizer values
                   tgt[[1]] <- tgt[[1]]+Fmedian

                   ## transform to native data
                   tgt[[1]] <- 2^(tgt[[1]])
                   tgt[[2]] <- 2^(tgt[[2]])
                }
   return(tgt)
   }
   
