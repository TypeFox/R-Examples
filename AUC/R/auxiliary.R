#labels need to be a factor, predictons need to be numeric
.confusionMatrix <- function(predictions,labels,perc.rank, measure) {
              
               
              #This function is a scaled down faster version of parts of code of the ROCR package.
              #It has less functionality and less error handling and is focused on speed.
              #For more functionality (e.g., averaging cross validaton runs) see ROCR
              
              if (measure=='sensitivity') {
                if (perc.rank==TRUE) predictions <- rank(predictions,ties.method="min")/length(predictions)
              } else if (measure =='specificity' || measure =='accuracy') {
                if (perc.rank==TRUE) predictions <- rank(predictions,ties.method="max")/length(predictions)
              }
  
              levels <- sort(levels(labels))
  
              labels <- ordered(labels,levels=levels)
                

              #  if (length(levels) != 2) {
              #    message <- paste("Number of classes should be 2.\n")
              #    stop(message)
              #  }
              
              
              ## compute cutoff/fp/tp data
              
              cutoffs <- numeric()
              fp <- numeric()
              tp <- numeric()
              fn <- numeric()
              tn <- numeric()
              n.pos <- numeric()
              n.neg <- numeric()
              n.pos.pred <- numeric()
              n.neg.pred <- numeric()
   
              n.pos <-  sum( labels == levels[2] )
              n.neg <-  sum( labels == levels[1] )
                        
              pos.label <- levels(labels)[2]
              neg.label <- levels(labels)[1]
              
              pred.order <- order(predictions, decreasing=TRUE)
              predictions.sorted <- predictions[pred.order]
              tp <- cumsum(labels[pred.order]==pos.label) #predicted to be positive, and in reality they are not
              fp <- cumsum(labels[pred.order]==neg.label) #predicted to be positive, but in reality they are not (i.e., negative since there are 2 classes)
              
              ## remove fp & tp for duplicated predictions
              ## as duplicated keeps the first occurrence, but we want the last, two rev are used.
              dups <- rev(duplicated(rev(predictions.sorted)))
              tp <- c(0, tp[!dups])
              fp <- c(0, fp[!dups])
              cutoffs <- c(1, predictions.sorted[!dups])
              

              fn <- n.pos - tp
              tn <- n.neg - fp
              n.pos.pred <- tp + fp
              n.neg.pred <- tn + fn

              list(cutoffs=cutoffs,tp=tp, fp=fp, tn=tn, fn=fn, n.pos=n.pos, n.neg=n.neg, n.pos.pred=n.pos.pred, n.neg.pred=n.neg.pred)
            
}