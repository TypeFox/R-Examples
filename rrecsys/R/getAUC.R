#Global AUC(area under the ROC curve).

# Reference: T. Fawcett. ROC Graphs: Notes and Practical Considerations for Data Mining Researchers.


setMethod("getAUC", signature = c(model = "evalModel"), function(model, alg, ...) {
    if (missing(alg)) {
        stop("Evaluation on recommendations can not proceed if argument alg and is not specified.")
    }
    
    auc <- c()
    
    for (i in 1:model@folds) {
        
        ptm <- Sys.time()
        
        copy_data <- model@data
        copy_data@data[model@fold_indices[[i]]] <- 0
        
        copy_data@data <- matrix(copy_data@data, nrow(model@data))
        
        r <- rrecsys(copy_data, alg = alg, ...)
        
        p <- predict(r)
        
        auc <- c(auc, 0)
        
        temp <- apply(model@data@data, 1, function(x) which(x != 1))
        for (u in 1:nrow(model@data)) {
            # if there are no pairs then the default auc is 0.5
            if ((length(temp[[u]]) == 0) | (length(model@fold_indices_x_user[[i]][[u]]) == 0)) {
                auc_on_user <- 0.5
            } else {
                auc_on_user <- 0
                # extract and compair all the items in the folds and compare their predicted value with the predicted value of an urated item.
                for (m in model@fold_indices_x_user[[i]][[u]]) {
                  for (n in temp[[u]]) {
                    
                    if (p[u, m] > p[u, n]) {
                      auc_on_user <- auc_on_user + 1
                    }
                    
                  }
                }
                # divede by the number of total checked pairs to get the real auc for the user.
                auc_on_user <- auc_on_user/(length(model@fold_indices_x_user[[i]][[u]]) * length(temp[[u]]))
                
            }
            auc[i] <- auc[i] + auc_on_user
            
        }
        # determine the average auc.
        auc[i] <- auc[i]/nrow(model@data)
        cat("\nFold:", i, "/", model@folds, " elapsed. Time:", as.numeric(Sys.time() - ptm, units = "secs"), "\n\n")
    }
    
    auc <- as.data.frame(list(AUC = auc))
    
    auc <- rbind(auc, colMeans(auc))
    
    row.names(auc) <- c(paste0(1:model@folds, rep("-fold", model@folds)), "Average")
    
    auc
    
    
}) 
