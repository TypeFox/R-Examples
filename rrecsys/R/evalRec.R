setMethod("evalRec", signature = c(model = "evalModel"), function(model, alg = NULL, topN = 3, goodRating = NULL, alpha = 10, ...) {
    
    
    
    if (missing(alg)) {
        stop("Evaluation on recommendations can not proceed if argument alg and is not specified.")
    }
    
    if (model@data@binary) {
        if (!is.null(goodRating)) {
            cat("NOTE: The \"goodRating\" attribute is not needed for binary dataset.\n")
        }
        goodRating <- 1
    }
    
    if (missing(goodRating)) {
        stop("Evaluation on recommendations can not proceed if argument goodRating and is not specified.")
    }
    
    cat("Evaluating top-", topN, " recommendation on ", rrecsysRegistry$get_entry(alg = alg)$alg, ".\n")
    
    
    nusers <- nrow(model@data)
    # initialize empty containers
    res <- NULL
    nDCG <- rep(0, model@folds)
    rankscore <- rep(0, model@folds)
    
    # iterations on ofolder
    for (i in 1:model@folds) {
        
        ptm <- Sys.time()
        # geta a copy of the rating matrix
        copy_data <- model@data
        # generate the training set
        copy_data@data[model@fold_indices[[i]]] <- 0
        copy_data@data <- matrix(copy_data@data, nusers)
        # train the train set
        
        
        r <- rrecsys(copy_data, alg = alg, ...)
        # get the recommended indices
        rec <- recommend(r, topN)@indices
        
        # item and user coverage calculation
        tot_rec <- 0
        item_coverage <- rep(FALSE, ncol(model@data))
        
        res_on_fold <- list(TP = 0, FP = 0, TN = 0, FN = 0, precision = 0, recall = 0, F1 = 0)
        
        for (m in 1:nusers) {
            #item coverage  
            for (n in rec[[m]]) {
                item_coverage[n] <- TRUE
                tot_rec <- tot_rec + 1
            }
          
            #determine results on user. 
            res_user <- getPrecRecall(model@data@data[m, ], rec[[m]], model@fold_indices_x_user[[i]][[m]], goodRating)
            
            for(j in 1:length(res_on_fold)) {
              res_on_fold[[j]] <- res_on_fold[[j]] + res_user[[j]]
            }
            
            nDCG[i] <- nDCG[i] + 
              getnDCG(model@data@data[m, ], rec[[m]], model@fold_indices_x_user[[i]][[m]], goodRating)
            
            rankscore[i] <- rankscore[i] + 
              getrankscore(model@data@data[m, ], rec[[m]], model@fold_indices_x_user[[i]][[m]], goodRating, alpha)
        }
        
        res_on_fold <- lapply(res_on_fold, function(x) x/nusers)
        
        # find average values for the confusion matrix.
        res <- rbind(res, as.data.frame(res_on_fold))
        
        #get averages
        nDCG[i] <- nDCG[i]/nusers
        rankscore[i] <- rankscore[i]/nusers
        
        cat("\nFold:", i, "/", model@folds, " elapsed. Time:", as.numeric(Sys.time() - ptm, units = "secs"), "\n\n")

        
    }
    
    res <- cbind(res, nDCG, rankscore)
    
    res <- rbind(res, colMeans(res))
    
    row.names(res) <- c(paste0(1:model@folds, rep("-fold", model@folds)), "Average")
    
    show(r)
    
    cat("\nItem coverage:", 100 * sum(item_coverage)/ncol(model@data), "%.\n")
    cat("\nUser coverage:", 100 * tot_rec/(topN * nrow(model@data)), "%.\n")
    res
    
})


