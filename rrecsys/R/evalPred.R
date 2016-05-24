setMethod("evalPred", signature = c(model = "evalModel"), function(model, alg = NULL, ...) {
    
    if (is.null(alg)) {
        stop("Evaluation on recommendations cannot proceed if argument alg and is not specified.")
    }
    
    nusers <- nrow(model@data@data)
    
    #user based RMSE & MAE
    uRMSE <- c()
    uMAE <- c()
    #global RMSE & MAE
    gRMSE <- c()
    gMAE <- c()
    
    for (i in 1:model@folds) {
        
        ptm <- Sys.time()
        
        copy_data <- model@data
        copy_data@data[model@fold_indices[[i]]] <- 0
        
        copy_data@data <- matrix(copy_data@data, nusers)
        
        r <- rrecsys(copy_data, alg = alg, ...)
        
        # mae & rmse on user
        predictions <- predict(r, Round = FALSE)
        
        # calculation on MAE and RMSE on each user
        users_rmse <- 0
        users_mae <- 0

        for (n in 1:nrow(model@data)) {
          temp <- abs(model@data@data[n, model@fold_indices_x_user[[i]][[n]]] - predictions[n, model@fold_indices_x_user[[i]][[n]]])

          if(length(temp) == 0) next
          
          users_rmse <- users_rmse + sqrt(sum((temp)^2)/length(model@fold_indices_x_user[[i]][[n]]))
          users_mae <- users_mae + sum((temp))/length(model@fold_indices_x_user[[i]][[n]])
        }
        # derivate an average MAE and RMSE on the whole rating matrix
        uRMSE <- c(uRMSE, users_rmse/nrow(model@data))
        uMAE <- c(uMAE, users_mae/nrow(model@data))
        
        # calculation on global MAE and RMSE
        temp <- model@data@data[model@fold_indices[[i]]] - predictions[model@fold_indices[[i]]]
        
        mae_i <- sum(abs(temp)) / length(model@fold_indices[[i]])
        
        rmse_i <- sqrt(sum(temp^2) / length(model@fold_indices[[i]]))
        
        gMAE <- c(gMAE, mae_i)
        gRMSE <- c(gRMSE, rmse_i)
        
        
        cat("\nFold:", i, "/", model@folds, " elapsed. Time:", as.numeric(Sys.time() - ptm, units = "secs"), "\n\n")

    }
    
    # average on folds
    uRMSE <- c(uRMSE, sum(uRMSE)/model@folds)
    uMAE <- c(uMAE, sum(uMAE)/model@folds)
    gRMSE <- c(gRMSE, sum(gRMSE)/model@folds)
    gMAE <- c(gMAE, sum(gMAE)/model@folds)
    
    
    names(uRMSE) <- c(paste0(1:model@folds, rep("-fold", model@folds)), "Average")
    names(uMAE) <- c(paste0(1:model@folds, rep("-fold", model@folds)), "Average")
    names(gRMSE) <- c(paste0(1:model@folds, rep("-fold", model@folds)), "Average")
    names(gMAE) <- c(paste0(1:model@folds, rep("-fold", model@folds)), "Average")
    
    show(r)
    cat("\n")
    as.data.frame(list(MAE = uMAE, RMSE = uRMSE, globalMAE = gMAE, globalRMSE = gRMSE))
    
}) 
