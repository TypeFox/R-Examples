flimList <-
function(fo) {
  fits <- fo$fit
  info <- fo$info
  times <- fo$times
  est.list <- list()
  sum.list <- list()
  models <- info$reg.fmlas
  if(is.null(fo$lambda)) {
    for(k in 1:length(models)) {
      mt <- length(times)-1
      est.mat <- matrix(0, mt, length(fits[[mt]][[k]]$coef))
      est.mat <- as.data.frame(est.mat)
      names(est.mat) <- names(fits[[mt]][[k]]$coef)
      rownames(est.mat) <- times[1:(mt)]   
      est.names.k <- names(est.mat)
      sum.list[[k]] <- matrix(0, nrow = 0, ncol = 4)
      for(i in 1:mt) {
        if(!is.null(fits[[i]][[k]])) {
          est.mat[i, ] <- round(fits[[i]][[k]]$coef, 3)
          sum.list[[k]] <- rbind(sum.list[[k]], summary(fits[[i]][[k]])$coef)
        } else {est.mat[i, ] <- NA}
      }
      est.list[[k]] <- est.mat
    }
    retme <- list(flim.obj = fo,
                  est.list = est.list,
                  sum.list = sum.list,
                  models = models,
                  coef = est.list)
    class(retme) <- "flimList"
    retme
  } else {
    for(k in 1:length(models)) {
      mt <- length(times)-1
      est.mat <- matrix(0, mt, length(coef(fits[[mt]][[k]])))
      est.mat <- as.data.frame(est.mat)
      names(est.mat) <- names(coef(fits[[mt]][[k]]))
      if(fits[[mt]][[k]]$Inter==1) names(est.mat)[1] <-"(Intercept)"
      rownames(est.mat) <- times[1:(mt)]   
      est.names.k <- names(est.mat)
      for(i in 1:mt) {
        if(!is.null(fits[[i]][[k]])) {
          est.mat[i, ] <- round(coef(fits[[i]][[k]]), 3)
        } else {est.mat[i, ] <- NA}
      }
      est.list[[k]] <- est.mat
    }
    retme <- list(flim.obj = fo,
                  est.list = est.list,
                  #sum.list = sum.list,
                  models = models)
    class(retme) <- "flimList"
    retme
  }
}
