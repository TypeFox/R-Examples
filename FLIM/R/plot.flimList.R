plot.flimList <-
function(x, response, ...) {
  if(is.null(x$flim.obj$lambda)) {
    par(mfrow=c(2,2))
    times <- x$flim.obj$times
    k <- which(x$flim.obj$info$responses == response)
    fits <- x$flim.obj$fit
    models <- x$flim.obj$info$reg.fmlas
    cat("  Model: ")
    print(models[[k]], showEnv=F)
    par(ask = TRUE)
    for(i in 1:(length(times)-1)) {
      fit <- fits[[i]][[k]]
      plot(fit, ...)
      cat(paste0("Time: ", times[i]))
    }
    par(mfrow=c(1,1))
    par(ask = FALSE)
  } else {print("No plot options for flimList with non NULL lambda")} 
}
