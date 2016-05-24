plotINFO <- function(ermobject, type = "both", theta = seq(-6, 6, length.out = 1001L), ...){

  ddd <- list(...)
  get_dots <- function(element, default){ ifelse(!is.null(ddd[[element]]), ddd[[element]], default) }
  extraVars <- c("mainI", "mainT", "ylabI", "ylabT", "xlab", "legpos")
  if(any(!(names(ddd) %in% extraVars))){ warning("additional argument(s) ignored.") }

  type <- match.arg(type, c("item", "test", "both"))

  if(type == "both"){
    old_pars <- par(mfrow=c(2L, 1L), no.readonly = TRUE)
    on.exit(par(old_pars))
  }

  if(type %in% c("item", "both")){
    iinfo <- item_info(ermobject, theta)
    info <- lapply(iinfo, function(x) x$i.info)
    pltinfo <- matrix(unlist(info), ncol = ncol(ermobject$X))
    matplot(x = theta, y = pltinfo, type = "l", main = get_dots("mainI", "Item Information"), xlab = get_dots("xlab", "Latent Trait"), ylab = get_dots("ylabI", "Information"))
    itmnames <- paste("Item", seq_len(ncol(ermobject$X)))
    legend(x = get_dots("legpos", "topright"), legend = itmnames, pch = NULL, lty = 1:5, col = 1:6)
  }

  if(type %in% c("item", "both")){
    tinfo <- test_info(ermobject, theta)
    plot(x = theta, y = tinfo, type = "l", main = get_dots("mainT", "Test Information"), xlab = get_dots("xlab", "Latent Trait"), ylab = get_dots("ylabT", "Scale Information"))
  }

}
