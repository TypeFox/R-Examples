`gq_diagnostics` <- 
function (gq_output, GR = TRUE, Geweke = TRUE, Heidel = TRUE,
    acf = TRUE, lags = c(1:6, 10, 13, 16, 20), silent = TRUE,
    transform = FALSE)
{
    if (!is.null(gq_output$draws)) {
        gq_output <- list(gq_output)
    }
    if (is.null(gq_output[[1]]$draws)) {
        stop("gq_output must be either a list of or a single return object from 'Analyze'")
    }
    gq_draws <- mcmc.list(mcmc(1))
    p <- ncol(gq_output[[1]]$draws$mu)
    R <- gq_output[[1]]$numrows_pt
    C <- gq_output[[1]]$numcols_pt
    dlist <- list(NULL)
    lname <- NULL
    sofar <- 1
    if (GR) {
        cat("Computing Gelman-Rubin diagnostic...")
        gq_subdraws <- mcmc.list(mcmc(1))
        NNs_subset <- matrix(1:(R*C),nrow=R,ncol=C,byrow=TRUE)[1:(R-1),1:(C-1)]
        for (i in 1:length(gq_output)) {
            tmp_chain <- cbind(gq_output[[i]]$draws$mu, 
			       gq_output[[i]]$draws$Sigma,
                               gq_output[[i]]$draws$NNs[,NNs_subset])

              # Ensure correct handling of the column names:
  	      colnames(tmp_chain) <- 
		c(colnames(gq_output[[i]]$draws$mu),
		  colnames(gq_output[[i]]$draws$Sigma),
		  colnames(gq_output[[i]]$draws$NNs)[NNs_subset])
            gq_subdraws[[i]] <- mcmc(tmp_chain)
        }
        tmp_out <- try(gelman.diag(gq_subdraws, transform = transform),
            silent = silent)
        if (class(tmp_out) != "try-error") {
            dlist[[sofar]] <- tmp_out
            cat("successful.\n")
        }
        else {
            dlist[[sofar]] <- NULL
            cat("failed.\n")
        }
        sofar <- sofar + 1
        lname <- c(lname, "gelman.diag")
        rm(gq_subdraws)
    }
    for (i in 1:length(gq_output)) {
        tmp_chain <- cbind(gq_output[[i]]$draws$mu, 
			   gq_output[[i]]$draws$Sigma,
            		   gq_output[[i]]$draws$NNs, 
			   gq_output[[i]]$draws$LAMBDA,
            		   gq_output[[i]]$draws$TURNOUT, 
			   gq_output[[i]]$draws$GAMMA)
        gq_draws[[i]] <- mcmc(tmp_chain)
    }
    if (Geweke) {
        cat("Computing Geweke diagnostic...")
        tmp_out <- try(geweke.diag(gq_draws), silent = silent)
        if (class(tmp_out) != "try-error") {
            dlist[[sofar]] <- tmp_out
            cat("successful.\n")
        }
        else {
            dlist[[sofar]] <- NULL
            cat("failed.\n")
        }
        sofar <- sofar + 1
        lname <- c(lname, "geweke.diag")
    }
    if (Heidel) {
        cat("Computing Heidelberger-Welch diagnostic...")
        tmp_out <- try(heidel.diag(gq_draws), silent = silent)
        if (class(tmp_out) != "try-error") {
            dlist[[sofar]] <- tmp_out
            cat("successful.\n")
        }
        else {
            dlist[[sofar]] <- NULL
            cat("failed.\n")
        }
        sofar <- sofar + 1
        lname <- c(lname, "heidel.diag")
    }
    if (acf) {
        cat("Computing autocorrelations...")
        tmp_out <- try(autocorr.diag(gq_draws, lags = lags),
            silent = silent)
        if (class(tmp_out) != "try-error") {
            dlist[[sofar]] <- tmp_out
            cat("successful.\n")
        }
        else {
            dlist[[sofar]] <- NULL
            cat("failed.\n")
        }
        sofar <- sofar + 1
        lname <- c(lname, "autocorr")
    }
    names(dlist) <- lname
    return(dlist)
}
