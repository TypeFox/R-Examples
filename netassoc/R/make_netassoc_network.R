make_netassoc_network <- function(obs, nul=vegan::permatfull(obs)$perm[[1]], method="partial_correlation", args=list(method="shrinkage",verbose=FALSE), p.method="fdr", alpha=0.05, numnulls=1000, plot=TRUE,plot.legend=TRUE, plot.title=TRUE, verbose=TRUE)
{	
  if (is.matrix(obs) & is.null(dimnames(obs)))
  {
    dimnames(obs) <- list(paste("Species", 1:nrow(obs)),paste("Site", 1:ncol(obs)))
  }
  if (is.matrix(nul) & is.null(dimnames(nul)))
  {
    dimnames(nul) <- list(paste("Species", 1:nrow(nul)),paste("Site", 1:ncol(nul)))
  }
  
  obs <- as.matrix(obs,rownames.force=TRUE)
  nul <- as.matrix(nul,rownames.force=TRUE)
  
  if(!all(dim(obs)==dim(nul)))
  {
    stop("obs and nul must be same dimensionalities.")
  }

  problemspecies <- which(rowSums(obs)==0 | rowSums(nul)==0 | is.na(rowSums(nul)) | is.na(rowSums(obs)))
  if (length(problemspecies) > 0)
  {
	  stop(sprintf("Some species do not occur in any sites: %s",paste(names(problemspecies),collapse=", ")))
  }
  
  problemsites <- which(colSums(obs)==0 | colSums(nul)==0 | is.na(colSums(nul)) | is.na(colSums(obs)))
  if (length(problemsites) > 0)
  {
  	stop(sprintf("Some sites have no individuals of any species: %s",paste(names(problemsites),collapse=", ")))
  }
  
  nsp <- nrow(obs)
  
  if (plot==TRUE)
  {
    plot_netassoc_matrix(obs, colors=colorRampPalette(c('white','green'))(51),onesided=TRUE,main="Observed sp x site",legend=plot.legend,title=plot.title)
    plot_netassoc_matrix(nul, colors=colorRampPalette(c('white','black'))(51),onesided=TRUE,main="Null sp x site",legend=plot.legend,title=plot.title)
   }
  
  if (verbose==TRUE) { cat('Calculating observed co-occurrence scores...\n') }
  pcor_obs <- do.call(method, args=c(list(mat=obs), args))
  #pcor_obs <- partial_correlation(obs, method, verbose=FALSE)
  
  if (plot==TRUE)
  {
     plot_netassoc_matrix(pcor_obs, colors=colorRampPalette(c('red','white','blue'))(51),onesided=FALSE,main="Obs sp x sp",legend=plot.legend,title=plot.title)
  }

  # create null distributions
  pcor_nul_all <- array(dim=c(nsp, nsp, numnulls))
  
  for (k in 1:numnulls)
  {
    # Generating null
    if (verbose==TRUE) { cat(sprintf('Generating null replicate %d...\n', k)) }
    else
    {
      cat(".")
    }
    
    nul_resampled <- generate_nul_resample(nul, obs)
    
    #pcor_nul <- partial_correlation(nul_resampled, method, verbose=FALSE)
    pcor_nul <- do.call(method, args=c(list(mat=nul_resampled), args))
    
    pcor_nul_all[,,k] <- pcor_nul
  }
  if (verbose==FALSE)
  {
    cat("\n")
  }
  
  # calculate how often obs values are less than nul
  if (verbose==TRUE) { cat('Calculating standardized effect sizes...\n') }
  
  absobs_minus_absnul <- -1*sweep(abs(pcor_nul_all),c(1,2), abs(pcor_obs), "-")
  
  pcor_nul_mean <- apply(pcor_nul_all, c(1,2), mean)
  dimnames(pcor_nul_mean) <- dimnames(pcor_obs)
  
  pcor_nul_sd <- apply(pcor_nul_all, c(1,2), sd)
  dimnames(pcor_nul_sd) <- dimnames(pcor_obs)
  
  if (plot==TRUE)
  {
    plot_netassoc_matrix(pcor_nul_mean, colors=colorRampPalette(c('red','white','blue'))(51),onesided=FALSE,main="Nul mean sp x sp",legend=plot.legend,title=plot.title)
    plot_netassoc_matrix(pcor_nul_sd, colors=colorRampPalette(c('white','gray'))(51),onesided=TRUE,main="Nul sd sp x sp",legend=plot.legend,title=plot.title)
    
  }
  
  pcor_ses <- (pcor_obs - pcor_nul_mean) / pcor_nul_sd
  
  pcor_pvalues <- apply(absobs_minus_absnul, c(1,2), function(x) { mean(x<0) })
  diag(pcor_pvalues) <- NA
  
  # control false discovery rate
  if (verbose==TRUE) { cat('Adjusting p-values for multiple comparisons... ') }

  pcor_pvalues_adjusted <- matrix(p.adjust(pcor_pvalues, method=p.method),nrow=nsp,ncol=nsp)
  dimnames(pcor_pvalues_adjusted) <- dimnames(pcor_obs)
  
  pcor_ses_trimmed <- pcor_ses
  pcor_ses_trimmed[pcor_pvalues_adjusted > alpha] <- NA
  
  pcor_ses_trimmed_pos <- pcor_ses_trimmed
  pcor_ses_trimmed_pos[pcor_ses_trimmed_pos <= 0] <- NA
  
  pcor_ses_trimmed_neg <- pcor_ses_trimmed
  pcor_ses_trimmed_neg[pcor_ses_trimmed_neg >= 0] <- NA
  pcor_ses_trimmed_neg <- -1 * pcor_ses_trimmed_neg
  
  if (plot==TRUE)
  {
    plot_netassoc_matrix(pcor_pvalues_adjusted, gray(1-1/1:100),main="P-values (corrected)",onesided=TRUE,legend=plot.legend,title=plot.title)
    plot_netassoc_matrix(pcor_ses, colorRampPalette(c('red','white','blue'))(51),main="S.E.S. co-occurrence score for sp x sp (raw)",legend=plot.legend,title=plot.title)    
    plot_netassoc_matrix(pcor_ses_trimmed,colorRampPalette(c('red','white','blue'))(51),main="S.E.S. co-occurrence score for sp x sp (trimmed)",legend=plot.legend,title=plot.title)
  }

  
  # convert to network representation
  if (verbose==TRUE) { cat('Building network...\n') }
  # overwrite NAs to zeros
  pcor_ses_trimmed_net <- pcor_ses_trimmed
  pcor_ses_trimmed_net[is.na(pcor_ses_trimmed_net)] <- 0
  pcor_ses_trimmed_pos_net <- pcor_ses_trimmed_pos
  pcor_ses_trimmed_pos_net[is.na(pcor_ses_trimmed_pos_net)] <- 0
  pcor_ses_trimmed_neg_net <- pcor_ses_trimmed_neg
  pcor_ses_trimmed_neg_net[is.na(pcor_ses_trimmed_neg_net)] <- 0
  
  network_all <- graph.adjacency(pcor_ses_trimmed_net,mode='directed',weighted=T)
  network_pos <- graph.adjacency(pcor_ses_trimmed_pos_net,mode='directed',weighted=T)
  network_neg <- graph.adjacency(pcor_ses_trimmed_neg_net,mode='directed',weighted=T)

  if (plot==TRUE)
  {
    plot_netassoc_network(network_all,legend=plot.legend)
    box()
    if (plot.title==TRUE)
    {
      title('Association network')
    }
  }
  

  return(list(
    matrix_spsite_obs=obs,
    matrix_spsite_nul=nul,
    matrix_spsp_obs=pcor_obs,
    matrix_spsp_ses_all=pcor_ses,
    matrix_spsp_ses_thresholded=pcor_ses_trimmed,
    matrix_spsp_pvalue=pcor_pvalues_adjusted,
    network_all=network_all,
    network_pos=network_pos,
    network_neg=network_neg
    ))

}






