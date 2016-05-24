summary.simmr_output =
function(object,type=c('diagnostics','quantiles','statistics','correlations'),group=1,...) {

  # Get the specified type
  type=match.arg(type,several.ok=TRUE)

  out_bgr = out_quantiles = out_statistics = out_cor = vector('list',length=length(group))
  names(out_bgr) = paste0('group_',group)
  names(out_quantiles) = paste0('group_',group)
  names(out_statistics) = paste0('group_',group)
  names(out_cor) = paste0('group_',group)
  
  # Loop through groups
  for(i in 1:length(group)) {
    
    if(length(group)>1) cat(paste("\nSummary for group",group[i],'\n'))
    out_all = do.call(rbind,object$output[[group[i]]])
    
    # Get objects
    out_bgr[[i]] = coda::gelman.diag(object$output[[group[i]]],multivariate=FALSE)$psrf
    out_quantiles[[i]] = t(apply(out_all,2,'quantile',probs=c(0.025,0.25,0.5,0.75,0.975)))
    #  coda:::summary.mcmc.list(object$output)$quantiles
    out_statistics[[i]] = t(apply(out_all,2,function(x) {return(c(mean=mean(x),sd=stats::sd(x)))}))
    # coda:::summary.mcmc.list(object$output)$statistics[,1:2]
    out_cor[[i]] = stats::cor(out_all)
    
    if ('diagnostics'%in%type) {
      # Print out gelman diagnostics of the output
      cat('Gelman diagnostics - these values should all be close to 1.\n')
      cat('If not, try a longer run of simmr_mcmc.\n')
      print(round(out_bgr[[i]],2))
    }
    
    if ('quantiles'%in%type) {
      # Print out quantiles argument
      print(round(out_quantiles[[i]],3))
    }
    
    if ('statistics'%in%type) {
      # Print out quantiles argument
      print(round(out_statistics[[i]],3))
    }
    
    if ('correlations'%in%type) {
      # Print out quantiles argument
      print(round(out_cor[[i]],3))
    }
    
  }  
  
  if(object$input$n_groups==1) {
    invisible(list(gelman=out_bgr[[1]],quantiles=out_quantiles[[1]],statistics=out_statistics[[1]],correlations=out_cor[[1]]))
  } else {
    invisible(list(gelman=out_bgr,quantiles=out_quantiles,statistics=out_statistics,correlations=out_cor))
  }
  
}
