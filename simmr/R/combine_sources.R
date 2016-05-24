combine_sources = function(simmr_out,to_combine=simmr_out$input$source_names[1:2],new_source_name='combined_source') {
  # A posteriori combining of sources
  
  # Check only two sources to be combined
  if(length(to_combine)!=2) stop("Currently only two sources can be combined")
  
  # Check class
  if(class(simmr_out)!='simmr_output') stop("Only objects of class simmr_output can be run through this function")

  # Find which columns to combine by number 
  to_combine_cols = match(to_combine,simmr_out$input$source_names)
  if(any(is.na(to_combine_cols))) stop('1 or more source names not found')

  simmr_new_out = simmr_out
  
  # 1 combine the chosen source means
  old_source_means = simmr_out$input$source_means 
  simmr_new_out$input$source_means = old_source_means[-to_combine_cols[2],]
  simmr_new_out$input$source_means[to_combine_cols[1],] = apply(old_source_means[to_combine_cols,],2,'mean')
  
  # 2 combine the source sds
  old_source_sds = simmr_out$input$source_sds
  simmr_new_out$input$source_sds = old_source_sds[-to_combine_cols[2],]
  simmr_new_out$input$source_sds[to_combine_cols[1],] = apply(old_source_sds[to_combine_cols,],2,function(x) sqrt(sum(x^2)))
  
  # 3 combine the correction means
  old_correction_means = simmr_out$input$correction_means 
  simmr_new_out$input$correction_means = old_correction_means[-to_combine_cols[2],]
  simmr_new_out$input$correction_means[to_combine_cols[1],] = apply(old_correction_means[to_combine_cols,],2,'mean')

  # 4 combine the correction sds
  old_correction_sds = simmr_out$input$correction_sds
  simmr_new_out$input$correction_sds = old_correction_sds[-to_combine_cols[2],]
  simmr_new_out$input$correction_sds[to_combine_cols[1],] = apply(old_correction_sds[to_combine_cols,],2,function(x) sqrt(sum(x^2)))
  
  # 5 combine the concentraion means
  old_concentration_means = simmr_out$input$concentration_means 
  simmr_new_out$input$concentration_means = old_concentration_means[-to_combine_cols[2],]
  simmr_new_out$input$concentration_means[to_combine_cols[1],] = apply(old_concentration_means[to_combine_cols,],2,'mean')
  
  # 6 change the source names
  old_source_names = simmr_out$input$source_names
  simmr_new_out$input$source_names = old_source_names[-to_combine_cols[2]]
  simmr_new_out$input$source_names[to_combine_cols[1]] = new_source_name
  
  # 7 Change n_sources
  simmr_new_out$input$n_sources = simmr_new_out$input$n_sources - 1
  
  # 8 Sum across all the output values
  old_pars = simmr_out$output
  for(j in 1:length(old_pars)) {
    for(i in 1:length(old_pars[[j]])) {
      simmr_new_out$output[[j]][[i]] = old_pars[[j]][[i]][,-to_combine_cols[2]]
      sum_post = old_pars[[j]][[i]][,to_combine_cols[1]]+old_pars[[j]][[i]][,to_combine_cols[2]]
      simmr_new_out$output[[j]][[i]][,to_combine_cols[1]] = sum_post
      colnames(simmr_new_out$output[[j]][[i]]) = c(simmr_new_out$input$source_names,colnames(simmr_new_out$input$mixtures))
    }  
  }
  
  class(simmr_new_out) = 'simmr_output'
  return(simmr_new_out)
  
}