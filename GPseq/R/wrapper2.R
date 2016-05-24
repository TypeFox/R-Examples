estimate_differential_expression<-function(reads,exons,genes,norm_gp,norm_p,do_permute)
{
  num_reads = dim(reads)[1];
  num_conditions = dim(reads)[2]-3;

  num_exons = dim(exons)[1];
  num_genes = dim(genes)[1];

  compare_gp = array(list(NULL),c(num_genes,num_conditions,num_conditions));
  compare_p = array(list(NULL),c(num_genes,num_conditions,num_conditions));

  for(i in 1:num_genes)
  {
    n = genes[i,4];
    y = matrix(0,num_conditions,n);
    exon_starts=exons[genes[i,2]:genes[i,3],2]
    exon_ends = exons[genes[i,2]:genes[i,3],3]
    n_exons = length(exon_starts);
    prev = 0;
  
    for(j in 1:n_exons)
    {
      y[1:num_conditions,(1+prev):((exon_ends[j]-exon_starts[j]+1)+prev)] = t(reads[exon_starts[j]:exon_ends[j],4:(3+num_conditions)]);  
      prev = (exon_ends[j]-exon_starts[j]+1)+prev;
    }
    out = try(apply(y,1,generalized_poisson_likelihood),silent=TRUE);

#Comparing across all conditions
    for(j in 1:(num_conditions-1))
    {
      for(k in (j+1):num_conditions)
      {
        if(is.list(out))
        {
          if(out[[j]]$mark==1 && out[[k]]$mark==1)
          {
            norm = norm_gp[k]/norm_gp[j];
            compare_j_k = try(likelihood_ratio_tissue_generalized_poisson(y[j,],out[[j]]$lambda,out[[j]]$theta,y[k,],out[[k]]$lambda,out[[k]]$theta,norm),silent=TRUE);
            if(do_permute == 0)
            {
              gamma =list(shape = -1,scale=-1);
            }
            else
            {
              gamma = try(perform_permutation_de(y[j,],y[k,],1000),silent=TRUE);
              if(!is.list(gamma))
              {
                gamma = list(shape=-1,scale=-1);
              }
            }
            if(is.list(compare_j_k))
            {
              final = c(compare_j_k,gamma);
              compare_gp[[i,j,k]] = final;
            }
          }
        }
        normp = norm_p[k]/norm_p[j];
        compare_j_k_poisson = likelihood_ratio_tissue_poisson(y[j,],mean(y[j,]),y[k,],mean(y[k,]),normp);
        compare_p[[i,j,k]] = compare_j_k_poisson;
      }
    }
  }
  return(list(gp_comparison = compare_gp,p_comparison = compare_p));
}

    
