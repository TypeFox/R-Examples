estimate_differential_splicing<-function(reads,exons,genes,do_permute)
{
  num_reads = dim(reads)[1];
  num_conditions = dim(reads)[2]-3;

  num_exons = dim(exons)[1];
  num_genes = dim(genes)[1];

  exon_per_gene = genes[,3]-genes[,2]+1;
  max_exon_per_gene = max(exon_per_gene);

  out_lr_gp = array(list(NULL),c(num_genes,max_exon_per_gene,num_conditions,num_conditions));
  out_lr_p = array(list(NULL),c(num_genes,max_exon_per_gene,num_conditions,num_conditions));

  for(i in 1:num_genes)
  {
    n = genes[i,4];
    y = matrix(0,num_conditions,n);
    exon_starts=exons[genes[i,2]:genes[i,3],2]
    exon_ends = exons[genes[i,2]:genes[i,3],3]
    n_exons = length(exon_starts);
    exon_start_map = rep(0,n_exons);
    exon_end_map = rep(0,n_exons);
    prev = 0;
    for(j in 1:n_exons)
    {
      y[1:num_conditions,(1+prev):((exon_ends[j]-exon_starts[j]+1)+prev)] = t(reads[exon_starts[j]:exon_ends[j],4:(3+num_conditions)]);  
      #Making a map of the exons in the genes
      exon_start_map[j] = 1+prev;
      exon_end_map[j] = (exon_ends[j]-exon_starts[j]+1)+prev;
      prev = (exon_ends[j]-exon_starts[j]+1)+prev;
    }
    out_gene = try(apply(y,1,generalized_poisson_likelihood),silent=TRUE);
    
    for(e in genes[i,2]:genes[i,3])
    {
      n = exons[e,4];
      z = matrix(0,num_conditions,n);
      exon_id = e-genes[i,2]+1;

      z[1:num_conditions,1:(exons[e,3]-exons[e,2]+1)] = t(reads[exons[e,2]:exons[e,3],4:(3+num_conditions)]);
     
      out_exon = try(apply(z,1,generalized_poisson_likelihood),silent=TRUE);
      
      for(r in 1:(num_conditions-1))
      {
        for(tr in (r+1):num_conditions)
        {
          if(is.list(out_gene) && is.list(out_exon))
          {
            if(out_gene[[r]]$mark == 1 && out_gene[[tr]]$mark == 1 && out_exon[[r]]$mark == 1 && out_exon[[tr]]$mark ==1)
            {
              out_lr = try(likelihood_ratio_generalized_poisson_exon_gene(z[r,],out_exon[[r]]$theta,out_exon[[r]]$lambda,y[r,],out_gene[[r]]$theta,out_gene[[r]]$lambda,z[tr,],out_exon[[tr]]$theta,out_exon[[tr]]$lambda,y[tr,],out_gene[[tr]]$theta,out_gene[[tr]]$lambda),silent=TRUE);
              out_lrp = likelihood_ratio_poisson_exon_gene(z[r,],y[r,],z[tr,],y[tr,]);
              gamma = list(shape = -1,scale = -1);
              if(is.list(out_lr))
              {
                if(out_lr$mark == 1 && do_permute == 1)
                {
                  gamma = try(perform_permutation_ds(y[r,],y[tr,],exon_start_map[exon_id],exon_end_map[exon_id],n,1000),silent=TRUE);
                  if(!is.list(gamma))
                  {
                    gamma = list(shape = -1,scale = -1);
                  }
                  cat("Mark = ",out_lr$mark,"Test = ",out_lr$Gptest,"Shape = ",gamma$shape,"Scale = ",gamma$scale,"\n");
                }
                final = c(out_lr,gamma);
                out_lr_gp[[i,e-genes[i,2]+1,r,tr]] = final;
              }
              out_lr_p[[i,e-genes[i,2]+1,r,tr]] = out_lrp;
            }
          }
        }
      }
    }
  }
  return(list(out_gp=out_lr_gp,out_p = out_lr_p));
}
