estimate_exon_gene_expression<-function(reads,exons,genes)
{
  num_reads = dim(reads)[1];
  num_conditions = dim(reads)[2]-3;

  num_exons = dim(exons)[1];
  num_genes = dim(genes)[1];

  exon_out = array(list(NULL),c(num_exons,num_conditions));
  chisq_exon_out = array(list(NULL),c(num_exons,num_conditions));

  gene_out = array(list(NULL),c(num_genes,num_conditions));
  chisq_gene_out = array(list(NULL),c(num_genes,num_conditions));

  norm_gp = rep(0,num_conditions);
  norm_p = rep(0,num_conditions);

  for(i in 1:num_genes)
  {
#Calculating Parameters and Goodness of Fit Statistics for exons
    for(e in genes[i,2]:genes[i,3])
    {
      n = exons[e,4];
      y = matrix(0,num_conditions,n);

      y[1:num_conditions,1:(exons[e,3]-exons[e,2]+1)] = t(reads[exons[e,2]:exons[e,3],4:(3+num_conditions)]);
      
      out = try(apply(y,1,generalized_poisson_likelihood),silent=TRUE);
      for(j in 1:num_conditions)
      {
        if(is.list(out))
        {
          if(out[[j]]$mark==1)
          {
            chisq=try(calc_chisq_statistic(y[j,],out[[j]]$lambda,out[[j]]$theta),silent=TRUE); 
            if(is.list(chisq))
            {
              chisq_exon_out[[e,j]] = chisq;
            }
          }
          exon_out[[e,j]] = out[[j]];
        }
        else
        {
          exon_out[[e,j]] = list(mark=0,theta=-1,lambda=-1,y_bar=mean(y[j,]),length=n);
        }
      }
    }
    
#Calculating Parameters and Goodness of Fit Statistics for genes

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
    for(j in 1:num_conditions)
    {
      if(is.list(out))
      {
        if(out[[j]]$mark==1)
        {
          chisq=try(calc_chisq_statistic(y[j,],out[[j]]$lambda,out[[j]]$theta),silent=TRUE); 
          norm_gp[j] = norm_gp[j]+ (out[[j]]$theta*out[[j]]$length);
          if(is.list(chisq))
          {
            chisq_gene_out[[i,j]] = chisq;
          }
        }
        norm_p[j] = norm_p[j] + ((out[[j]]$y_bar)*(out[[j]]$length));
        gene_out[[i,j]] = out[[j]];
      }
      else
      {
        gene_out[[i,j]] = list(mark=0,theta=-1,lambda=-1,y_bar=mean(y[j,]),length=n);
      }
    }
  }
    
 
  return(list(exon_out = exon_out,chisq_exon_out = chisq_exon_out,gene_out = gene_out,chisq_gene_out = chisq_gene_out,norm_gp = norm_gp, norm_p = norm_p));
}
