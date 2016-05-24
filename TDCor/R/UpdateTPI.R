UpdateTPI <-
function(TPI,dataset,l_genes,l_prior)

{

uptodate=T

names(l_prior)=l_genes

l_genes_tot=unique(c(TPI$input$l_genes,l_genes))



prob_TPI=TPI$prob_TPI

prob_TPI_domain=TPI$prob_TPI_domain



prob_TPI_ind=rep(0,length(l_genes_tot))

names(prob_TPI_ind)=l_genes_tot

prob_TPI_ind[match(TPI$input$l_genes,l_genes_tot)]=TPI$prob_TPI_ind



new_genes=substract(l_genes,TPI$input$l_genes)

new_prior=l_prior[match(new_genes,l_genes)]



l_prior_tot=rep(0,length(l_genes_tot))

l_prior_tot[match(TPI$input$l_genes,l_genes_tot)]=TPI$input$l_prior

l_prior_tot[match(new_genes,l_genes_tot)]=new_prior

names(l_prior_tot)=l_genes_tot



common_genes=intersect(l_genes,TPI$input$l_genes)

if (length(common_genes)>0)

{

gene_changed_prior=common_genes[which(sign(abs(l_prior[common_genes]))>sign(abs(TPI$input$l_prior[common_genes])))]

}



if (sum(abs(new_prior))>0)

{

uptodate=F

TPI_new=CalculateTPI(dataset,new_genes,new_prior,TPI$input$times,TPI$input$time_step,TPI$input$N,TPI$input$ks_int,TPI$input$kd_int,TPI$input$delta_int,TPI$input$noise,TPI$input$delay)



kin=seq(1,length(new_genes))

for (k in kin)

{

if (new_prior[k]!=0)

{

prob_TPI_ind[match(new_genes[k],l_genes_tot)]=length(prob_TPI)+1

prob_TPI_domain[[length(prob_TPI)+1]]=TPI_new$prob_TPI_domain[[TPI_new$prob_TPI_ind[new_genes[k]]]]

prob_TPI[[length(prob_TPI)+1]]=TPI_new$prob_TPI[[TPI_new$prob_TPI_ind[new_genes[k]]]]

}

}

}



if (length(gene_changed_prior)>0)

{

uptodate=F

TPI_new=CalculateTPI(dataset,gene_changed_prior,rep(1,length(gene_changed_prior)),TPI$input$times,TPI$input$time_step,TPI$input$N,TPI$input$ks_int,TPI$input$kd_int,TPI$input$delta_int,TPI$input$noise,TPI$input$delay)

kin=seq(1,length(gene_changed_prior))

for (k in kin)

{

prob_TPI_ind[match(gene_changed_prior[k],l_genes_tot)]=length(prob_TPI)+1

prob_TPI_domain[[length(prob_TPI)+1]]=TPI_new$prob_TPI_domain[[k]]

prob_TPI[[length(prob_TPI)+1]]=TPI_new$prob_TPI[[k]]

}

}





if (uptodate)

{

message("The TPI database is up-to-date.")

}else

{

message("The TPI database has been updated.")

}





input=list(l_genes=l_genes_tot,l_prior=l_prior_tot,times=TPI$input$times,time_step=TPI$input$time_step,N=TPI$input$N,ks_int=TPI$input$ks_int,kd_int=TPI$input$kd_int,delta_int=TPI$input$delta_int,noise=TPI$input$noise,delay=TPI$input$delay)

output=list(prob_TPI_ind=prob_TPI_ind,prob_TPI=prob_TPI,prob_TPI_domain=prob_TPI_domain,input=input)

return(output)

}
