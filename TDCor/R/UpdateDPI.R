UpdateDPI <-
function(DPI,dataset,l_genes,l_prior)

{

uptodate=T

names(l_prior)=l_genes

l_genes_tot=unique(c(DPI$input$l_genes,l_genes))



prob_DPI=DPI$prob_DPI

prob_DPI_domain=DPI$prob_DPI_domain



prob_DPI_ind=rep(0,length(l_genes_tot))

names(prob_DPI_ind)=l_genes_tot

prob_DPI_ind[match(DPI$input$l_genes,l_genes_tot)]=DPI$prob_DPI_ind



new_genes=substract(l_genes,DPI$input$l_genes)

new_prior=l_prior[match(new_genes,l_genes)]



l_prior_tot=rep(0,length(l_genes_tot))

l_prior_tot[match(DPI$input$l_genes,l_genes_tot)]=DPI$input$l_prior

l_prior_tot[match(new_genes,l_genes_tot)]=new_prior

names(l_prior_tot)=l_genes_tot



common_genes=intersect(l_genes,DPI$input$l_genes)

if (length(common_genes)>0)

{

gene_changed_prior=common_genes[which(sign(abs(l_prior[common_genes]))>sign(abs(DPI$input$l_prior[common_genes])))]

}





if (sum(abs(new_prior))>0)

{

uptodate=F

DPI_new=CalculateDPI(dataset,new_genes,new_prior,DPI$input$times,DPI$input$time_step,DPI$input$N,DPI$input$ks_int,DPI$input$kd_int,DPI$input$delta_int,DPI$input$noise,DPI$input$delay)



kin=seq(1,length(new_genes))

for (k in kin)

{

if (new_prior[k]!=0)

{

prob_DPI_ind[match(new_genes[k],l_genes_tot)]=length(prob_DPI)+1

prob_DPI_domain[[length(prob_DPI)+1]]=DPI_new$prob_DPI_domain[[DPI_new$prob_DPI_ind[new_genes[k]]]]

prob_DPI[[length(prob_DPI)+1]]=DPI_new$prob_DPI[[DPI_new$prob_DPI_ind[new_genes[k]]]]

}

}

}



if (length(gene_changed_prior)>0)

{

uptodate=F

DPI_new=CalculateDPI(dataset,gene_changed_prior,rep(1,length(gene_changed_prior)),DPI$input$times,DPI$input$time_step,DPI$input$N,DPI$input$ks_int,DPI$input$kd_int,DPI$input$delta_int,DPI$input$noise,DPI$input$delay)



kin=seq(1,length(gene_changed_prior))

for (k in kin)

{

prob_DPI_ind[match(gene_changed_prior[k],l_genes_tot)]=length(prob_DPI)+1

prob_DPI_domain[[length(prob_DPI)+1]]=DPI_new$prob_DPI_domain[[k]]

prob_DPI[[length(prob_DPI)+1]]=DPI_new$prob_DPI[[k]]

}



}



if (uptodate)

{

message("The DPI database is up-to-date.")

}else

{

message("The DPI database has been updated.")

}



input=list(l_genes=l_genes_tot,l_prior=l_prior_tot,times=DPI$input$times,time_step=DPI$input$time_step,N=DPI$input$N,ks_int=DPI$input$ks_int,kd_int=DPI$input$kd_int,delta_int=DPI$input$delta_int,noise=DPI$input$noise,delay=DPI$input$delay)

output=list(prob_DPI_ind=prob_DPI_ind,prob_DPI=prob_DPI,prob_DPI_domain=prob_DPI_domain,input=input)

return(output)

}
