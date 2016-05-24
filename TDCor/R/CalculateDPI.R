CalculateDPI <-
function(dataset,l_genes,l_prior,times,time_step,N,ks_int,kd_int,delta_int,noise,delay)

{



# Checking parameters



if (N<5000)

{warning("The N parameter is too small (N<5000). The estimated DPI distributions may not be reliable.")}



if (min(ks_int)<=0)

{stop("The ks_int parameters must be positive.")}



if (min(kd_int)<=0)

{stop("The kd_int parameters must be positive.")}



if (min(delta_int)<0)

{stop("The delta_int parameters must be positive.")}



if (noise>0.25)

{warning("The noise parameter is very high. Please check.")}



if (noise<0)

{stop("The noise parameter must be positive or equal to 0.")}



if (length(clean.at(dataset,l_genes))!=length(l_genes))

{stop("Some genes in l_genes do not have any entry in the dataset. Please check.")}



if (sum(duplicated(l_genes))>0)

{stop("Some genes are duplicated in l_genes. Please remove duplicates.")}



rd_sub=dataset[l_genes,]

dim_net=length(l_genes)



if (length(l_genes)==1)

{

rd_sub=dataset[rep(l_genes,2),]

}



# Generating the spline functions of the profiles



rd_sub_norm=norm.data(rd_sub)



SplineList=list()

times2=c(seq(times[1]-20*time_step,times[1]-time_step,time_step),times)
times3=seq(times[1]-20*time_step,times[length(times)],time_step)

for (i in 1:dim_net)

{ 

h=splinefun(times2,c(rep(rd_sub_norm[i,1],20),rd_sub_norm[i,]))

SplineList[[i]]=splinefun(times3,norm.data(h(times3)))

}





# Modelling each two topologies for all the regulators in the list



prob_DPI_ind=rep(0,dim_net)

names(prob_DPI_ind)=l_genes

prob_DPI=list()

prob_DPI_domain=list()



iin=which(l_prior!=0)



for (i in iin)

{



number_of_cores=parallel::detectCores()

if  (number_of_cores==1)

{ 

mat_dpi=dpi.bs(gene_expr=rd_sub_norm[i,],ks_int,kd_int,delta_int,times,times2,time_step,noise,delay,NI=N)

}else

{



# Recruiting all the cores for the parallel computation

cl <- parallel::makeCluster(rep("localhost",number_of_cores), type = "SOCK")



# Initializing independent seeds for the random number generators of the different cores

parallel::clusterSetRNGStream(cl)



# Distributing the work between all the cores

parallel::clusterEvalQ(cl,library(TDCor))

NI=ceiling(N/number_of_cores)

mat_dpi_part=parallel::clusterCall(cl, dpi.bs, gene_expr=rd_sub_norm[i,],ks_int,kd_int,delta_int,

times,times2,time_step,noise,delay,NI=NI)



# End of the parallel processing. Merging all results into one matrix

parallel::stopCluster(cl)

mat_dpi=matrix(0,2,NI*number_of_cores)

for (n in 1:number_of_cores)

{

mat_dpi[,seq((n-1)*NI+1,n*NI)]=mat_dpi_part[[n]]

}

}





# Calculating conditional probabilities



d_diamond=density(mat_dpi[1,])

d_coreg1=density(mat_dpi[2,])

d_coreg2=density(-mat_dpi[2,])



s_diamond=splinefun(d_diamond$x,d_diamond$y)

s_coreg1=splinefun(d_coreg1$x,d_coreg1$y)

s_coreg2=splinefun(d_coreg2$x,d_coreg2$y)



pin=seq(min(d_diamond$x,d_coreg1$x,d_coreg2$x),max(d_diamond$x,d_coreg1$x,d_coreg2$x),length.out=round(N/10))



p_diamond=splinefun(pin,minto0(s_diamond(pin))/(minto0(s_diamond(pin))+minto0(s_coreg1(pin))+minto0(s_coreg2(pin))))

p_coreg1=splinefun(pin,minto0(s_coreg1(pin))/(minto0(s_diamond(pin))+minto0(s_coreg1(pin))+minto0(s_coreg2(pin))))

p_coreg2=splinefun(pin,minto0(s_coreg2(pin))/(minto0(s_diamond(pin))+minto0(s_coreg1(pin))+minto0(s_coreg2(pin))))



prob_DPI_ind[i]=length(prob_DPI)+1

prob_DPI_domain[[length(prob_DPI)+1]]=c(min(pin),max(pin))

prob_DPI[[length(prob_DPI)+1]]=list(diamond=p_diamond,coreg1=p_coreg1,coreg2=p_coreg2)

}



input=list(l_genes=l_genes,l_prior=l_prior,times=times,time_step=time_step,N=N,ks_int=ks_int,kd_int=kd_int,delta_int=delta_int,noise=noise,delay=delay)

output=list(prob_DPI_ind=prob_DPI_ind,prob_DPI=prob_DPI,prob_DPI_domain=prob_DPI_domain,input=input)

return(output)

}
