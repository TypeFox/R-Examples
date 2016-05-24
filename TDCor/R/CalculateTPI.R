CalculateTPI <-
function(dataset,l_genes,l_prior,times,time_step,N,ks_int,kd_int,delta_int,noise,delay)

{



# Checking parameters



if (N<5000)

{warning("The N parameter is too small (N<5000). The estimated TPI distributions may not be reliable.")}



if (min(ks_int)<=0)

{stop("The ks_int parameters must be positive.")}



if (min(kd_int)<=0)

{stop("The kd_int parameters must be positive.")}



if (min(delta_int)<0)

{stop("The delta_int parameters must be positive.")}



if (noise>0.25)

{warning("The noise parameter is very high. Please check.")}



if (noise<0)

{stop("The noise parameter must be positive.")}



if (length(clean.at(dataset,l_genes))!=length(l_genes))

{stop("Some genes in l_genes do not have an entry in the dataset. Please check.")}



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





# Modelling each three topologies for all the regulators in the list



prob_TPI_ind=rep(0,dim_net)

names(prob_TPI_ind)=l_genes

prob_TPI=list()

prob_TPI_domain=list()



iin=which(l_prior!=0)



for (i in iin)

{

number_of_cores=parallel::detectCores()

if  (number_of_cores==1)

{ 

mat_tpi=tpi.bs(gene_expr=rd_sub_norm[i,],ks_int,kd_int,delta_int,times,times2,time_step,noise,delay,NI=N)

}else

{



# Recruiting all the cores for the parallel computation

cl <- parallel::makeCluster(rep("localhost",number_of_cores), type = "SOCK")



# Initializing independent seeds for the random number generators of the different cores

parallel::clusterSetRNGStream(cl)



# Distributing the work between all the cores

parallel::clusterEvalQ(cl,library(TDCor))

NI=ceiling(N/number_of_cores)

mat_tpi_part=parallel::clusterCall(cl, tpi.bs, gene_expr=rd_sub_norm[i,],ks_int,kd_int,delta_int,

times,times2,time_step,noise,delay,NI=NI)



# End of the parallel processing. Merging all results into one matrix

parallel::stopCluster(cl)

mat_tpi=matrix(0,3,NI*number_of_cores)

for (n in 1:number_of_cores)

{

mat_tpi[,seq((n-1)*NI+1,n*NI)]=mat_tpi_part[[n]]

}



}





# Calculating conditional probabilities



d_cascade=density(mat_tpi[2,])

d_coreg=density(mat_tpi[3,])

d_ffl=density(mat_tpi[1,])



s_cascade=splinefun(d_cascade$x,d_cascade$y)

s_coreg=splinefun(d_coreg$x,d_coreg$y)

s_ffl=splinefun(d_ffl$x,d_ffl$y)



pin=seq(min(d_coreg$x,d_coreg$x,d_ffl$x),max(d_coreg$x,d_coreg$x,d_ffl$x),length.out=round(N/10))



p_cascade=splinefun(pin,minto0(s_cascade(pin))/(minto0(s_cascade(pin))+minto0(s_coreg(pin))+minto0(s_ffl(pin))))

p_coreg=splinefun(pin,minto0(s_coreg(pin))/(minto0(s_cascade(pin))+minto0(s_coreg(pin))+minto0(s_ffl(pin))))

p_ffl=splinefun(pin,minto0(s_ffl(pin))/(minto0(s_cascade(pin))+minto0(s_coreg(pin))+minto0(s_ffl(pin))))



prob_TPI_ind[i]=length(prob_TPI)+1

prob_TPI[[length(prob_TPI)+1]]=list(cascade=p_cascade,coreg=p_coreg,ffl=p_ffl)

prob_TPI_domain[[length(prob_TPI)+1]]=c(min(pin),max(pin))

}



input=list(l_genes=l_genes,l_prior=l_prior,times=times,time_step=time_step,N=N,ks_int=ks_int,kd_int=kd_int,delta_int=delta_int,noise=noise,delay=delay)

output=list(prob_TPI_ind=prob_TPI_ind,prob_TPI=prob_TPI,prob_TPI_domain=prob_TPI_domain,input=input)

return(output)

}
