TDCOR <-
function(dataset,l_genes,TPI,DPI,...)

{

## TDCOR v.1.2



TDCor_Version="v.1.2"





# check for additional function arguments



    Lst <- list(...)



    if( !is.null(Lst$tol) )

    { tol <- Lst$tol}else

    { tol <- 0.15} # tolerance used by the function potential.interactions()



    if( !is.null(Lst$delayspan) )

    { delayspan<- Lst$delayspan}else

    { delayspan <- 12} # span of delays that will be screed in  potential.interactions()



    if( !is.null(Lst$thr_cor) )

    { thr_cor <- Lst$thr_cor

if (length(thr_cor)==1){thr_cor=rep(thr_cor,2)}

    }else

    { thr_cor <-c(0.7,0.9)} # threshold for correlation 



    if( !is.null(Lst$delaymax) )

    { delaymax <- Lst$delaymax

if (length(delaymax)==1){delaymax=rep(delaymax,2)}

    }else

    { delaymax <- c(3,3)} # threshold for time delay maximum for a direct interaction (in hours) -included



    if( !is.null(Lst$delaymin) )

    { delaymin <- Lst$delaymin

if (length(delaymin)==1){delaymin=rep(delaymin,2)}

    }else

    { delaymin <-c(0,0)} # threshold for time delay minimum for a direct interaction (in hours) -excluded



    if( !is.null(Lst$delay) )

    { delay <- Lst$delay}else

    { delay <-  3} # delay used to calculate the index of directness among other things



    if( !is.null(Lst$thr_overlap) )

    { thr_overlap <- Lst$thr_overlap

if (length(thr_overlap)==1){thr_overlap=rep(thr_overlap,2)}

    }else

    { thr_overlap <- c(0.5,0.6)} # threshold for overlap (negative interaction pruning)



     if( !is.null(Lst$thrpTPI) )

    { thrpTPI<- Lst$thrpTPI

if (length(thrpTPI)==1){thrpTPI=rep(thrpTPI,2)}

    }else

    { thrpTPI <-c(0.5,0.75)}    # threshold for pruning cascade errors in df_TDCOR



    if( !is.null(Lst$thrpDPI) )

    { thrpDPI<- Lst$thrpDPI

if (length(thrpDPI)==1){thrpDPI=rep(thrpDPI,2)}

    }else

    { thrpDPI<-c(0.8,0.9)}    # threshold for pruning order 2 cascade errors in df_TDCOR



    if( !is.null(Lst$thr_isr) )

    { thr_isr <- Lst$thr_isr

      if (length(thr_isr)==1){thr_isr=rep(thr_isr,2)}

    }else

    { thr_isr <- c(3,6)}# threshold for self-regulation 



    if( !is.null(Lst$times) )

    { times <- Lst$times}else

    { times <- c(0,6+seq(0,16)*3)}# sequence of the times of sampling needed for self-regulations analysis



    if( !is.null(Lst$thr_ind1) )

    { thr_ind1 <- Lst$thr_ind1

      if (length(thr_ind1)==1){thr_ind1=rep(thr_ind1,2)}

    }else

    { thr_ind1 <- c(0.5,0.5)}# threshold for pruning using index of directness



    if( !is.null(Lst$thr_ind2) )

    { thr_ind2 <- Lst$thr_ind2

      if (length(thr_ind2)==1){thr_ind2=rep(thr_ind2,2)}

    }else

    { thr_ind2 <- c(4.5,4.5)}# threshold for pruning using index of directness



    if( !is.null(Lst$time_step) )

    { time_step <- Lst$time_step}else

    { time_step <- 1}# time step in time unit for analysing the fan-in errors



    if( !is.null(Lst$search.EP) )

    { search.EP <- Lst$search.EP}else

    { search.EP <- TRUE}# activate or not the search for the signal "entry point" (EP) in the GRN



    if( !is.null(Lst$thr_bool_EP) )

    { thr_bool_EP <- Lst$thr_bool_EP}else

    { thr_bool_EP <- 0.8}# threshold for converting normalised into boolean profile during analysis of steady states for predicting EP



    if( !is.null(Lst$MinTarNumber) )

    { MinTarNumber <- Lst$MinTarNumber}else

    { MinTarNumber <- 5}# minimum number of targets that are not at steady state at t=0 for the regulator to be listed as potential EP



    if( !is.null(Lst$MinProp) )

    { MinProp <- Lst$MinProp}else

    { MinProp <- 0.75}   # minimum proportion of targets in [0,1] that must no be at steady state at t=0 for the regulator to be listed as potential EP



    if( !is.null(Lst$MaxEPNumber) )

    { MaxEPNumber <- Lst$MaxEPNumber}else

    { MaxEPNumber <- 1}   # Maximum number of possible EP



    if( !is.null(Lst$regmax) )

    { regmax <- Lst$regmax}else

    { regmax <- 6}# maximum number of regulator allowed for each targets



    if( !is.null(Lst$n0) )

    { n0 <- Lst$n0}else

    { n0 <-1000}# number of iteration for bootstraping df_TDCOR (external loop)



    if( !is.null(Lst$n1) )

    { n1 <- Lst$n1}else

    { n1 <-10}# number of iteration for bootstraping inside df_TDCOR (interal loop)



    if( !is.null(Lst$l_names) )

    { l_names <- Lst$l_names}else

    { l_names <-l_genes}# names of the genes to be returned



    if( !is.null(Lst$l_prior) )

    { l_prior <- Lst$l_prior}else

    { l_prior=rep(2,length(l_genes))} # prior information about the nature of the regulators



    if( !is.null(Lst$outfile_name) )

    { outfile_name <- Lst$outfile_name}else

    { outfile_name="TDCor_output.txt"} # Name of the text file to print the predictions in.



# Checking parameters



if (min(thr_cor)<0.6)

{warning("The thr_cor parameter is very low. Please check.")}



if (min(thrpTPI)<0.5)

{stop("The thrpTPI parameter is too low. The minimum should be greater than or equal to 0.5.")}



if (min(thrpDPI)<0.5)

{stop("The thrpDPI parameter is too low. The minimum should be greater than or equal to 0.5.")}



if (max(thrpTPI)>1)

{stop("The thrpTPI parameter is too high. The maximum should be smaller than or equal to 1.")}



if (max(thrpDPI)>1)

{stop("The thrpDPI parameter is too high. The maximum should be smaller than or equal to 1.")}



if (length(substract(levels(factor(l_prior)),c("-1","0","1","2")))>0)

{stop("The prior contains impossible values (It should contains only -1, 0, 1 and 2). Please check.")}



if (regmax<3)

{warning("regmax is very small. Please check.")}



if (regmax<1)

{stop("regmax should be greater than or equal to 1.")}



if (length(clean.at(dataset,l_genes))!=length(l_genes))

{stop("Some genes in l_genes do not have an entry in the dataset. Please check or use clean.at() to remove them.")}



if (sum(duplicated(l_genes))>0)

{stop("Some genes are duplicated in l_genes. Please remove duplicates.")}



# Starts here !



rd_sub=dataset[l_genes,]

rownames(rd_sub)=l_names

dim_net=dim(rd_sub)[1]

dim_time=dim(rd_sub)[2]





##

message("Generating the interpolating spline functions...")

##



rd_sub_norm=norm.data(rd_sub)



SplineList=list()

times2=c(seq(times[1]-20*time_step,times[1]-time_step,time_step),times)
times3=seq(times[1]-20*time_step,times[length(times)],time_step)

for (i in 1:dim_net)

{ 

h=splinefun(times2,c(rep(rd_sub_norm[i,1],20),rd_sub_norm[i,]))
SplineList[[i]]=splinefun(times3,norm.data(h(times3)))

}





### Calculating the matrices that will feed df_TDCOR





##

message("Computing of the delays and associated correlations...")

##



pot_int=estimate.multi.delay.cor(SplineList=SplineList,l_prior=l_prior,delayspan=delayspan,delaymax=max(delaymax),times=times,time_step=time_step,thr_cor=thr_cor,tol=tol)

mat_cor=pot_int$mat_cor

mat_delay=pot_int$mat_delay



##

message("Computing the indices of directness and the indices of self-regulation...")

##



mat_ind_wp=matrix(0,dim_net,dim_net)

mat_isr=matrix(0,dim_net,dim_net)





iin=which(apply(mat_delay,1,sum)>0)

for (i in iin)

{



jin=which(mat_cor[i,]>0)

if (length(jin)>0)

{

mat_ind_wp[i,jin]=round(sapply(SplineList[jin],directness.index,

ht=SplineList[[i]],time_l=times[1], time_u=times[length(times)],

delay=delay,type=1),3)

}



jin=which(mat_cor[i,]<0)

if (length(jin)>0)

{

mat_ind_wp[i,jin]=round(sapply(SplineList[jin],directness.index,

ht=SplineList[[i]],time_l=times[1], time_u=times[length(times)],

delay=delay,type=-1),3)

}





if (l_prior[i]!=0)

{

jin=which(mat_cor[i,]!=0)

mat_isr[i,jin]=round(sapply(SplineList[jin],isr.index,

ht=SplineList[[i]],time_l=times[1],time_u=times[length(times)],

delay=delay),3)

}

}



##

message("Computing the indices of overlap...")

##



mat_overlap=round(overlap.index(net_pot=sign(mat_cor),SplineList=SplineList,times=times,time_step=time_step,mat_delay=mat_delay),3)





#####################

##  Bootstrapping  ##

#####################



##

message("Bootstrapping...")

##





number_of_cores=parallel::detectCores()

if  (number_of_cores==1)

{ 

df_TDCOR_output=df_TDCOR(l_genes=l_genes,l_names=l_names,l_prior=l_prior,rd_sub_norm=rd_sub_norm,times=times,mat_cor=mat_cor,mat_delay=mat_delay,mat_overlap=mat_overlap,delaymin=delaymin,delaymax=delaymax,

mat_ind_wp=mat_ind_wp,thr_cor=thr_cor,thr_overlap=thr_overlap,thrpTPI=thrpTPI,thrpDPI=thrpDPI,TPI=TPI,DPI=DPI,

N=n0,n1=n1,SplineList=SplineList,time_step=time_step,thr_ind1=thr_ind1,thr_ind2=thr_ind2,thr_bool_EP=thr_bool_EP,

MinTarNumber=MinTarNumber,MinProp=MinProp,MaxEPNumber=MaxEPNumber,search.EP=search.EP,thr_isr=thr_isr,mat_isr=mat_isr)



bootstrap=data.frame(df_TDCOR_output$net)
names(bootstrap)=l_names
rownames(bootstrap)=l_names
bootstrap=as.matrix(bootstrap)


EP=df_TDCOR_output$EP
names(EP)=l_names


}else

{



# Recruiting all the cores for the parallel computation

cl <- parallel::makeCluster(rep("localhost",number_of_cores), type = "SOCK")


# Initializing independent seeds for the random number generators of the different cores

parallel::clusterSetRNGStream(cl)


# Distributing the work between all the cores

parallel::clusterEvalQ(cl,library(TDCor))



df_TDCOR_output_parts=parallel::clusterCall(cl, df_TDCOR,
l_genes=l_genes,l_names=l_names,l_prior=l_prior,rd_sub_norm=rd_sub_norm,times=times,
mat_cor=mat_cor,mat_delay=mat_delay,mat_overlap=mat_overlap,delaymin=delaymin,delaymax=delaymax,
mat_ind_wp=mat_ind_wp,thr_cor=thr_cor,thr_overlap=thr_overlap,thrpTPI=thrpTPI,thrpDPI=thrpDPI,
TPI=TPI,DPI=DPI,N=ceiling(n0/number_of_cores),n1=n1,SplineList=SplineList,time_step=time_step,
thr_ind1=thr_ind1,thr_ind2=thr_ind2,thr_bool_EP=thr_bool_EP,MinTarNumber=MinTarNumber,
MinProp=MinProp,MaxEPNumber=MaxEPNumber,search.EP=search.EP,thr_isr=thr_isr,mat_isr=mat_isr)



# End of the parallel processing. Merging all results



parallel::stopCluster(cl)

bootstrap=df_TDCOR_output_parts[[1]]$net

EP=df_TDCOR_output_parts[[1]]$EP



for (n in 2:number_of_cores)

{

bootstrap=bootstrap+df_TDCOR_output_parts[[n]]$net

EP=EP+df_TDCOR_output_parts[[n]]$EP

}

bootstrap=data.frame(bootstrap)
names(bootstrap)=l_names
rownames(bootstrap)=l_names
bootstrap=as.matrix(bootstrap)

names(EP)=l_names

n0=ceiling(n0/number_of_cores)*number_of_cores

}







#######################################

## PRUNING  max number of regulators ##

#######################################



##

message(paste("Selecting the ", regmax," best regulators...",sep=""))

##



self_reg=diag(bootstrap)

diag(bootstrap)=0

iin=which(rowSums(abs(bootstrap))>regmax)



for (i in iin)

{ 

 k=rev(sort(abs(bootstrap[i,])))

 o=which(abs(bootstrap[i,])<k[regmax])

   bootstrap[i,o]=0

}



diag(bootstrap)=self_reg



###############

## Finishing ##

###############





mat_cor=data.frame(mat_cor) ; rownames(mat_cor)=l_names ; names(mat_cor)=l_names ;mat_cor=as.matrix(mat_cor)

mat_delay=data.frame(mat_delay) ; rownames(mat_delay)=l_names ; names(mat_delay)=l_names ;mat_delay=as.matrix(mat_delay)

mat_isr=data.frame(mat_isr) ; rownames(mat_isr)=l_names ; names(mat_isr)=l_names ;mat_isr=as.matrix(mat_isr)

mat_overlap=data.frame(mat_overlap) ; rownames(mat_overlap)=l_names ; names(mat_overlap)=l_names ;mat_overlap=as.matrix(mat_overlap)



intermediate=list(mat_cor=mat_cor,mat_isr=mat_isr,mat_overlap=mat_overlap)



input=list(version=TDCor_Version,l_genes=l_genes,l_prior=l_prior,tol=tol,thr_cor=thr_cor,thr_ind1=thr_ind1,thr_ind2=thr_ind2,thr_overlap=thr_overlap,n0=n0,n1=n1,delaymin=delaymin,delaymax=delaymax,delayspan=delayspan,delay=delay,

thrpDPI=thrpDPI,thrpTPI=thrpTPI,thr_bool_EP=thr_bool_EP,MinTarNumber=MinTarNumber,MinProp=MinProp,MaxEPNumber=MaxEPNumber,thr_isr=thr_isr,regmax=regmax,TPI_input=TPI$input,DPI_input=DPI$input)



output=list(input,intermediate,bootstrap,mat_ind_wp,mat_delay,EP)

names(output)=c("input","intermediate","network","ID","delay","EP")



tbl=table.nw(output,outfile=outfile_name,n=n0)



print(tbl)

message("Done !")

message(paste("TDCor predictions have been printed in the '", outfile_name,"' file.",sep=""))

message(paste("TDCor input parameters have been printed in the 'Param-", outfile_name,"' file.",sep=""))



output[[7]]=tbl

names(output)[7]="predictions"



return(output)

}
