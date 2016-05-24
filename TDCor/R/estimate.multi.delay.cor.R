estimate.multi.delay.cor <-
function(SplineList,l_prior,delayspan,delaymax,times,time_step,thr_cor,tol)

{

dim_net=length(l_prior)

iin=1:dim_net

jin=which(l_prior!=0)



inter_time=seq(times[1],times[length(times)],time_step) # vector of time points analysed

delayin=seq(-delayspan,delayspan,time_step)# vector of delays to be screened

mat_cor=matrix(0,dim_net,dim_net)# mat_cor will contain the final measures of correlation (output by the function)

mat_delay=matrix(0,dim_net,dim_net)# mat_delay will contain the final measures of delay (output by the function)



number_of_cores=parallel::detectCores()

if  (number_of_cores==1)

{ 

for (j in jin)

{

tested_targets=substract(iin,jin[jin<=j])

if (length(tested_targets)>0)

{

cad=sapply(SplineList[tested_targets],estimate.delay.cor,hr=SplineList[[j]], inter_time=inter_time,times=times,time_step=time_step,thr_cor=thr_cor,tol=tol,delaymax=delaymax,delayin=delayin)


cad[c(2,4),l_prior[tested_targets]==0]=0

mat_cor[tested_targets,j]=cad[1,]

mat_cor[j,tested_targets]=cad[2,]

mat_delay[tested_targets,j]=cad[3,]

mat_delay[j,tested_targets]=cad[4,]

}

}

}else

{



# Recruiting all the cores for the parallel computation

cl <- parallel::makeCluster(rep("localhost",number_of_cores), type = "SOCK")



# Distributing the work between all the cores

parallel::clusterEvalQ(cl,library(TDCor))



for (j in jin)

{

tested_targets=substract(iin,jin[jin<=j])

if (length(tested_targets)>0)

{

cad=parallel::parSapplyLB(cl,SplineList[tested_targets],estimate.delay.cor,hr=SplineList[[j]], inter_time=inter_time,times=times,

time_step=time_step,thr_cor=thr_cor,tol=tol,delaymax=delaymax,delayin=delayin)



cad[c(2,4),l_prior[tested_targets]==0]=0

mat_cor[tested_targets,j]=cad[1,]

mat_cor[j,tested_targets]=cad[2,]

mat_delay[tested_targets,j]=cad[3,]

mat_delay[j,tested_targets]=cad[4,]

}

}

# Stop cluster

parallel::stopCluster(cl)

}



return(list(mat_cor=mat_cor,mat_delay=mat_delay))

}
