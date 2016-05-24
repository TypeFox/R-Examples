df_TDCOR <-
function(l_genes,l_names,l_prior,rd_sub_norm,times,mat_cor,mat_delay,mat_overlap,mat_ind_wp,N,n1,SplineList,mat_isr,TPI,DPI,

thr_cor,delaymin,delaymax,thr_ind1,thr_ind2,thr_overlap,thrpTPI,thrpDPI,thr_isr,delay,time_step,search.EP,

thr_bool_EP,MinTarNumber,MinProp,MaxEPNumber)

{

## df_TDCOR for TDCOR multicore





dim_net=dim(rd_sub_norm)[1]



bootstrap=matrix(0,dim_net,dim_net)

EP=rep(0,dim_net)

names(EP)=l_names





############################

# EXTERNAL BOOOTSTRAP LOOP #

############################



for (counter in 1:N)

{

# Choice of random parameters



thr_corR=runif(1,min=thr_cor[1],max=thr_cor[2])

delayminR=runif(1,min=delaymin[1],max=delaymin[2])

delaymaxR=runif(1,min=delaymax[1],max=delaymax[2])

thr_ind1R=runif(1,min=thr_ind1[1],max=thr_ind1[2])

thr_ind2R=runif(1,min=thr_ind2[1],max=thr_ind2[2])

thr_overlapR=runif(1,min=thr_overlap[1],max=thr_overlap[2])

thrpTPIR=runif(1,min=thrpTPI[1],max=thrpTPI[2])

thrpDPIR=runif(1,min=thrpDPI[1],max=thrpDPI[2])

thr_isrR=runif(1,min=thr_isr[1],max=thr_isr[2])





##################################

## BUILDING PRELIMINARY NETWORK ##

##################################



# First condition: high correlation



net=ceiling(sign(abs(mat_cor)-thr_corR)/2)*sign(mat_cor)



# Second condition : existence of a significant delay but not too long...



filt_ID=1-ceiling(sign(mat_ind_wp-thr_ind2R)/2)

filt_ID[mat_ind_wp<thr_ind1R]=0

mat_filt1=matrix(1,dim_net,dim_net)

mat_filt1[mat_delay<=delayminR|mat_delay>delaymaxR]=0



net_shortdelay=net*mat_filt1*filt_ID



mat_filt2=matrix(1,dim_net,dim_net)

mat_filt2[mat_delay<=delaymaxR|mat_delay>2*delaymaxR]=0

mat_filt2[apply(abs(net_shortdelay),1,sum)>0,]=0

net_longdelay=net*mat_filt2





######################

## PRUNING PHASE #1 ##

######################





# prior-based pruning



prior_filt=sign(t(matrix(rep(l_prior,dim_net),dim_net,dim_net)))*net_shortdelay

prior_filt[prior_filt<0]=0

prior_filt[,l_prior==2]=1



# overlap-based pruning of negative interactions



overlap_filt=matrix(1,dim_net,dim_net)

overlap_filt[mat_overlap!=0 & mat_overlap < thr_overlapR]=0





# merging matrices



net_f=(net_shortdelay*prior_filt + net_longdelay)*overlap_filt





################################

##   INTERNAL BOOTSTRAP LOOP##

################################



net_fs=net_f

net_bs=matrix(0,dim_net,dim_net)



for (counter in 1:n1)

{

net_f=net_fs



###

# Pruning of the triangles



net_f=triangles.filter(net_f=net_f,TPI=TPI,l_genes=l_genes,thrpTPI=thrpTPIR,SplineList=SplineList,times=times,time_step=time_step)





###

# Pruning of diamonds motives



net_f=diamonds.filter(net_f=net_f,DPI=DPI,l_genes=l_genes,thrpDPI=thrpDPIR,SplineList=SplineList,times=times,time_step=time_step)



###################################################

## Looking for MRST (entry point) of the signal  ##

###################################################



if (search.EP)

{

EP_pred=entry.points(rd_sub_norm=rd_sub_norm,net_f=net_f,thr_bool_EP=thr_bool_EP,MinTarNumber=MinTarNumber,MinProp=MinProp,MaxEPNumber=MaxEPNumber)

EP=EP+EP_pred$EP/n1



# Filter based on steady state at t=0 and rapid change upon signal



if (max(EP_pred$EP)>0)

{

# SSC=0 if at steady state, SSC=1 if not (SSC Steady State comparison)



SSC=EP_pred$SSC 



# SC=0 : difference between state at first and second time point ; 



CS=EP_pred$CS



# Targets of EP must not be at their predicted steady state (with tolerance)

epin=which(EP_pred$EP>0)

for (ep in epin)

{

kin=which(net_f[,ep]!=0)

for (k in kin)

{

if (SSC[k]==0 & abs(rd_sub_norm[ep,1]-rd_sub_norm[k,1])<(1-thr_bool_EP))

{

net_f[k,ep]=0

}

if (SSC[k]==1 & CS[k]==0 & mat_delay[k,ep]<=delaymaxR)

{

net_f[k,ep]=0

}

}

}

}

}





############################################

## Looking for potential self-regulations ##

############################################





iin=which(l_prior!=0 & apply(abs(net_f),1,sum)>0)

if (length(iin)>1)

{iin=sample(iin,length(iin))}



for (i in iin)

{

kin=which(net_f[i,]!=0)

isrv=mat_isr[i,kin]



if (prod(isrv<= 1/thr_isrR)==1)

{ net_f[i,i]=1 }



if (prod(isrv>= thr_isrR)==1)

{ net_f[i,i]=-1 }

}



# Updating bootstrap



net_bs=net_bs+net_f/n1



}



bootstrap=bootstrap+net_bs



}



##################

## Finishing... ##

##################





# Returning the various matrices to the main function



output=list(net=bootstrap,EP=EP)

return(output)

}
