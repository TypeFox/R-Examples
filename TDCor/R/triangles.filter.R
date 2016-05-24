triangles.filter <-
function(net_f,l_genes,TPI,thrpTPI,SplineList,times,time_step)

{

sn=sum(abs(net_f))

sn0=0

delay=TPI$input$delay



while(sn0!=sn)

{sn0=sn



iin=which(rowSums(abs(net_f))>1)

if (length(iin)>1)

{iin=sample(iin,length(iin))}



for (i in iin)

{ kin=which(net_f[i,]!=0)

if (length(kin)>1)

  {kin=sample(kin,length(kin))}

  for (k in kin)

  {

creg<-which(net_f[k,]*net_f[i,]!=0)

if (length(creg)>1)

   {

creg=sample(creg,length(creg))

}



    for (u in creg)

    { 

predictor=tpi.index(SplineList[[u]],SplineList[[k]],SplineList[[i]],time_l=times[1]+delay,time_u=times[length(times)],time_step=time_step,delay=delay)

 

if ( TPI$prob_TPI[[TPI$prob_TPI_ind[l_genes[u]]]][[1]](predictor) >= thrpTPI) # pruning cascade errors

     {

net_f[i,u]=0

}

if  (TPI$prob_TPI[[TPI$prob_TPI_ind[l_genes[u]]]][[2]](predictor) >= thrpTPI) # pruning fan-out errors

      {

net_f[i,k]=0

}

    }



}

}

sn=sum(abs(net_f))

}



return(net_f)

}
