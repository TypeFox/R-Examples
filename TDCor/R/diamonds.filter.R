diamonds.filter <-
function(net_f,DPI,l_genes,thrpDPI,SplineList,times,time_step)

{

sn=sum(abs(net_f))

sn0=0

delay=DPI$input$delay



while(sn0!=sn)

{sn0=sn



iin=which(rowSums(net_f)>1)

if (length(iin)>1)

{ iin=sample(iin,length(iin))}



for (i in iin)

{ 

kin=which(net_f[i,]!=0)

if (length(kin)>1)

  {

kin=sample(kin,length(kin))

}

  for (k1 in kin)

 { 

for (k2 in substract(kin,k1))

{ 

crt=which(net_f[k1,]*net_f[k2,]!=0)

  if (length(crt)>0)

  { 

u=rep(0,length(crt))

v=rep(0,length(crt))



predictor=dpi.index(hr1=SplineList[[k1]],hr2=SplineList[[k2]],ht=SplineList[[i]],time_l=times[1]+delay,time_u=times[length(times)],time_step=time_step,delay=delay)

for (k3 in  crt)

{ 

u[which(crt==k3)]=DPI$prob_DPI[[DPI$prob_DPI_ind[l_genes[k3]]]][[2]](predictor)

v[which(crt==k3)]=DPI$prob_DPI[[DPI$prob_DPI_ind[l_genes[k3]]]][[3]](predictor)

}



          if (sum(u>thrpDPI)> length(crt)/2) 

      {

net_f[i,k2]=0

}

          if (sum(v>thrpDPI)> length(crt)/2) 

      {

net_f[i,k1]=0

}

         }

      }

  }  

}



sn=sum(abs(net_f))

}

return(net_f)

}
