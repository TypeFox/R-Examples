posterior.gpd.exp <-
function(data,threshold,int)
  {burnin=thin=10;int*thin/2  
   ajuste=gpd(data,threshold);data=ajuste$data;n=length(data)
   taumcb=rgamma(int,n+0.001,0.001+n*mean(data-threshold))
   sigmamcb=1/taumcb
   sigmamcb}
