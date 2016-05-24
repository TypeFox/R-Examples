retlevel.gpd <-
function(vector,t,data,threshold)
    {
     p=1-1/t
     Nu=length(data[data>threshold])
     N=length(data)
     amostra=qgpd(1-(1/t)*N/Nu,vector[,2],threshold,vector[,1])
     res=quantile(amostra,0.5)
     res}
