retlevel.gev <-
function(vector,t)
    {
     amostra=qgev(1-1/t,vector[,3],vector[,1],vector[,2])
     res=quantile(amostra,0.5)
     res}
