retlevel.gumbel <-
function(vector,t)
    {
     amostra=vector[,1]-vector[,2]*log(-log(1-1/t))
     res=quantile(amostra,0.5)
     res}
