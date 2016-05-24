"simpExp" <-
function (sumT, par, t, tc) 
{
    sumexp <- exp(-par[2] * (t - tc))
    suma <- 1 
    if (length(par) > 2) { 
       for(i in seq(3, length(par), by = 2)){
	     sumexp <- sumexp + (par[i] * exp(-par[i+1] * (t - tc)))
	     suma <- suma + par[i]
       }     
        
    }
    resT <- sumT + (par[1] * (sumexp/suma))
    resT
}
