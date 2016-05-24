local.min <-
function(x)

{

xl=length(x)



local_opt=which(sign(x[2:(xl-1)]-x[1:(xl-2)])==sign(x[2:(xl-1)]-x[3:xl]))+1

local_min=intersect(local_opt,which(x[2:(xl-1)]-x[1:(xl-2)]<0)+1)



return(local_min)

}
