get.LODallshare <- function(vec,pshare)
{
if (any(pshare$ped.tocompute.vec%in%vec)) sum(pshare$mlog10pshare[pshare$ped.tocompute.vec%in%vec])
else NA
}
