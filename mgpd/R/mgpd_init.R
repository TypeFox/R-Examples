mgpd_init <-
function(xdat)
{
p              = NULL
for(i in 1:dim(xdat)[2]) p = c(p,fgev(xdat[,i],std.err = FALSE)[[1]])
p
}
