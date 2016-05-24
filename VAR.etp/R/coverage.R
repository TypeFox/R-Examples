coverage <-
function(intv,truehl)
{
cover <- 0
if(intv[1] < truehl & intv[2] > truehl) cover <- 1
return(cover)
}
