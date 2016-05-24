`getmoday` <-
function(jul, iyear)
{
if(length(iyear)<length(jul))
{
iyear = rep(iyear, length(jul))

}
inine =   tojul(iyear,1,1);
ijul =    inine + jul - 1;
MD = fromjul( ijul, iyear);

return(list(mo=MD$mo, dom=MD$dom))

}

