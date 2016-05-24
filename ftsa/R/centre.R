centre <-
function(x, type)
{
switch(type, mean = func.mean(t(x)), var = func.var(t(x)),
median = depth.FM(t(x))$median,
trimmed = depth.FM(t(x))$mtrim)
}
