
reflect=function(f,...)
{
    function(x,...) ifelse(x>0, f(x,...), f(-x,...))/2
}
fold=function(f,...)
{
    function(x,...) ifelse(x>=0, f(x,...)+f(-x,...), 0)
}
