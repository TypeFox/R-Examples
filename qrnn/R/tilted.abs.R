tilted.abs <-
function(x, tau)
{
    tabs <- x
    tabs[x>=0] <- tabs[x>=0]*tau
    tabs[x<0] <- tabs[x<0]*(tau-1)
    tabs
}

