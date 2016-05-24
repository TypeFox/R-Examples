`Blocksr` <-
function (x) 
{
y<-Blocks(x)

y <- c(y, rev(y))
        n <- length(y)
    y[seq(from=1, to=n, by=2)]





}

