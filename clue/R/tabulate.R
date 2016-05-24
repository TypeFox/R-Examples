cl_tabulate <-
function(x)
{
    values <- unique(x)
    counts <- tabulate(match(x, values))
    ## Still a bit tricky to create a data frame with a list "column"
    ## which is not protected by I(); otherwise, we oculd simply do
    ##   data.frame(values = I(values), counts = counts)
    out <- data.frame(values = double(length(values)), counts = counts)
    out$values <- values
    out
}
