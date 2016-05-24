"count.points.id" <- function(xy, id, w)
{
    ## Verifications
    if (ncol(xy)!=2)
        stop("xy should have 2 columns")
    if (length(id) != nrow(xy))
        stop("id should be of the same length as xy")

    ## Prepares the data
    x<-xy[,1]
    y<-xy[,2]
    id<-factor(id)
    lx<-split(x, id)
    ly<-split(y, id)
    output<-list()

    ## Use of the function count.points for each animal
    for (i in 1:length(levels(id)))
        output[[levels(id)[i]]]<-count.points(cbind(lx[[i]], ly[[i]]), w)

    ## output
    output<-as.kasc(output)
}

