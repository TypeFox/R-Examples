
matrixplot <- function(x,colorbar=FALSE,xlab="",ylab="",title=NULL,rowlabels=NULL,collabels=NULL)
{
    if(is.null(rowlabels))
    {
        rowlabels <- rownames(x)
    }

    if(is.null(collabels))
    {
        collabels <- colnames(x)
    }

    if(is.null(rowlabels))
    {
        rowlabels<-1:nrow(x)
    }

    if(is.null(collabels))
    {
        collabels<-1:ncol(x)
    }

    inds <- nrow(x):1

    image(1:ncol(x),1:nrow(x),as.matrix(t(x)[,inds]),xlab=xlab,ylab=ylab,axes=FALSE,col=hsv(0,0,1-(0:99)/99,1))

    if(!is.null(title))
    {
        title(main=title)
    }

    axis(BELOW<-1, at=1:ncol(x),labels=collabels)
    axis(LEFT<-2, at=1:nrow(x),labels=rowlabels[inds])
}   

order_round_matrix <- function(x)
{
    floor(10*as.matrix(x[,order(colSums(x),decreasing=T)]))
}

cumtime <- function()
{
    sum(
        sapply(
            strsplit(
                system(
                    'cat /proc/$PPID/stat | cut -d " " -f 14-17',
                    intern=TRUE
                ),
                ' '
            ),
            strtoi
        )
    )
}
