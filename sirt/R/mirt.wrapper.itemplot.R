


mirt.wrapper.itemplot <- function( mirt.obj , ask=TRUE , ...){
    I <- ncol( mirt.obj@Data$data )
    for (ii in 1:I){
        # ii <- 1
        main <- paste0("Trace Lines of Item " , colnames( mirt.obj@Data$data )[ii] )
        print( mirt::itemplot(mirt.obj,item=ii , main=main , ...) )
        graphics::par(ask=ask)
                }
        }