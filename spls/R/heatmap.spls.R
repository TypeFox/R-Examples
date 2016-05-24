
# heatmap-type MSPE plot (called by cv.spls)

"heatmap.spls" <-
function( mat, coln=16, as='n', ... )
{
    # heatmap-type MSPE plot
    
    par( fig=c(0,7/10,0,1) )
    image( mat, col=cm.colors(coln), axes=FALSE, ... )
    
    # axis label = row / col names
    
    if ( as=='n' )
    {
        axis( 1, seq(0,1,length=nrow(mat)), rownames(mat) )
        axis( 2, seq(0,1,length=ncol(mat)), colnames(mat), las=1 )
    }
    
    # axis label = row / col indices
    
    if ( as=='i' )
    {
        axis( 1, seq(0,1,length=nrow(mat)), c(1:nrow(mat)), tck=FALSE )
        axis( 2, seq(0,1,length=ncol(mat)), c(1:ncol(mat)), tck=FALSE )
    }   
    
    # if as!='n' & as!='i', then no axis
    
    # color bar
    
    par( fig=c(7/10,1,0,1), new=TRUE )
    colvec = matrix( seq( min(mat), max(mat), length=coln ) )
    image( t(colvec), col=cm.colors(coln), axes=FALSE )
    axis( 2, seq(0.0625, 1, length=8),
        format( colvec[seq(2,16,2)], digits=3, nsmall=2 ), las=1 )
        
    par( fig=c(0,1,0,1) )
}
