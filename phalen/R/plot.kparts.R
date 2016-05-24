plot.kparts <-
function(x,...) {
    plotcol = c("#0094D3","#4CB244","#FCCE66","#F38A1E","#DC3410","#9486AC","#77D5FF")
    repx = ceiling(length(x$partitions$parts)/length(plotcol))
    kpcol = as.character(
      merge(data.frame('parts' = x$partitions$parts,
                       'colors' = rep(plotcol,repx)[1:length(x$partitions$parts)]),
            x$data,
            by.x = 'parts',
            by.y = 'parts')$color)
    
    plot(x$data[,2:3], col=kpcol,pch=16)
}
