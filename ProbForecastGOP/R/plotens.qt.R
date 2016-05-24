"plotens.qt" <-
function(x,y,grid,number,lims,qt.vector){
   x.lim <- c(min(x,na.rm=TRUE),max(x,na.rm=TRUE))
   y.lim <- c(min(y,na.rm=TRUE),max(y,na.rm=TRUE))
   for(i in 1:number){
     title.1 <- qt.vector[i]*100
     title <- paste(title.1,"-th Percentile")
     ens.plot(grid[,,i],lims,x.lim,y.lim,title)
   }
}
