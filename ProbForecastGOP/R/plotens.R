"plotens" <-
function(x,y,grid,number,lims,title,n.pages.4){
   x.lim <- c(min(x,na.rm=TRUE),max(x,na.rm=TRUE))
   y.lim <- c(min(y,na.rm=TRUE),max(y,na.rm=TRUE))
   count <- NULL
   for(j in 1:number){
     count <- n.pages.4 + j
     title.n <- paste(title,count)
     ens.plot(grid[,,j],lims,x.lim,y.lim,title.n)
   }
}
