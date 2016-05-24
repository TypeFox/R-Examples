plot.monte<-function(x,...){
#library(lattice)
z<-x   #monte.object
super.sym <- lattice::trellis.par.get("superpose.symbol")
grp.labs<-as.factor(paste("Group",z$data[,1],sep=" "))
lattice::splom(~z$data[,2:(z$nvar+1)], groups = grp.labs,  
      # panel = panel.superpose,
       key = list(title = "Monte: Scatter Plot Matrices",
                  columns = 3, 
                  points = list(pch = super.sym$pch[1:z$nclus],
                  col = super.sym$col[1:3]),text=list(as.vector(unlist(unique(grp.labs))))))

}
