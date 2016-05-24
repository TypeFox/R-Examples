plot.TCE <-
function(x,align = FALSE,...){
  gcinfo(FALSE)
	
	if(length(x$lambda) == 1)	par(mfrow = c(1, 2))
	if(length(x$lambda) == 2)	par(mfrow = c(1, 3))
	if(length(x$lambda) >= 3)	par(mfrow = c(1, 4))
	
	if(length(x$lambda) <= 3)	z.final = 1:length(x$lambda)
	
	if(length(x$lambda) >=4){
		z.max = max(x$sparsity)
		z.min = min(x$sparsity)
		z = z.max - z.min
		z.unique = unique(c(which(x$sparsity>=(z.min + 0.03*z))[1],which(x$sparsity>=(z.min + 0.07*z))[1],which(x$sparsity>=(z.min + 0.15*z))[1]))

		
		if(length(z.unique) == 1){
			if(z.unique<(length(x$lambda)-1))	z.final = c(z.unique,z.unique+1,z.unique+2)
			if(z.unique==(length(x$lambda)-1)) z.final = c(z.unique-1,z.unique,z.unique+1)
			if(z.unique==length(x$lambda)) 	z.final = c(z.unique-2,z.unique-1,z.unique)
		}
		
		if(length(z.unique) == 2){
			if(diff(z.unique)==1){
				if(z.unique[2]<length(x$lambda)) z.final = c(z.unique,z.unique[2]+1) 
				if(z.unique[2]==length(x$lambda)) z.final = c(z.unique[1]-1,z.unique)
			}
			if(diff(z.unique)>1) z.final = c(z.unique[1],z.unique[1]+1,z.unique[2])
		}
		
		if(length(z.unique) == 3) z.final = z.unique
		
		rm(z.max,z.min,z,z.unique)
		gc()
		
	}
	plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Sparsity vs. Regularization")
	
	lines(x$lambda[z.final],x$sparsity[z.final],type = "p")
	
	if(align){
		layout.grid = layout.fruchterman.reingold(graph.adjacency(as.matrix(x$path[[z.final[length(z.final)]]]), mode="undirected", diag=FALSE))
		for(i in z.final){
			g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
			plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
			rm(g)
			gc()
		}
	rm(layout.grid)
	}
	if(!align){
		for(i in z.final){
			g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
			layout.grid = layout.fruchterman.reingold(g)
			plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
			rm(g,layout.grid)
			gc()
		}
	}
	if(align) cat("Three plotted graphs are aligned according to the third graph\n")   
}
