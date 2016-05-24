# based on equal "work" per task
defineTasks <- function(n, nSlaves, tasksPerSlave=1, plot=FALSE) {
	
	# calculate how much "work" each value of x will need
	# this is equivilant to nCm(n-m+1,3) - nCm(n-m,3) where n = number of genes and m is in the set [1:n-2]
	m <- 1:(n-2)
	work <- ((n-m)*(n-m-1))/2
	
	slice <- function(v, n){
		subtot <- floor(sum(v)/n)
		cumtot <- cumsum(v)
		p <- rep(0,n-1)
		for(i in 1:(n-1)) {
			p[i] <- max(which(cumtot < (subtot*i)))
		}
		p <- append(p, 0, after=0)
		p <- append(p, length(work))
		
		return(p)
	}
	
	
	slice_boundaries <- slice(work, nSlaves*tasksPerSlave)
	
	tasks <- vector('list')
	colours <- vector()
	for(i in 1:(length(slice_boundaries)-1)) {
		start <- slice_boundaries[i]+1
		end <- slice_boundaries[i+1]
		
		tasks[[length(tasks)+1]] <- list(x=(start:end))
		colours <- append(colours,rep(i,end-start+1))
	}
	
	if (plot) {
		par(mfrow=c(2:1))
		# Colours indicate the sets of gene trios assigned to the same task. By default, there is 1 task per slave CPU
		plot(work, col=colours, pch=20, ylab="Work", xlab="Gene Trio Set (m)", main="Work per Set of Gene Trios")
		
		plot(cumsum(work)/max(cumsum(work))*100, col=colours, pch=20, ylab="Cumlative Work (%)", xlab="Gene Trio Set (m)", main="Cumlative Work for Sets of Gene Trios")
		abline(h=c(cumsum(work)[c(slice_boundaries[2:(length(slice_boundaries)-1)])]/max(cumsum(work))*100))
	}
	
	return(tasks)
}
