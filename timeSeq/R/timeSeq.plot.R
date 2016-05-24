timeSeq.plot = function(timeSeq.obj, N) {  

	if (N < 1) {
		cat("N must be a positive integer.\n")
		return(cat("ERROR!"))
	}	
	if (N > 10) {
		cat("N should be no greater than 10.\n")
		return(cat("WARNING!"))
	}	
	graphics.off()
	obj = timeSeq.obj
	ans = obj$sorted
    group.length = obj$group.length
    group1.length = obj$group1.length
    group2.length = obj$group2.length   
    reads = obj$reads   
    exon.length = obj$exon.length
  
    p = round(sqrt(N * 2))
    if (p %% 2 == 1) {p = p + 1}
    q = ceiling((2 * N) / p)
    p = p / 2
    l = c(1/p, 1/q)
    now = c(0, 0)
  
    for (i in c(1 : N)) {
        gene = obj$data[obj$gene.names == ans$NPDE_list$genenames[i], ]
        n = ans$NPDE_list$count[i]
        dataset = data.frame(matrix(0, nrow = n *  group.length, ncol = 6))
        dataset[,1] = as.vector((gene[, c(1 : group.length)]))
        dataset[,2] = factor(c(rep( c(1: group1.length), each = n ), rep( c(1: group2.length), each = n )))
        dataset[,3] = factor(rep( c(1,2), c(group1.length * n, group2.length * n) ))
        dataset[,4] = factor(rep( c(1 : n), group.length))
        dataset[,5] = rep(exon.length[obj$gene.names == ans$NPDE_list$genenames[i]], group.length)
        colnames(dataset) = c("Count", "Time", "Group", "ID", "Length")
        ##print(i)
        ##print(dataset)
        ## xyplot(Count/Length ~ Time|Group, data = dataset, groups = ID, type = "l")      
		## for pass the CRAN check
		ID = NULL
        Plot = xyplot(Count/Length ~ Time|Group, 
                      data = dataset, 
                      #main = paste("Top ",ranks[1]," significant NPDE genes",sep = ""),
                      groups = ID, 
                      type = "l",
                      main = paste(ans$NPDE_list$genenames[i])
                      #layout = c(p , q),
                      #scales = list(y = list(relation = 'free'))
                      )
        print(Plot, position=c(now, l + now), more = TRUE)
        if (i %% p == 0) {
            now[2] = now[2] + l[2]
            now[1] = 0
        } else now[1] = now[1] + l[1]
    }
    ##print(data)     
}