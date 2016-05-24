recluster.dist <-function(mat, phylo=NULL, dist="simpson") {
                                                   
	mat<-as.matrix(mat)
	if(is.numeric(mat) == "FALSE") 	stop("ERROR: not numeric matrix")
	if(prod(rowSums(mat)) == 0) stop("ERROR: empty site(s) in matrix")
	#if(prod(colSums(mat)) == 0) stop("ERROR: species with 0 occurrence in matrix")

	options<-c("simpson","sorensen","nestedness","beta3","richness","jaccard","jturnover","jnestedness","phylosor","phylosort","phylosorpd","unifrac","unifract","unifracpd")
	if (is.na(pmatch(dist[1],options))) {
			res <- designdist(mat, method = dist,terms = c("binary"))
			attr(res,"dist") <- dist
			return(res)
	}

	dist <- match.arg(dist,options)

	res <- switch(dist,
		simpson = {
			res <- designdist(mat, method = "(pmin(A,B)-J)/(pmin(A,B))",terms = c("binary"))
			attr(res,"dist") <- "simpson"
			res
		},
		beta3 = {
			res <- designdist(mat, method = "2*((pmin(A,B)-J)/(A+B-J))",terms = c("binary"))
			attr(res,"dist") <- "beta3"
			res
		},
		jaccard = {
			res <- designdist(mat, method = "(A+B-2*J)/(A+B-J)",terms = c("binary"))
			attr(res,"dist") <- "jaccard"
			res
		},
		jturnover = {
			res <- designdist(mat, method = "2*(pmin(A,B)-J)/(J+2*(pmin(A,B)-J))",terms = c("binary"))
			attr(res,"dist") <- "jturnover"
			res
		},
		jnestedness = {
			res <- designdist(mat, method = "(((pmax(A,B)-J)-(pmin(A,B)-J))/(A+B-J))*(J/(J+2*(pmin(A,B)-J)))",terms = c("binary"))
			attr(res,"dist") <- "jnestedness"
			res
		},
		sorensen = {
			res <- designdist(mat, method = "(A+B-2*J)/(A+B)",terms = c("binary"))
			attr(res,"dist") <- "sorensen"
			res
		},
		richness = {
			res <- designdist(mat, method = "(abs((A-J)-(B-J)))/(A+B-J)",terms = c("binary"))
			attr(res,"dist") <- "richness"
			res
		},
		nestedness = {
			res <- designdist(mat, method = "(abs(A-B)/(A+B))*(J/pmin(A,B))",terms = c("binary"))
			attr(res,"dist") <- "nestedness"
			res
		},
		phylosor = {
			res <- phylosor (mat, phylo)
			attr(res,"dist") <- "phylosor"
			res
		},
		phylosort = {
			res <- phylosort (mat, phylo)
			attr(res,"dist") <- "phylosort"
			res
		},
		phylosorpd = {
			res <- phylosorpd (mat, phylo)
			attr(res,"dist") <- "phylosorpd"
			res
		},
		unifrac = {
			res <- unifrac (mat, phylo)
			attr(res,"dist") <- "unifrac"
			res
		},
		unifract = {
			res <- unifract (mat, phylo)
			attr(res,"dist") <- "unifract"
			res
		},
		unifracpd = {
			res <- unifracpd (mat, phylo)
			attr(res,"dist") <- "unifracpd"
			res
		}
	)
}


