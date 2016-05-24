#' Read in the members of a gene group
#'
#' This function reads in the names of the members of a group of genes
#' and stores them to a list. The input file needs to be formatted
#' appropriately, with one group per line:
#' groupname<tab>gene1<tab>gene2...geneK
#' 
#' @param filename the file containing the list of gene groups
#' @param sep indicate the appropriate separator if not tab
#'
#' @export
#'
#' @return A list whose names are gene group names and whose elements
#' are vectors of genes
read.groups = function(filename, sep="\t") {
	groups = readLines(filename)
	group_names = vector()
	group_list = list()
	for (i in 1:length(groups)) {
		splitLine = unlist(strsplit(groups[i],sep))
		group_names[i] = splitLine[1]
		group_list[[i]] = splitLine[2:length(splitLine)]
	}	
	names(group_list)=group_names
	return(group_list)
}

#' Find groups represented in the data
#'
#' This function takes in the list of genes for which
#' gene expression data are available as well as the
#' list of gene groups produced by read.groups and a 
#' minimum size (min_size) and returns those genes that 
#' have at least min_size genes with data available
#'
#' @param genes a vector containing the names of each gene for which expression data is available
#' @param groups a list of gene groups, in the same format as the output of read.groups
#' @param min_size the minimum size of groups to be considered
#' 
#' @export
#'
#' @return A vector of group names that had at least min_size genes represented in the data
#'
#' @examples
#' data(yeast)
#' length(GO.groups)
#' GO.groups.pruned = good.groups(colnames(yeast.hybrid),GO.groups,10)
#' length(GO.groups.pruned)
good.groups = function(genes,groups,min_size=2) {
	#genes is just a list of names!
	good = c()
	for (i in 1:length(groups)) {
		if (sum(groups[[i]]%in%genes)>min_size)
		good = c(good,i)
	}
	return(names(groups)[good])
}

#' Remove genes from a specified group from the data
#'
#' This function will remove data corresponding to genes that 
#' are members of groups.to.remove.
#'
#' @param dat gene expression data, rows are species, columns are gene, including colnames.
#' @param groups a list of gene groups, in the same format as the output of read.groups
#' @param groups.to.remove a vector of group names, the members of which will be removed from dat
#'
#' @export
#'
#' @return A new matrix of data with the genes that belonged to groups.to.remove gone.
#'
#' @examples
#' data(yeast)
#' GO.groups.pruned = good.groups(colnames(yeast.hybrid),GO.groups,10)
#' dim(yeast.hybrid)
#' to_remove = setdiff(names(GO.groups),GO.groups.pruned)
#' yeast.hybrid.pruned = remove.from.genome(yeast.hybrid,GO.groups,to_remove)
#' dim(yeast.hybrid.pruned)
remove.from.genome = function(dat,groups,groups.to.remove) {
	dat.removed = dat
	for (i in 1:length(groups.to.remove)) {
		to_remove = which(colnames(dat.removed)%in%groups[[groups.to.remove[i]]])
		if (length(to_remove)!=0) dat.removed = dat.removed[,-to_remove]
	}
	return(dat.removed)
}

#a function to read the expression levels and put them in the right format
#rows are species, columns are genes, and it needs to be in the same order as the tip labels

#' Read gene expression data from a file, sorted by a phylogenetic tree
#' 
#' This funciton will read in gene expression data from a text file
#' and format it for downstream analysis. In particular, it will ensure
#' that species are sorted appropriately given the phylogenetic tree
#' and that the data is appropriately normalized by one species
#' 
#' @param filename expression data file. Ideally, rows correspond to species
#' and columns correspond to genes. First column should be species names, first row
#' should be gene names. If the file is formatted in a transpose (i.e. genes are
#' rows and species are columns) then see the transpose argument
#' @param phy an ape-format phylogenetic tree containing all the species that have
#' gene expressiond data. The tip labels of phy should correspond to the species
#' names in the gene expression data file
#' @param transpose a logical indicating whether to transpose the data from the input file.
#' Only necessary if the input file has rows corresponding to genes and columns corresponding to species
#' @param normalize indicates which species to normalize by. If normalization is undesired (unlikely)
#' set to 0.
#' @param sep the character that separates entries in the expression file
#'
#' @export
#'
#' @return A matrix containing gene expression data, in which rows correspond to species
#' and columns correspond to genes
read.exp = function(filename,phy,transpose=FALSE,normalize=1,sep="\t") {
	if (transpose) {
		genes = t(read.table(filename,row.names=1,header=T))
	} else {
		genes = read.table(filename,row.names=1,header=T)
	}
	genes.reorder = matrix(nrow=0,ncol=ncol(genes))
	tipLabels = phy$tip.label
	for (i in 1:length(tipLabels)) {
		if (normalize) {
			genes.reorder = rbind(genes.reorder,genes[tipLabels[i],]-genes[tipLabels[normalize],])
		} else {
			genes.reorder = rbind(genes.reorder,genes[tipLabels[i],])
		}
	}
	colnames(genes.reorder)=colnames(genes)
	return(genes.reorder)
}


#' Draw normal-inverse-gamma distributed random variates
#'
#' This function will return normal-inverse-gamma distributed random variates
#'
#' @param n the number of variates to generate
#' @param alpha shape parameter of the inverse gamma distribution
#' @param beta scale parameter of the inverse gamma distribution
#'
#' @export
#'
#' @return A vector of random variates arising form the normal-inverse-gamma distribution
rnorminvgamma = function(n,alpha,beta) {
	#draw gamma variances
	vars = 1/rgamma(n,alpha,beta)
	#draw the rnorms
	return(rnorm(n,0,sqrt(vars)))
}

#' Compute the pdf for normal-inverse-gamma random variates
#'
#' This function returns the normal-inverse-gamma density evaluated at specific points
#'
#' @param x a vector of points at which to evaluat the density
#' @param alpha shape parameter of the inverse gamma distribution
#' @param beta scale parameter of the inverse gamma distribution
#' @param log a logical indicated whether to return the log of the pdf
#'
#' @export
#'
#' @return A vector densities corresponding to the entries of x
dnorminvgamma = function(x,alpha,beta,log=FALSE) {
	res = -1/2*log(2*pi)+alpha*log(beta)-(alpha+1/2)*log(x^2/2+beta)+lgamma(alpha+1/2)-lgamma(alpha)
	if (log) {
		return(res)
	} else {
		return(exp(res))
	}
}


#' Compute the pdf for multi-variate normal-inverse-gamma random variates
#'
#' This function returns the multi-variate normal-inverse-gamma density evaluated at specific points
#'
#' @param x a matrix where each row is a sample and each column is a dimension. 
#' @param mu a vector indicating the mean in each dimension
#' @param alpha shape parameter of the inverse gamma distribution
#' @param beta scale parameter of the inverse gamma distribution
#' @param T the variance-covariance matrix
#' @param log a logical indicated whether to return the log of the pdf
#'
#' @export
#'
#' @return A vector densities corresponding to the rows of x
dmvnorminvgamma = function(x,mu,alpha,beta,T,log=FALSE) {
	n = ncol(T)
	distval = mahalanobis(x,center=mu,cov=T)
	res = -1/2*log((2*pi)^n*det(T))+alpha*log(beta)-(alpha+n/2)*log(distval/2+beta)+lgamma(alpha+n/2)-lgamma(alpha)
	if (log) {
		return(res)
	} else {
		return(exp(res))
	}
}

#' Simulate phylogenetic comparative data as Brownian motions with inverse gamma distributed rates
#'
#' This function simulates the evolution of a group of traits evolving
#' as independent Brownian motions with inverse gamma distributed rates.
#' 
#' @param phy an ape format phylogeny on which to simulate
#' @param alpha the shape parameter of the inverse gamma distribution
#' @param beta the scale parameter of the inverse gamma distribution
#' @param rates a vector of rates, with each entry corresponding to an edge
#' of phy. Rates should be in the same order as edges in phy$edge
#' @param n the number of traits to simulate
#'
#' @export
#'
#' @return A matrix with each row corresponding to a species
#' and each column corresponding to a trait
#'
#' @examples
#' data(yeast)
#' norminvgamma.shift.sim.group(yeast.tree,2,2,rep(1,6),10)
norminvgamma.shift.sim.group = function(phy,alpha,beta,rates, n) {
	#create a matrix of nodes
	#columns are nodes, rows are genes
	nodes = matrix(ncol=(phy$Nnode+length(phy$tip.label)),nrow=n)
	nodes[,phy$edge[1,1]] = 0 #set the root for all genes to be zero
	sigma2 = 1/rgamma(n,alpha,beta) #draw n rates
	#construct a covariance matrix for the brownian motions as an outer product
	bmvcv = 0*sqrt(sigma2)%o%sqrt(sigma2)
	diag(bmvcv) = sigma2
	for (i in 1:length(phy$edge[,1])) {
		curLen = phy$edge.length[i]*rates[i]
		nodes[,phy$edge[i,2]] = nodes[,phy$edge[i,1]] + rmvnorm(1,sigma=curLen*bmvcv)
	}
	return(t(nodes[,1:length(phy$tip.label)]))
}

#' Simulate phylogenetic comparative data as OUs with inverse gamma distributed rates
#'
#' This function simulates the evolution of a group of traits evolving
#' as independent Ornstein-Uhlenbeck processes with inverse gamma distributed rates.
#' 
#' @param phy an ape format phylogeny on which to simulate
#' @param theta the strength of constraint
#' @param alpha the shape parameter of the inverse gamma distribution
#' @param beta the scale parameter of the inverse gamma distribution
#' @param n the number of traits to simulate
#'
#' @export
#'
#' @return A matrix with each row corresponding to a species
#' and each column corresponding to a trait
#'
#' @examples
#' data(yeast)
#' OU.invgamma.sim.group(yeast.tree,2,2,2,10)
OU.invgamma.sim.group = function(phy,theta,alpha,beta,n) {
	nodes = matrix(ncol=(phy$Nnode+length(phy$tip.label)),nrow=n)
	sigma2=1/rgamma(n,alpha,beta)
	nodes[,phy$edge[1,1]]=rnorm(n,mean=0,sd=sqrt(sigma2/(2*theta)))
	for (i in 1:length(phy$edge[,1])) {
		curLen = phy$edge.length[i]
		curVar = sigma2/(2*theta)*(1-exp(-2*theta*curLen))
		curMean = nodes[,phy$edge[i,1]]*exp(-theta*curLen)
		curVCV = diag(curVar)
		nodes[,phy$edge[i,2]] = rmvnorm(1,mean=curMean,sigma=diag(curVar))
	}
	return(t(nodes[,1:length(phy$tip.label)]))
}

norminvgamma_shift_like = function(phy,dat,alpha,beta,rates) {
	vcv = updatevcv(phy,rates)
	dat = as.matrix(dat)
	dmvnorminvgamma(t(dat),0,alpha,beta,vcv,log=TRUE)
}


#' Calculate the likelihood of normalized comparative data as Brownian motions with inverse gamma distributed rates
#'
#' This function calculates the likelihood of the observed trait data 
#' assuming that each trait evolves according to an independent Brownian
#' motion with inverse gamma distributed rates. The data are normalized relative
#' to the trait values in a specified species.
#' 
#' @param phy an ape format phylogeny on which to simulate
#' @param dat a matrix of comparative data, in which rows correspond to species
#' and columns correspond to traits
#' @param alpha the shape parameter of the inverse gamma distribution
#' @param beta the scale parameter of the inverse gamma distribution
#' @param rates a vector of rates for each branch of the phylogeny. The order 
#' of elements in rates shoud correspond to the order of phy$branches
#' @param norm the species by which all the data is normalized
#'
#' @export
#'
#' @return A vector, with the likelihood of each gene the observed data
#'
#' @examples
#' data(yeast)
#' sim.dat = norminvgamma.shift.sim.group(yeast.tree,2,2,rep(1,6),10)
#' norminvgamma.shift.like.norm(yeast.tree,sim.dat,2,2,rep(1,6))
norminvgamma.shift.like.norm = function(phy,dat,alpha,beta,rates,norm=1) {
	#make sure the data is normalized!
	dat = apply(dat,2,function(x){x-x[norm]})
	vcv = updatevcv(phy,rates)
	vcvnorm = normalize.vcv(vcv)[-norm,-norm]
	dat = as.matrix(dat)
	dmvnorminvgamma(t(dat[-norm,]),0,alpha,beta,vcvnorm,log=TRUE)
}


#' Calculate the likelihood of normalized comparative data as OUs with inverse gamma distributed rates
#'
#' This function calculates the likelihood of the observed trait data 
#' assuming that each trait evolves according to an independent Ornstein
#' Ulenbeck processes with inverse gamma distributed rates. The data are normalized relative
#' to the trait values in a specified species.
#' 
#' @param phy an ape format phylogeny on which to simulate
#' @param dat a matrix of comparative data, in which rows correspond to species
#' and columns correspond to traits
#' @param alpha the shape parameter of the inverse gamma distribution
#' @param beta the scale parameter of the inverse gamma distribution
#' @param theta the constraint parameter of the Ornstein-Uhlenbeck process
#' @param norm the species by which all the data is normalized
#'
#' @export
#'
#' @return A vector, with the likelihood of each gene the observed data
#'
#' @examples
#' data(yeast)
#' sim.dat = OU.invgamma.sim.group(yeast.tree,2,2,2,10)
#' OU.invgamma.like.norm(yeast.tree,sim.dat,2,2,2)
OU.invgamma.like.norm = function(phy,dat,alpha,beta,theta,norm=1) {
	#make sure the data is normalized
	dat = apply(dat,2,function(x){x-x[norm]})
	vcv = OU.vcv(phy,theta)
	vcvnorm = normalize.vcv(vcv)[-norm,-norm]
	dat = as.matrix(dat)
	dmvnorminvgamma(t(dat[-norm,]),0,alpha,beta,vcvnorm,log=TRUE)
}

one_shift_invgamma_norm_beta = function(pars,phy,dat,branches,norm=1) {
	#make sure that the data is normalized
	dat = apply(dat,2,function(x){x-x[norm]})
	alpha = pars[1]
	beta = pars[2]
	shift_par = pars[3]
	rates = rep(1,nrow(phy$edge))
	rates[branches] = shift_par
	return(-sum(norminvgamma.shift.like.norm(phy,dat,alpha,beta,rates)))
}

OU_invgamma_norm_beta = function(pars,phy,dat,norm=1) {
	#make sure that the data is normalized
	dat = apply(dat,2,function(x){x-x[norm]})
	alpha = pars[1]
	beta = pars[2]
	theta = pars[3]
	return(-sum(OU.invgamma.like.norm(phy,dat,alpha,beta,theta)))
}

#' The mean of the inverse gamma distribution
#'
#' This function returns the mean of inverse gamma distributed random variables
#'
#' @param alpha the shape parameter of the inverse gamma distribution
#' @param beta the scale parameter of the inverse gamma distribution
#'
#' @export
#'
#' @return the mean value of the parameterized inverse gamma disribution,
#' given by beta/(alpha-1) if alpha > 1
mean_invgamma = function(alpha,beta) {
	if (alpha > 1) {
		beta/(alpha-1)
	} else {
		Inf
	}
}

#' The variance of the inverse gamma distribution
#'
#' This function returns the variance of inverse gamma distributed random variables
#'
#' @param alpha the shape parameter of the inverse gamma distribution
#' @param beta the scale parameter of the inverse gamma distribution
#'
#' @export
#'
#' @return the variance of the parameterized inverse gamma disribution,
#' given by beta^2/((alpha-1)^2*(alpha-2)) if alpha > 2
var_invgamma = function(alpha,beta) {
	if (alpha > 2) {
		beta^2/((alpha-1)^2*(alpha-2))	
	} else {
		Inf
	}
}

#' The variance covariance matrix for an Ornstein-Uhlenbeck process
#'
#' This function returns the variance-covariance matrix corresponding to an
#' Ornstein-Uhlenbeck process run along a phylogeny
#'
#' @param phy an ape format phylogeny
#' @param theta the constraint parameter of the Ornstein-Uhlenbeck process
#'
#' @export
#'
#' @return a matrix with nrow = ncol = length(phy$tip.label), where the i,jth entry corresponds
#' to the covariance between species i and j. NB: this is computed assuming that the rate
#' of evolution is equal to one, and can be rescaled simply by multiplying
OU.vcv = function(phy,theta) {
	#theta is STRENGTH OF SELECTION
	times = updatevcv(phy,1)
	V = matrix(nrow=nrow(times),ncol=ncol(times))
	for (i in 1:nrow(V)) {
		for (j in 1:ncol(V)) {
			V[i,j] = 1/(2*theta)*exp(-2*theta*(times[i,i]-times[i,j]))*(1-exp(-2*theta*times[i,j]))
		}
	}
	return(V)
}

#' Compute the variance-covariance matrix of normalized phylogenetic data
#'
#' This function takes a variance-covariance matrix corresponding to some
#' model of trait evolution along a phylogeny and returns the modified
#' variance-covariance matrix that results from normalizing the data by
#' the trait value in one of the species
#'
#' @param vcv origin variance-covariance matrix
#' @param which.norm the species by which the data are normalized
#'
#' @export
#'
#' @return a matrix with nrow = ncol = ncol(vcv).
normalize.vcv = function(vcv,which.norm=1) {
	#which.norm is the index of the taxon that is used to normalize
	new.vcv = vcv
	for (i in 1:nrow(vcv)) { 
		for (j in 1:nrow(vcv)) {
			if (i == j) {
				new.vcv[i,i] = vcv[i,i]+vcv[which.norm,which.norm]-2*vcv[i,which.norm]
			} else {
				new.vcv[i,j] = vcv[i,j] - vcv[i,which.norm]-vcv[j,which.norm]+vcv[which.norm,which.norm]
			}
		} 
	} 
	return(new.vcv)
}

compute_wAIC_once = function(LL,k,n) {
	#LL is a vector with LLs for each model
	#k is a vector with the number of parameters in each model
	#n is the sample size
	weight = vector(length=length(LL))
	aic = 2*LL+2*k#+(2*k*(k+1))/(n-k-1)
	delta = aic - min(aic)
	weight = exp(-delta/2)
	weight = weight/sum(weight)
	return(weight)	
}

compute_wAIC = function(LL,k,n) {
	#LL is a matrix of likelihoods---each row corresponds to a dataset, each column a model
	#k is a vector indicating the number of parameters in each model
	#n is a vector indicating the sample size of each dataset
	weight = LL
	aic = LL
	delta = LL
	for (i in 1:nrow(LL)) {
		print(k)
		print(n[i])
		aic[i,] = 2*LL[i,] + 2*k#+2*k*(k+1)/(n[i]-k-1)
		delta[i,] = aic[i,] - min(aic[i,])
		weight[i,] = exp(-delta[i,]/2)/sum(exp(-delta[i,]/2))
	}
	#aic <<- aic
	#delta <<- delta
	return(weight)
}


#' Compute the vector of branch scaling parameters
#'
#' This function computes the square root of the phylogenetic
#' distance between each species in the tree and one other specified species.
#' Assuming the true model is 
#' Brownian motion with no rate shifts, the distributions of trait 
#' change from a given species to any other species
#' should be identical when divided by the square root of the distance between
#' the two species. 
#'
#' @param phy ape format phylogeny
#' @param norm the species to which distances are computed
#' @param species a vector of species names for which to compute distances
#'
#' @export
#'
#' @return A vector of square root of the distance between each species and norm
compute.sqrt.dist = function(phy,norm = 1,species = phy$tip.label[-norm]) {
	#species is a vector of names of tips
	#norm is a name of a tip
	mrca_mat = mrca(phy)
	#go through each species and get the right number
	normalize = c()
	for (i in 1:length(species)) {
		cur_mrca = mrca_mat[species[i],norm]
		#trace back to the mrca from each one
		cur_branch = which(phy$edge[,2]==which(phy$tip.label==species[i]))
		distance_1 = phy$edge.length[cur_branch]
		while(phy$edge[cur_branch,1]!=cur_mrca) {
			cur_branch = which(phy$edge[,2]==phy$edge[cur_branch,1])
			distance_1 = distance_1 + phy$edge.length[cur_branch]
		}
		cur_branch = which(phy$edge[,2]==norm)
		distance_2 = phy$edge.length[cur_branch]
		while(phy$edge[cur_branch,1]!=cur_mrca) {
			cur_branch = which(phy$edge[,2]==phy$edge[cur_branch,1])
			distance_2 = distance_2 + phy$edge.length[cur_branch]
		}
		
		normalize = c(normalize,sqrt(distance_1+distance_2))
	}
	names(normalize) = species
	return(normalize)
}

plot_density = function(data_matrix,color_vec,names.arg = 1:ncol(data_matrix),from=0,xlab="rate",plot.legend=T,...) {
  dens = list()
  cur_max_dens = 0
  #cur_max_x = 0
  for (i in 1:ncol(data_matrix)) {
  	if (from > -Inf){ 
    	dens[[i]] = density(data_matrix[,i],from=from)
  	} else {
  		dens[[i]] = density(data_matrix[,i])
  	}
    if (max(dens[[i]]$y) > cur_max_dens) {
      cur_max_dens = max(dens[[i]]$y)
    }
    #if (max(dens[[i]]$x) > cur_max_x) {
    #  cur_max_x = max(dens[[i]]$x)
    #}
  }
  plot(dens[[1]],xlab=xlab,ylim=c(0,cur_max_dens),ylab="Density",col=color_vec[1],...)
  for (i in 2:length(dens)) {
    lines(dens[[i]],col=color_vec[i],...)
  }
  if (plot.legend) {
  	legend("topright",legend=names.arg,lty=1,col=color_vec)
  }
}

#' Plot densities of expression differences, possibly normalized
#'
#' This function plots the densities of log fold change in each species
#' relative to a single species. Expression differences may be normalized,
#' to assess the fit of the phylogentic model
#'
#' @param dat a matrix of comparative data, in which rows correspond to species
#' and columns correspond to traits 
#' @param groups a list of gene groups, in the same format as the output of read.groups
#' @param group_name the name of the gene group for which to plot density. Leave as default to use
#' all genes
#' @param normalize a vector of normalizations, corresponding to the order of speices in
#' dat
#' @param remove.row a species of the data matrix to remove. Use this option to ensure that the
#' normalizing species is not plotted
#' @param color a vector of colors in which to plot densities
#' @param main the title of the plot. Default is group name.
#' @param names.arg the name to assign to each density. Default is species name.
#' @param plot.legend a logical indicating whether to plot a legend
#' @param lwd the line width for each density
#'
#' @export
#'
#' @return Nothing
#'
#' @examples
#' data(yeast)
#' sqrt.dist = compute.sqrt.dist(yeast.tree)
#' par(mfrow=c(1,2))
#' test_group = "GO:0007346|regulation of mitotic cell cycle"
#' plot_logfoldchange(yeast.homozygote,GO.groups,test_group)
#' plot_logfoldchange(yeast.homozygote,GO.groups,test_group,normalize=sqrt.dist)
plot_logfoldchange = function(dat,groups,group_name="",normalize=rep(1,nrow(dat)),remove.row=1,color=1:nrow(dat),main=group_name,names.arg=rownames(dat)[-remove.row],plot.legend=T,lwd=1) {
	names.arg=names.arg
	if (remove.row) {
		dat = dat[-remove.row,];
	}
	if (group_name != "") {
		cur_dat = t(dat[,which(colnames(dat)%in%groups[[group_name]])]/normalize)
	} else {
		main = "All genes"
		cur_dat = t(dat/normalize)
	}
	plot_density(cur_dat,names.arg=names.arg,color_vec=color,from=-Inf,main=main,xlab=expression(log[2]~"fold change"),plot.legend=plot.legend,lwd=lwd)
}

rmse = function(truth,estimate) {
	est_mean = mean(estimate)
	est_rmse = sqrt(mean((estimate-truth)^2))
	return(list(mean=est_mean,rmse=est_rmse))
}

descendant_branches = function(node,phy) {
	#returns the rows that are descendants of node
	#node is the node number as listed in phy$edge	
	if (node <= length(phy$tip.label)) {
		return(NULL)
	} else {
		cur_desc = which(phy$edge[,1]==node)
		return(c(cur_desc,descendant_branches(phy$edge[cur_desc[1],2],phy),descendant_branches(phy$edge[cur_desc[2],2],phy)))
	}
}

#' Test all possible single-rate shift Brownian motion models and an Ornstein-Uhlenbeck model
#' 
#' This function will find the maximum likelihood estimate of the parameters
#' of every single rate shift model that is compatible with the phylogeny phy, 
#' as well as the likelihood and wAIC for each model. The procedure is described 
#' in Schraiber et al (2013).
#'
#' @param phy an ape format phylogeny
#' @param dat a matrix of gene expression data. Rows of dat correspond to species
#' and columns of dat correspond to genes
#' @param norm the species by which data should be normalized
#' 
#' @export
#'
#' @return A list of several elements: res is the full output of the optim runs
#' used maximize the likelihood, branches are lists of the branches that have a
#' rate shift for each model, LL is the log likelihood for each model, wAIC is the
#' Akaike information criterion weight for each model, alpha are maximum likelihood estimates
#' of the shape parameter of the inverse gamma distribution for each model, beta are maximum
#' likelihood estimates of the scale paramter for each model and shift are maximum likelihood
#' estimates of the rate shift parameter for each model (except for Ornstein-Uhlenbeck, in which
#' shift is an estimate of the constraint parameter of the OU process).
test.subtrees = function(phy,dat,norm=1) {
	#ROWS OF DAT ARE SPECIEs
	#COLUMNS OF DAT ARE GENES
	#ALSO TESTS THE OU MODEL
	#normalize the data
	dat = apply(dat,2,function(x){x-x[norm]})
	dat = dat-rowMeans(dat)
	internal_nodes = unique(phy$edge[,1]) 
	res = list()
	branches = list()
	#get the descendants of both sides of the root
	root = phy$edge[1,1]
	desc_branch_of_root = which(phy$edge[,1]==root)
	num_branches = length(phy$edge[,1])
	for (node in internal_nodes) {
		#get descendants AND the branch leading to that node
		cur_branches = c(descendant_branches(node,phy), which(phy$edge[,2]==node))
		#if what's left is just one leaf, ignore it and count it in leaves
		if (length(cur_branches)==(num_branches-1)) {
			next
		}
		#if it's the other side of the root, ignore it
		if (desc_branch_of_root[2] %in% cur_branches && length(cur_branches)!=num_branches) {
			next
		}
		branches[[length(branches)+1]] = cur_branches
		res[[length(res)+1]] = optim(c(2,2,2),one_shift_invgamma_norm_beta,gr=NULL,phy,dat,cur_branches,method="L-BFGS-B",lower=.0001,upper=1000)
	}
	all_nodes = unique(phy$edge[,2])
	leaves = setdiff(all_nodes,internal_nodes)
	for (leaf in leaves) {
		cur_branch = which(phy$edge[,2]==leaf)
		branches[[length(branches)+1]] = cur_branch
		res[[length(res)+1]] = optim(c(2,2,2),one_shift_invgamma_norm_beta,gr=NULL,phy,dat,cur_branch,method="L-BFGS-B",lower=.0001,upper=1000)
	}
	#test the OU model
	res[[length(res)+1]] = optim(c(1,1,1),OU_invgamma_norm_beta,gr=NULL,phy,dat,method="L-BFGS-B",lower=.0001,upper=1000)
	LL = unlist(lapply(res,function(x){x$value}))
	k = c(2,rep(3,length(LL)-1)) #NB: OU still has only 3 free parameters, so this is fine!
	n = ncol(dat)
	#print(c(k,n))
	wAIC = compute_wAIC_once(LL,k,n)
	alpha = unlist(lapply(res,function(x){x$par[1]}))
	beta = unlist(lapply(res,function(x){x$par[2]}))
	shift = unlist(lapply(res,function(x){x$par[3]})) #NB: shift is actually theta for the OU
	return(list(res=res,branches=branches,LL=LL,wAIC = wAIC,alpha = alpha, beta = beta, shift = shift))
}

#' Test all possible single-rate shift Brownian motion models and an Ornstein-Uhlenbeck model
#' for an arbitrary number of predefined gene groups.
#' 
#' This function will find the maximum likelihood estimate of the parameters
#' of every single rate shift model that is compatible with the phylogeny phy, 
#' as well as the likelihood and wAIC for each model and for each gene group. The procedure is described 
#' in Schraiber et al (2013).
#'
#' @param phy an ape format phylogeny
#' @param dat a matrix of gene expression data. Rows of dat correspond to species
#' and columns of dat correspond to genes
#' @param groups a list of gene groups, formatted like the output of read.groups
#' @param norm the species by which data should be normalized
#' @param print_names a logical indicating whether to print the name of the group currently
#' being analyzed. Useful to keep track of the progress of the function
#' 
#' @export
#'
#' @return A list of several elements. wAIC is a matrix of Akaike weights for each
#' model for each group (rows are groups, columns are models), alpha is the maximum likelihood shape
#' parameter of the inverse gamma distribution for each model and group, beta is the maximum likelihood
#' scale parameter of the inverse gamma distribution for each model and group, and shift is the maximum
#' likelihood rate shift parameter for each model and each group, except for the final model which is
#' Ornstein-Uhlenbeck, in which case it corresponds to the constraint parameter. Branches indicates
#' which branches of the tree experience a rate shift.
#'
#' @examples
#' \dontrun{
#' data(yeast)
#' GO.groups.pruned = good.groups(colnames(yeast.homozygote),GO.groups,30)
#' to_test = GO.groups[GO.groups.pruned[1:100]]
#' yeast.test = test.groups(yeast.tree,yeast.homozygote,to_test,print_names=T)
#' }
test.groups = function(phy,dat,groups,norm=1,print_names=F) {
	group_res = list()
	for (i in 1:length(groups)) {
		if (print_names) print(names(groups)[i])
		cur_dat = dat[,colnames(dat)%in%groups[[i]]]
		cur_dat = cur_dat-rowMeans(cur_dat)
		group_res[[i]] = test.subtrees(phy,dat[,colnames(dat)%in%groups[[i]]],norm=1)
	}
	wAIC_matrix = matrix(ncol=length(group_res[[1]]$shift),nrow=length(groups))
	rownames(wAIC_matrix) = names(groups)
	alpha_matrix = wAIC_matrix
	beta_matrix = wAIC_matrix
	shift_matrix = wAIC_matrix
	for (i in 1:length(groups)) {
		wAIC_matrix[i,] = group_res[[i]]$wAIC
		alpha_matrix[i,] = group_res[[i]]$alpha
		beta_matrix[i,] = group_res[[i]]$beta
		shift_matrix[i,] = group_res[[i]]$shift
	}
	return(list(wAIC = wAIC_matrix, alpha = alpha_matrix, beta = beta_matrix, shift = shift_matrix, branches = group_res[[1]]$branches))
}


number_of_models = function(phy) {
	num_models = 0
	branches = list()
	internal_nodes = unique(phy$edge[,1])
	root = phy$edge[1,1]
	desc_branch_of_root = which(phy$edge[,1]==root)
	num_branches = length(phy$edge[,1])
	for (node in internal_nodes) {
		#get descendants AND the branch leading to that node
		cur_branches = c(descendant_branches(node,phy), which(phy$edge[,2]==node))
		#if this is just the other side of the root, ignore it
		if (desc_branch_of_root[2] %in% cur_branches && length(cur_branches)!=num_branches) {
			next
		}
		branches[[length(branches)+1]] = cur_branches
		num_models = num_models + 1
	}
	all_nodes = unique(phy$edge[,2])
	leaves = setdiff(all_nodes,internal_nodes)
	#make sure that there isn't a duplicate
	for (i in 1:length(branches)) {
		if (length(branches[[i]]) == (num_branches-1)) {
			#find the child nodes of these branches
			children = c()
			for (branch in branches[[i]]) {
				children = c(children,phy$edge[branch,2])
			}
			missing_leaf = setdiff(all_nodes,children)
			leaves = leaves[-which(leaves==missing_leaf)]
		} 
	}
	num_models = num_models + length(leaves)
	return(num_models)
}

#' Plot a barplot of AIC weights for each model and each gene group
#' 
#' This function will plot a barplot where each bar corresponds to a gene group
#' and the proportion of each bar that is filled with a certain color corresponds to the
#' AIC weight for a given model. Bars are sorted according to which model has the
#' highest weight.
#'
#' @param dat a wAIC matrix, rows are be gene groups
#' and columns are models
#' @param names a vector of names for each model
#' @param col a vector of colors, one color corresponding to each model
#' @param title name for the barplot
#' @param cex character expansion factor
#' @param border type of border around each bar in the barplots
#' @param space amount of space between bars in the barplot
#' @param names.arg names of bars in the barplot
#' 
#' @export
#'
#' @return invisibly returns the ordering of gene groups in the barplot.
#'
#' @examples
#' \dontrun{
#' data(yeast)
#' GO.groups.pruned = good.groups(colnames(yeast.homozygote),GO.groups,30)
#' test_groups = GO.groups[GO.groups.pruned[1:100]]
#' yeast.test = test.groups(yeast.tree,yeast.homozygote,test_groups,print_names=T)
#' plot_wAICbarplot(mytest$wAIC,1:7)	
#' }
plot_wAICbarplot = function(dat,names,col=2:(length(names)+1),title="",cex=1.4,border=par("fg"),space = NULL,names.arg=1:nrow(dat)) {
	#plots the results of the AIC scores. Sorts automatically and stuff
	#dat should be a LIST of wAIC matrices
	## each wAIC matrix should be rows are gene groups, columns are models
	#sort the data
	num_model = length(dat[1,])
	dat.sorted = list()
	orderings = vector(length=length(dat[,1]))

	max_support = apply(dat,1,which.max)
	temp_dat = matrix(nrow=nrow(dat),ncol=ncol(dat))
	cur_start = 1
	for (j in 1:num_model) {
		#figure out which ones have the most support for model j
		cur_groups = which(max_support == j)
		if (length(cur_groups)!=0) {
			cur_end = cur_start+length(cur_groups)-1
			#sort them
			cur_order = order(dat[cur_groups,j],decreasing=T)
			orderings[cur_start:cur_end] = cur_groups[cur_order]
			temp_dat[cur_start:cur_end,] = dat[cur_groups[cur_order],]
			cur_start = cur_end+1
		}
	}
	dat = temp_dat
		
	#get the layout ready
	layout_vec = c(1,1,2)
	layout(matrix(nrow=2,ncol=3,byrow=T,layout_vec),widths=c(1,1,.5),heights=1)
	#plot some things
	barplot(t(dat),names.arg=names.arg,col=col,ylab="Proportion of AIC weight",xlab="",main=title,mar=c(1,1,1,1),border=border,space=space,las=2)
	plot(1,type="n",axes=F,xlab="",ylab="")
	legend("center",legend=names,fill=col,cex=cex,bty="n")
	invisible(orderings)
}

#modified from Jonathan Eastman's package, auteur
updatevcv = function(ape.tre, new.rates) {
	vv=vcv(updatetree(ape.tre, new.rates))
	return(vv)
}

#borrowed from Jonathan Eastman's package, auteur
updatetree = function(ape.tre, new.rates){
	ape.tre$edge.length=ape.tre$edge.length*new.rates
	return(ape.tre)
}
