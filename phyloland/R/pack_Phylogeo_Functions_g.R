#Function that estimates model parameters, genealogies and internal locations through Bayesian Markov chain Monte Carlo (MCMC) algorithms.
#fileTREES: contain tree in nexus format, tip names must be the same as in fileDATA
#fileDATA: for each tip name, give the corresponding space coordinates
PLD_interface <- function(fileTREES, fileDATA, num_step=100000, freq=100, burnin=0, ess_lim=100, sigma=NA, lambda=NA, tau=NA, num_step_sigma=1, num_step_lambda=1, num_step_tau=1, id_filena=NA, pattern_trees_likelihood="treeLikelihood", names_locations=NA){
  if(missing(fileTREES)){
    stop("'fileTREES' is missing")	
  }
  if(missing(fileDATA)){
    stop("'fileDATA' is missing")	
  }
  verif_arg(fileTREES, type = "character", len = 1, name_function = "PLD_interface")
  verif_arg(fileDATA, type = "character", len = 1, name_function = "PLD_interface")
  if (sum(!(is.na(sigma))) != 0){verif_arg(sigma, type = "numeric", len = 2, name_function = "PLD_interface", positive = 2)}
  if (sum(!(is.na(lambda))) != 0){verif_arg(lambda, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)}
  if (sum(!(is.na(tau))) != 0){verif_arg(tau, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)}
  verif_arg(burnin, type = "numeric", len = 1, name_function = "PLD_interface", positive = 1)
  verif_arg(num_step, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  verif_arg(freq, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  verif_arg(num_step_sigma, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  verif_arg(num_step_lambda, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  verif_arg(num_step_tau, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  verif_arg(ess_lim, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  if (sum(!(is.na(id_filena))) != 0){verif_arg(id_filena, type = "character", len = 1, name_function = "PLD_interface")}
  verif_arg(pattern_trees_likelihood, type = "character", len = 1, name_function = "PLD_interface")
  
  for(i in 1:length(readLines(fileDATA))){
    if (length(strsplit(readLines(fileDATA)[i],split="\t")[[1]]) != 3){
      stop("Error in 'fileDATA', is it TAB delimited?")	
    }
  } 
  if (sum(is.na(as.matrix(read.table(fileDATA, header = FALSE, sep="\t"))))){
    stop("Error in 'fileDATA', is it TAB delimited?")	
  }
  
  # read trees phylo
  trees_phylo = read.nexus(fileTREES)
  if (class(trees_phylo) == "phylo"){
    sample_geneal = c(trees_phylo)
    #gtreelikelihood = 1
  }else{
    sample_geneal = trees_phylo[(burnin+1):length(trees_phylo)]
    #trees_par = readLines(fileTREES)[which(regexpr("STATE",readLines(fileTREES)) != -1)]
    #pattern = paste(pattern_trees_likelihood,"=",sep = "")
    # tree likelihood (in BEAST file)
    #gtreelikelihood = lapply(trees_par, treeLikeli, pattern = pattern)
  }
  
  #tips, watch out: because of the way ape read nexus file, the order of tips may be different in each tree
  #unless there is a "Translate" block in nexus file
  tips = sample_geneal[[1]]$tip.label
  space_loc = read_space(fileDATA, tips)	#read space file
  loc_tips = space_loc[[1]] #give for each element of tips the index of location in space
  space = space_loc[[3]]
  
  if (sum(is.na(names_locations)) != 0){
    #names_locations = as.character(c(1:dim(space)[2]))
    names_locations = as.character(tips)
  }else{
    names_locations = as.character(names_locations)
  }
  verif_arg(names_locations, len = length(tips), name_function = "PLD_interface")
  colnames_space = rep("",dim(space)[2])
  for (n in 1:dim(space)[2]){
    colnames_space[n] = names_locations[which(loc_tips==n)[1]]
  }
  colnames(space) = colnames_space
  
  # gtreel and locations_sim for mcmc_phyloland
  gtreel = converttree(sample(sample_geneal, 1)[[1]])
  locations_sim = cbind(rep(0,dim(gtreel$nodes)[1]),rep(0,dim(gtreel$nodes)[1]),rep(0,dim(gtreel$nodes)[1]))
  locations_sim[1:length(loc_tips),1] = loc_tips
  
  # files creation
  if (is.na(id_filena)){
    #id_filena = format(Sys.time(),"%a%d%b%Y_%H%M%S")
    filena = NA
    filena_param = NA
    filena_loc = NA
    filena_tracer = NA
    filena_tree = NA
  }else{
    filena = paste('phyloland_',id_filena,'.csv',sep='')
    filena_param = paste('parameters_',id_filena,'.csv',sep='')
    filena_loc = paste('loc_',id_filena,'.log',sep='')
    filena_tracer = paste("tracer_",id_filena,".log",sep="")
    filena_tree = paste("tree_",id_filena,".tre",sep="")
  }
  
  # no output filename given, should not be allowed...
  if (!is.na(filena_loc)){
    nodes = c(1:length(gtreel$nodes[,1]))
    nodes[1:length(tips)] = tips
    write.table(t(c("num_step",nodes)),row.names=FALSE,col.names=FALSE,file=filena_loc,sep="\t",eol="\n",append=FALSE)
  }
  
  # MCMC
  mcmco = mcmc_phyloland(space=space, gtreel=gtreel, simul_values=list(sigma,lambda,tau,locations_sim), treelikelihood=NA, est_sigma=c(1,1), est_lambda=1, est_tau=1, sample_geneal=sample_geneal, sample_loc=1, plot_phylo=0, plot_MCMC=0, save_MCMC=0, Nstep=num_step, freq=freq, Nstep_sigma=num_step_sigma, Nstep_lambda=num_step_lambda, Nstep_loc=1, Nstep_genealogy=1, Nstep_tau=num_step_tau, pchange=0, ess_lim=ess_lim, model=4, show_loc=0, filena_loc, filena_tracer, filena_tree)
  
  stat_res = stat_phyloland(mcmco[[1]])
  stat_loc = vector("numeric",7)
  
  trees = mcmco[[3]]
  #sampled_loc = read_loc(filena_loc)
  sampled_loc = mcmco[[2]][,2:ncol(mcmco[[2]])]
  ind = mcmco[[6]]
  rownames(space) = c("Dimension1", "Dimension2")
  x <- (list( trees = trees, locations = sampled_loc, tips = tips, space = space, index = ind, mcmc = mcmco, sigma_limit = round(unlist(mcmco[[8]]),3)))
  class(x) <- "phyloland"	
  # phyloland object
  #save(x, file = paste("phyloland_",id_filena,".Rdata",sep = ""))
  
  # files
  if (!is.na(filena_param)){
    cat(paste('num_location',sep=","),file=filena_param,sep=",",append=FALSE)
    for (ndim in 1:dim(space)[1]){
      cat(paste(',prior_sigma',ndim,sep=""),file=filena_param,sep=",",append=TRUE)
    }
    cat(paste(',prior_lambda', 'prior_tau', 'num_step', 'freq', 'num_step_sigma', 'num_step_lambda', 'num_step_tau', 'ess_lim', '\n', sep=","),file=filena_param,sep=",",append=TRUE)
    
    cat(dim(space)[2],file=filena_param,sep=",",append=TRUE)
    
    if (sum(is.na(sigma)) != 0){
      for (ndim in 1:dim(space)[1]){
        cat(paste(',unif',sep=","),file=filena_param,sep=",",append=TRUE)
      }
    }else{
      cat('',sigma,file=filena_param,sep=",",append=TRUE)
    }
    if (is.na(lambda)){
      cat(paste(',unif',sep=","),file=filena_param,sep=",",append=TRUE)
    }else{
      cat('',lambda,file=filena_param,sep=",",append=TRUE)
    }	
    if (is.na(tau)){
      cat(paste(',unif',sep=","),file=filena_param,sep=",",append=TRUE)
    }else{
      cat('',tau,file=filena_param,sep=",",append=TRUE)
    }		
    cat(c('',num_step, freq, num_step_sigma, num_step_lambda, num_step_tau, ess_lim,'\n'),file=filena_param,sep=",",append=TRUE)
  }
  
  if (!is.na(id_filena)){
    for (ndim in 1:dim(space)[1]){
      if (ndim == 1){
        cat(paste(paste('min_sigma',ndim,sep=""),'quant05','quant25','quant50','quant75','quant95','max','mode','mean','last',sep=","),file=filena,sep=",",append=FALSE)
      }else{
        cat(paste(paste(',min_sigma',ndim,sep=""),'quant05','quant25','quant50','quant75','quant95','max','mode','mean','last',sep=","),file=filena,sep=",",append=TRUE)
      }
    }
    cat(paste(',min_lambda','quant05','quant25','quant50','quant75','quant95','max_lambda','mode_lambda','mean_lambda','last_lambda',sep=","),file=filena,sep=",",append=TRUE)
    cat(paste(',min_tau','quant05','quant25','quant50','quant75','quant95','max_tau','mode_tau','mean_tau','last_tau',sep=","),file=filena,sep=",",append=TRUE)
    for (ndim in 1:dim(space)[1]){
      cat(paste(',sigma_upper_limit',ndim,sep=""), file=filena, sep=",", append=TRUE)
    }
    cat('\n',file=filena, sep=",",append=TRUE)
    cat(c(unlist(stat_res),round(unlist(mcmco[[8]]),3),'\n'),file=filena,sep=",",append=TRUE)
  }
  
  return(x)
}


# Function that plots sampled trees with their locations.
PLD_plot_trees <- function(x, sub_sample = 0, one_plot = FALSE){
  if (missing(x)){
    stop("'x' is missing")	
  }
  verif_arg(x, type = "phyloland", name_function = "plot_trees")
  if(class(x$trees) == "phylo"){
    x$trees = c(x$trees)	
  }	
  verif_arg(sub_sample, type = c("numeric","integer"), name_function = "plot_trees")
  verif_arg(one_plot, type = "logical", len = 1, name_function = "plot_trees")
  trees = x$trees
  verif_arg(x$trees, type = "multiPhylo" , name_function = "plot_trees")
  sampled_loc = x$locations
  
  if(length(sub_sample) == 1){
    if (sub_sample == 0){
      sub_sample = 1:length(trees)	
    }		
  }	
  if (sum(is.element(sub_sample,1:length(trees))==FALSE) != 0){
    stop(paste(sub_sample[which(is.element(sub_sample,1:length(trees))==FALSE)[1]]," : subscript out of bounds"))
  }
  
  s = sub_sample
  
  n = length(s)
  if (one_plot == TRUE ){
    nb_row = round(sqrt(n))
    nb_col = ceiling(sqrt(n))
    par(mfrow = c(nb_row,nb_col))
  }
  for (i in 1:n){
    print(paste("tree",s[i]),sep = "")
    plot_trees_tips(trees[[s[i]]],sampled_loc[s[i],],x$space)
    if(one_plot == FALSE){
      if (i < n){
        readline()
      }
    }
  } 
}

# Function that displays the density plots of the migration times for each location.
PLD_loc_mrca <- function(x, tips, sub_sample = 0, plot_distrib = FALSE, col = NA){
  if (missing(x)){
    stop("'x' is missing")	
  }
  verif_arg(x, type = c("phyloland"), name_function = "loc_mrca")
  if(class(x$trees) == "phylo"){
    x$trees = c(x$trees)	
  }
  trees = x$trees
  verif_arg(x$trees, type = "multiPhylo" , name_function = "loc_mrca")
  sampled_loc = x$locations
  corresp = x$tips
  space = x$space
  names_locations = colnames(space)
  verif_arg(plot_distrib, type = "logical", len = 1, name_function = "loc_mrca")
  verif_arg(sub_sample, type = c("numeric","integer"), name_function = "loc_mrca")
  if(length(sub_sample) == 1){
    if (sub_sample == 0){
      sub_sample = 1:length(trees)	
    }		
  }
  if (sum(is.element(sub_sample,1:length(trees))==FALSE) != 0){
    stop(paste(sub_sample[which(is.element(sub_sample,1:length(trees))==FALSE)[1]]," : subscript out of bounds"))
  }
  trees = trees[sub_sample]
  sampled_loc = sampled_loc[sub_sample,]
  
  if (missing(tips)){
    tips = corresp
  }else{
    tips = unique(tips)	
  }
  verif_arg(tips, len = c(1:length(corresp)), name_function = "loc_mrca")
  tips_ind = vector("numeric",length(tips))
  tips = as.character(tips)	
  for (i in 1:length(tips)){
    if (length(which(corresp == tips[i]))==0){
      stop(paste(" tip \"" , tips[i],"\" not found", sep = ""))
    }else{
      tips_ind[i] = which(corresp == tips[i])
    }
  }	
  if (sum(is.na(col))!=0){
    col = NULL	
  }
  loc_MRCA = vector("numeric",length(trees))	
  if (length(tips) == 1){
    mat = c(space[1,sampled_loc[1,tips_ind]], space[2,sampled_loc[1,tips_ind]], 1)
    loc_MRCA = rep(sampled_loc[1,tips_ind],length(trees))
    res = list(frequencies = mat , locationsMRCA = loc_MRCA)
    if (plot_distrib == TRUE){
      tab <- table(loc_MRCA)/sum(table(loc_MRCA))
      tab_complete <- vector("numeric",dim(space)[2])
      tab_complete[as.integer(rownames(tab))] = tab
      names(tab_complete) = names_locations
      #x11() #dev.new()
      barplot(tab_complete, main = "", ylab = "frequencies", xlab = "locations" , col = col)
    }
    return(res)	
  }else{
    for (i in 1:length(trees)){
      loc_MRCA[i] = loc_ancestor(converttree(trees[[i]]), tips_ind, sampled_loc[i,])
    }
    mat <- matrix(ncol = 3, nrow = length(table(loc_MRCA)), dimnames = list(as.character(names(sort(table(loc_MRCA),decreasing=TRUE))),c(rownames(space),"Frequency")))
    sort_loc = as.integer(names(sort(table(loc_MRCA),decreasing=TRUE)))
    mat[,1] = space[1,sort_loc]
    mat[,2] = space[2,sort_loc]
    mat[,3] = sort(table(loc_MRCA),decreasing=TRUE) / sum(table(loc_MRCA))
    rownames(mat) = names_locations[as.integer(rownames(mat))]
    if (plot_distrib == TRUE){				
      tab <- table(loc_MRCA)/sum(table(loc_MRCA))
      tab_complete <- vector("numeric",dim(space)[2])
      tab_complete[as.integer(rownames(tab))] = tab
      names(tab_complete) = names_locations
      #x11() #dev.new()
      barplot(tab_complete, ylab = "frequencies", xlab = "locations" , col = col)
    }
    #if (plot_distrib[2] == TRUE){
    #			x11()
    #			loc_lat = vector("numeric", length(loc_MRCA))
    #			loc_long = vector("numeric", length(loc_MRCA))
    #			for(i in 1:length(loc_MRCA)){
    #				loc_lat[i] = space[1,loc_MRCA[i]]
    #				loc_long[i] = space[2,loc_MRCA[i]]
    #			}
    #			kde = kde2d(loc_lat, loc_long, n = 300)
    #			image(kde)
    #			points(space[1,],space[2,],cex=2, pch = 16)
    #		}
    res = list(frequencies = mat, locationsMRCA = loc_MRCA)
    return(res)
  }
}

# Function that lists the migration events from the genealogies sampled by the MCMC.
PLD_stat_mig <- function(x, sub_sample = 0, first = FALSE){
  if (missing(x)){
    stop("'x' is missing")	
  }
  verif_arg(x, type = "phyloland", name_function = "stat_mig")
  if(class(x$trees) == "phylo"){
    x$trees = c(x$trees)	
  }
  verif_arg(sub_sample, type = c("numeric","integer"), name_function = "stat_mig")
  verif_arg(first, type = "logical", len = 1, name_function = "stat_mig")
  trees = x$trees
  verif_arg(x$trees, type = "multiPhylo" , name_function = "stat_mig")
  sampled_loc = x$locations
  space = x$space
  names_locations = colnames(space)
  
  if(length(sub_sample) == 1){
    if (sub_sample == 0){
      sub_sample = 1:length(trees)	
    }		
  }
  if (sum(is.element(sub_sample,1:length(trees))==FALSE) != 0){
    stop(paste(sub_sample[which(is.element(sub_sample,1:length(trees))==FALSE)[1]]," : subscript out of bounds"))
  }
  trees = trees[sub_sample]
  sampled_loc = sampled_loc[sub_sample,]
  
  migmat = matrix(0,nrow=dim(space)[2],ncol=dim(space)[2])
  if (first == FALSE){
    timemat = array(0,c(dim(space)[2],dim(space)[2],(sum(is.na(converttree(trees[[1]])[,2]))-1)*length(trees)))
  }else{
    timemat = array(0,c(dim(space)[2],dim(space)[2],length(trees)))
  }
  for (m in 1:length(sampled_loc[,1])){
    treel = reorder_treel(converttree(trees[[m]]))[[1]]
    history = treel2mig(treel,sampled_loc[m,],space)
    # record location at the root = ancestral location of MRCA
    timemat[history[1,1],history[1,1],m] = history[1,5+dim(space)[1]] # ! will be confounded with the first migration event
    for (n in 1:dim(history)[1]){
      migmat[history[n,1],history[n,2]] = migmat[history[n,1],history[n,2]] + 1
      if (first == FALSE){
        timemat[history[n,1],history[n,2],which(timemat[history[n,1],history[n,2],]==0)[1]] = history[n,5+dim(space)[1]]
      }else{
        if (sum(timemat[,history[n,2],m])==0){
          timemat[history[n,1],history[n,2],m] = history[n,5+dim(space)[1]]
        }
      }
    }
  }
  rownames(timemat) = names_locations
  colnames(timemat) = names_locations
  rownames(migmat) = names_locations
  colnames(migmat) = names_locations
  statmig = list(migmat = migmat, timemat = timemat)
  return(statmig)
}

#Function that displays the density plots of the migration times for each location.
PLD_plot_stat_mig <- function(timemat, color = NA, xy_legend = NA, group = 0){
  if (missing(timemat)){
    stop("'timemat' is missing")	
  }
  verif_arg(timemat, type = c("array"), name_function = "plot_stat_mig")
  if (sum(is.na(xy_legend)) == 0){
    verif_arg(xy_legend, type = "numeric", len = 2, name_function = "plot_stat_mig")
  }
  if (sum(is.na(color)) == 0){
    verif_arg(color, len = nrow(timemat), name_function = "plot_stat_mig")
  }
  if (group==1) {
    # Group locations according to row names
    timemat2 = array(0,c(length(unique(rownames(timemat))),length(unique(rownames(timemat))),length(timemat[1,1,])))
    rownames(timemat2) = unique(rownames(timemat))
    colnames(timemat2) = unique(colnames(timemat))
    for (n in 1:length(timemat[1,1,])){
      for (loc1 in unique(rownames(timemat))){
        for (loc2 in unique(colnames(timemat))){
          timemat2[loc1,loc2,n] = max( timemat[which(rownames(timemat[,,1])==loc1),which(colnames(timemat[,,1])==loc2),n] )
        }
      }
    }
    timemat = timemat2
  }
  junk.x = NULL
  junk.y = NULL
  for(i in 1:(nrow(timemat))){
    if (length(timemat[,i,][timemat[,i,]!=0])>0){
      junk.x = c(junk.x, density(timemat[,i,][timemat[,i,]!=0])$x)
      junk.y = c(junk.y, density(timemat[,i,][timemat[,i,]!=0])$y)
    }
  }
  xr <- range(junk.x)
  yr <- range(junk.y)
  if (sum(is.na(color)) != 0){
    color = rainbow(nrow(timemat))
  }	
  tmp = density(timemat[,1,][timemat[,1,]!=0])
  ymax = max(tmp$y)
  plot(tmp, col = color[1], xlim = xr, ylim = yr, main = "", ylab = "Density", xlab ="Time (BP)", lwd = 2)
  if (nrow(timemat) > 1){
    for(i in 2:nrow(timemat)){
      if (length(timemat[,i,][timemat[,i,]!=0])>0){
        tmp = density(timemat[,i,][timemat[,i,]!=0])
        if (max(tmp$y)>ymax) ymax = max(tmp$y)
        lines(tmp, col = color[i], lwd = 2)
      }
    }
  }
  if (sum(is.na(xy_legend))!=0){
    xy_legend = c(0.75*max(xr),max(yr))	
  }
  legend(xy_legend[1], xy_legend[2], colnames(timemat),color)
}


################ INTERNAL FUNCTIONS ####################

### convert louis tree structure to the history of migrations
treel2mig <- function(gtreel, location, space){
  check_treel(gtreel, location)
  history = matrix(0, 0, 5+dim(space)[1]) ## |dep|dest|proba|time2father|dist(|dist...)|height
  rootn = which(is.na(gtreel$nodes[,1]))
  ## get the absolute height of each internal node
  inode_height = internal_node_height(gtreel)
  ## need to treat the nodes from the oldest to the most recent
  node_id = sort(inode_height, decreasing=TRUE, index.return=TRUE)$ix
  if (node_id[1] != rootn){stop('treel2mig error')}
  for (n in 1:length(node_id)) {
    ploc = location[[node_id[n]]]
    siblings = which(gtreel$nodes[,1] == node_id[n])
    if (length(siblings) == 0) {next}
    if (location[[siblings[1]]] == ploc && location[[siblings[2]]] != ploc) {
      noddest = siblings[2]
      locdest = location[[siblings[2]]]
    } else if (location[[siblings[2]]] == ploc && location[[siblings[1]]] != ploc) {
      noddest = siblings[1]
      locdest = location[[siblings[1]]]
    } else if (location[[siblings[1]]] == ploc && location[[siblings[2]]] == ploc) {
      noddest = siblings[1]
      locdest = location[[siblings[1]]]
    }
    distmig = as.vector(abs(space[,ploc]-space[,locdest]))
    history = rbind(history, c(ploc, locdest, NA, gtreel$nodes[noddest,4], distmig, inode_height[node_id[n]]))
  }
  return(history)
}

### convert louis tree structure to parenthesis format
treel2par <- function(gtreel){
	tree_phylo = lconverttree(gtreel)
	return(write.tree(tree_phylo))
}

### convert ape phylo tree structure (edges) to node index:parent|left child|right child|branch length to parent|
converttree <- function(tree_phylo){
	#gtreel = matrix(NA, (2*tree_phylo$Nnode + 1), 4)
  gtreel <- list( nodes=matrix(NA,(2*tree_phylo$Nnode+1),4), nodes.label=seq(1,(2*tree_phylo$Nnode+1)) )
	for (n in 1:(2*tree_phylo$Nnode)) {
		gtreel$nodes[tree_phylo$edge[n,2],1] = tree_phylo$edge[n,1] # parent
		if (is.na(gtreel$nodes[tree_phylo$edge[n,1],2])) { # has that edge been visited yet
			gtreel$nodes[tree_phylo$edge[n,1],2] = tree_phylo$edge[n,2] # child left
		} else {
			gtreel$nodes[tree_phylo$edge[n,1],3] = tree_phylo$edge[n,2] # child right
		}
		gtreel$nodes[tree_phylo$edge[n,2],4] = tree_phylo$edge.length[n]
		if (tree_phylo$edge[n,2]<=(tree_phylo$Nnode+1)) { gtreel$nodes.label[tree_phylo$edge[n,2]] = tree_phylo$tip.label[tree_phylo$edge[n,2]] } # tip
		if (!is.null(tree_phylo$node.label)) { gtreel$nodes.label[tree_phylo$edge[n,1]] = tree_phylo$node.label[tree_phylo$edge[n,1]-(tree_phylo$Nnode+1)] }
	}
	return(gtreel)
}

### convert louis tree structure to ape phylo tree structure (edges), assuming tips appear first in gtreel
lconverttree <- function(gtreel,alphab=0){
  tipid = which(is.na(gtreel$nodes[,2]))
  if (alphab==1){	tipid = as.integer(sort(as.character(tipid))) } #tips sorted alphabetically
  n = length(tipid)
  Nnode = n - 1
  nodeid = rep(0, Nnode)
  nodeid[1] = which(is.na(gtreel$nodes[,1])) # root
  edge = matrix(NA, nrow(gtreel$nodes)-1, 2)
  edgelength = vector('numeric', nrow(gtreel$nodes)-1)
  edgeid = 1
  for (i in 1:nrow(gtreel$nodes)){
    if (is.na(gtreel$nodes[i,1])) {i = i + 1 ; next} # root
    x = which(nodeid == gtreel$nodes[i,1])
    if (length(x) == 0){
      x = match(0, nodeid) # next free nodelabel
      nodeid[x] = gtreel$nodes[i,1]
    }
    edge[edgeid,1] = x + n
    z = which(tipid == i)
    if (length(z) == 0){ # is it a tip?
      y = which(nodeid == i)
      if (length(y) == 0){
        y = match(0, nodeid) # next free nodelabel
        nodeid[y] = i
      }
      edge[edgeid,2] = y + n
    }else{
      edge[edgeid,2] = z
    }
    edgelength[edgeid] = gtreel$nodes[i,4]
    edgeid = edgeid + 1
    i = i + 1
  }
  #tree_phylo = list(edge = edge, tip.label = paste("t",tipid,sep=""), node.label = paste("t", nodeid, sep=""), Nnode = Nnode, edge.length = edgelength)
  tree_phylo = list(edge = edge, tip.label = gtreel$nodes.label[tipid], node.label = gtreel$nodes.label[nodeid], Nnode = Nnode, edge.length = edgelength)
  class(tree_phylo) = "phylo"
  return(tree_phylo)
}
### get all possible locations for internal nodes
internal_loc <- function(gtreel, loc_tips){
	locations = vector("list", length(gtreel$nodes[,1]))
	for (n in 1:length(gtreel$nodes[,1])) {
		if (is.na(gtreel$nodes[n,2]) & is.na(gtreel$nodes[n,3])) { ## tip
			locations[[n]] = loc_tips[n]
		} else {
			locations[[n]] = unique(c(ftips(gtreel$nodes[n,2], gtreel, loc_tips), ftips(gtreel$nodes[n,3], gtreel, loc_tips)))
		}
	}
	return(locations)
}

### return the set of location of the tips below a node
ftips <- function(node, gtreel, loc_tips){
	if (is.na(gtreel$nodes[node,2]) & is.na(gtreel$nodes[node,3])) {
		return(loc_tips[node])
	} else {
		return(unique(c(ftips(gtreel$nodes[node,2], gtreel, loc_tips), ftips(gtreel$nodes[node,3], gtreel, loc_tips))))
	}
}

### return the set of nodes (including tips) below a given node
children <- function(node, gtreel){
	if (is.na(gtreel$nodes[node,2]) & is.na(gtreel$nodes[node,3])) {
		return(node)
	} else {
		return(c(node, children(gtreel$nodes[node,2], gtreel), children(gtreel$nodes[node,3], gtreel)))
	}
}

### return the set of tips below a given node
childrentip <- function(node, gtreel){
	if (is.na(gtreel$nodes[node,2]) & is.na(gtreel$nodes[node,3])) {
		return(node)
	} else {
		return(c(childrentip(gtreel$nodes[node,2], gtreel), childrentip(gtreel$nodes[node,3], gtreel)))
	}
}

### get the node absolute heights of internal nodes from branch length to parent
internal_node_height <- function(gtreel){
	inode_height = vector("numeric", length(gtreel$nodes[,1]))
	for (n in 1:length(gtreel$nodes[,1])) {
		if (is.na(gtreel$nodes[n,2]) & is.na(gtreel$nodes[n,3])) { ## tip
			inode_height[n] = 0
		} else {
			inode_height[n] = gtreel$nodes[gtreel$nodes[n,2],4] + height(gtreel$nodes[n,2], gtreel)
		}
	}
	return(inode_height)
}

### return the height of a node
height <- function(node, gtreel){
	if (is.na(gtreel$nodes[node,2]) & is.na(gtreel$nodes[node,3])) {
		return(0)
	} else {
		return(gtreel$nodes[gtreel$nodes[node,2],4] + height(gtreel$nodes[node,2], gtreel))
	}
}

### Compute likelihood of history of migration given a tree
### Simulate internal locations
### gtreel: tree in format "louis"
### optional: a list of possible locations for each internal nodes and locations of tips
### getL: if ==1, return likelihood of a specific set of locations
### model: each node is 2 dispersals potentially on father's loc (1) or exactly one dispersal out of father's loc (2) or one mig per node (4)
### samunif: if ==1, draw the locations of internal nodes uniformly in the list of possible locations, otherwise use the model
sim_history <- function(gtreel, space, sigma, lambda, tau, locations = 0, getL = 0, model = 1, samunif = 0, numtip = 50){
	history_proba_val = 0
	if (model==4 && locations[[1]][1]>0) { ## only one migration per internal node possible
		rootn = which(is.na(gtreel$nodes[,1]))
		locations_sampled = matrix(0,nrow=length(gtreel$nodes[,1]),ncol=(1+(dim(space)[1]))) ## location(destination) | distance(s)
		## get the absolute height of each internal node
		inode_height = internal_node_height(gtreel)
		## need to treat the nodes from the oldest to the most recent
		node_id = sort(inode_height,decreasing=TRUE,index.return=TRUE)$ix
		if (node_id[1]!=rootn){stop('sim_history error root')}
		## keep track of colonised locations
		occupied = vector('numeric',ncol(space))
		if (length(locations[[rootn]])>1){ ## uniform to get the location at the root
			locations_sampled[rootn,1:2] = c( sample(locations[[rootn]],1), 0)
		}else{
			locations_sampled[rootn,1:2] = c( locations[[rootn]], 0)
		}
		occupied[locations_sampled[rootn,1]] = 1
		for (n in 1:length(node_id)) {
			if ( getL==1 || getL==2 ){ ploc = locations[[node_id[n]]]
			}else{ ploc = locations_sampled[node_id[n],1]} ## parent's location
			if ( getL==1 || getL==2 ){## return likelihood of a specific migration event out of the possible ones
				siblings = which(gtreel$nodes[,1]==node_id[n])
				if (length(siblings)==0) {next}
				if (getL==2 && length(which(gtreel$nodes[,1]==siblings[1]))==0 && length(which(gtreel$nodes[,1]==siblings[2]))==0) {next} ## hastings for loc sampling, ignore last migrations
				if (locations[[siblings[1]]]==ploc && locations[[siblings[2]]]!=ploc) {
					noddest = siblings[2]
					locdest = locations[[siblings[2]]]
				} else if (locations[[siblings[2]]]==ploc && locations[[siblings[1]]]!=ploc) {
					noddest = siblings[1]
					locdest = locations[[siblings[1]]]
				} else if (locations[[siblings[1]]]==ploc && locations[[siblings[2]]]==ploc) {
					noddest = siblings[1]
					locdest = locations[[siblings[1]]]
				} else if (locations[[siblings[1]]]!=ploc && locations[[siblings[2]]]!=ploc) { stop('sim_history error 3') }
				if (n==1){
					migration_event = new_migration( space, occupied, sigma, lambda, tau, 0, 0, c(ploc,locdest,0) )
				}else{
					migration_event = new_migration( space, occupied, sigma, lambda, tau, 0, 0, c(ploc,locdest,abs(inode_height[node_id[n]]-inode_height[node_id[n-1]])) )
				}
#if (migration_event[3]>1) print(migration_event)
				history_proba_val = history_proba_val + log(migration_event[3])
				occupied[locdest] = occupied[locdest] + 1
			}else{## sample a migration from the father's location if children nodes have more than 1 possible location, according to model (samunif!=1) or uniformly in the list of possible node location (samunif==1)
				migration_event = c(0,0,1,0,vector('numeric',dim(space)[1])) ## default: proba equals 1
				siblings = which(gtreel$nodes[,1]==node_id[n])
				loc_sib1 = locations[[siblings[1]]]
				loc_sib2 = locations[[siblings[2]]]
				if (length(loc_sib1)==1 && length(loc_sib2)==1){
					locations_sampled[siblings[1],1] = loc_sib1
					locations_sampled[siblings[2],1] = loc_sib2
					if (loc_sib1==ploc) {occupied[loc_sib2] = occupied[loc_sib2] + 1}
					else if (loc_sib2==ploc) {occupied[loc_sib1] = occupied[loc_sib1] + 1}
					else {stop('sim_history error 0')}
					## which one is the parent's
				} else if (length(loc_sib1)>1 && length(loc_sib2)==1){
					if (loc_sib2==ploc) {
						if (samunif==1){ migration_event[2] = sample(loc_sib1,1) }
						else {
							migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib1, 0 )
							}
						locations_sampled[siblings[1],] = c(migration_event[2],migration_event[5:length(migration_event)])
						occupied[migration_event[2]] = occupied[migration_event[2]] +1
					}else{ ## loc_sib1 must contain ploc
						if (sum(loc_sib1==ploc)==0) stop('sim_history error 1')
						locations_sampled[siblings[1],1] = ploc
						occupied[loc_sib2] = occupied[loc_sib2] + 1
					}
					locations_sampled[siblings[2],1] = loc_sib2
				} else if (length(loc_sib1)==1 && length(loc_sib2)>1){
					if (loc_sib1==ploc) {
						if (samunif==1){ migration_event[2] = sample(loc_sib2,1) }
						else {
							migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib2)
							}
						locations_sampled[siblings[2],] = c(migration_event[2],migration_event[5:length(migration_event)])
						occupied[migration_event[2]] = occupied[migration_event[2]] + 1
					}else{ ## loc_sib2 must contain ploc
						if (sum(loc_sib2==ploc)==0) stop('sim_history error 2')
						locations_sampled[siblings[2],1] = ploc
						occupied[loc_sib1] = occupied[loc_sib1] + 1
					}
					locations_sampled[siblings[1],1] = loc_sib1
				} else if (sum(loc_sib1==ploc)>0 && sum(loc_sib2==ploc)>0){ ## draw which daughter stays according to number of leaves below that are father's location
					numer1 = sum(loc_sib1==ploc)
					numer2 = sum(loc_sib2==ploc)
					denom = numer1 + numer2
					sib_stays = sample(1:2,1,prob=c(numer1/denom,numer2/denom))
					if (sib_stays==1){
						locations_sampled[siblings[1],1] = ploc
						if (samunif==1){ migration_event[2] = sample(loc_sib2,1) }
						else {
							migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib2 , 0 )
							}
						locations_sampled[siblings[2],] = c(migration_event[2],migration_event[5:length(migration_event)])
						occupied[migration_event[2]] = occupied[migration_event[2]] + 1
					}else{
						locations_sampled[siblings[2],1] = ploc
						if (samunif==1){ migration_event[2] = sample(loc_sib1,1) }
						else {
							migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib1 , 0  )
							}
						locations_sampled[siblings[1],] = c(migration_event[2],migration_event[5:length(migration_event)])
						occupied[migration_event[2]] = occupied[migration_event[2]] + 1
					}
				} else if (sum(loc_sib1==ploc)==0 && sum(loc_sib2==ploc)>0){
					locations_sampled[siblings[2],1] = ploc
					if (samunif==1){ migration_event[2] = sample(loc_sib1,1) }
					else {
						migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib1, 0 )
						}
					locations_sampled[siblings[1],] = c(migration_event[2],migration_event[5:length(migration_event)])
					occupied[migration_event[2]] = occupied[migration_event[2]] + 1
				} else if (sum(loc_sib1==ploc)>0 && sum(loc_sib2==ploc)==0){
					locations_sampled[siblings[1],1] = ploc
					if (samunif==1){ migration_event[2] = sample(loc_sib2,1) }
					else {
						migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib2, 0 )
						}
					locations_sampled[siblings[2],] = c(migration_event[2],migration_event[5:length(migration_event)])
					occupied[migration_event[2]] = occupied[migration_event[2]] + 1
				}
				history_proba_val = history_proba_val + log(migration_event[3])
			}
		}
		# calculate LogLikelihood
		history_proba_val = log(1/ncol(space)) + history_proba_val
	}
	return(list(locations_sampled,history_proba_val,gtreel))
}

## simulate a sequence of migration events
sim_mig <- function(space, sigma, lambda, tau, timelimit, poplimit=0){
	occupied = vector('numeric', ncol(space))
	occupied[sample(1:ncol(space),1)] = 1 ## root
	history = matrix(0, 0, 4+dim(space)[1]) ## |dep|dest|proba|time2father|dist(|dist...)|
	time_elapsed = 0
	stop = 0
	while (stop==0){
		mig_event = new_migration(space, occupied, sigma, lambda, tau , 0 , 0 , 0)
		time_elapsed = time_elapsed + mig_event[4]
		history = rbind(history, mig_event, deparse.level=0)
		occupied[mig_event[2]] = occupied[mig_event[2]] + 1
		if (timelimit > 0){
			if (time_elapsed > timelimit) {
				if (nrow(history) > 1){
					history = history[1:(nrow(history)-1),] ## remove last migration that goes above the limit
				}
				if (is.vector(history)){ ## control, at least 2 migrations needed (one is removed in mig2treel())
					return(list(history, time_elapsed))
				}
				history[nrow(history),4] = history[nrow(history),4] + (timelimit-sum(history[,4]))
				stop=1
			}
		} else if (poplimit>0) {
		  if (sum(occupied)>=poplimit) {stop=1}
		} else {
			if (sum(occupied>0) == length(occupied)) {stop=1}
		}
	}
	if (nrow(history) == 1){
		history = as.vector(history)
		}
	return(list(history, time_elapsed))
}

##
## run MCMC for estimating variance and proba of colonising a previously occupied location
# freq: how often the estimated values are recorded
# Nstep_sigma: how many draws of sigma values at each iteration
# Nstep_lambda: how many draws of lambda at each iteration
# Nstep_step: how many draws of internal locations at each iteration
# Nstep_step: how many draws of genealogy at each iteration
# est_sigma: do not estimate sigma (0), use a uniform prior (1) or a dual prior (2)
mcmc_phyloland<-function(space, gtreel, simul_values, treelikelihood, est_sigma=c(1,1), est_lambda=1, est_tau=1, sample_geneal=0, sample_loc=0, 
                         plot_phylo=0, plot_MCMC=0, save_MCMC=0,
                         Nstep=1000, freq=1, Nstep_sigma=1, Nstep_lambda=1, Nstep_loc=1, Nstep_genealogy=1, Nstep_tau=1, pchange=0, ess_lim=100,
                         model=4, show_loc=0, filena_loc=NA, file_tracer=NA, file_tree=NA){
  
  if (!is.na(file_tracer)){
    write.table(t(c(0, rep(0, dim(space)[1]), 0, 0, 0)), file = file_tracer, row.names = FALSE, col.names = c("Sample", paste("Sigma", 1:dim(space)[1]), "Lambda", "Tau", "LogLikelihood"), sep = "\t", append = FALSE)
  }
  if (!is.na(file_tree)){
    #write("", file = file_tree, append = FALSE)
    file.create(file = file_tree)
  }
  if (!is.na(filena_loc)){
    file.create(file = filena_loc)
  }

  sigma_simul = simul_values[[1]]
  lambda_simul = simul_values[[2]]
  tau_simul = simul_values[[3]]
  locations_sim = simul_values[[4]]
  
  Ndim = dim(space)[1]
  tnumi = sum(is.na(gtreel$nodes[,2])) ## number of tips
  no_metropolis = 0 ## debug: to not use metropolis ratio (no data)
  
  mat_Dists = space_dist(space)
  
  if (sum(est_sigma)>0) sigma = rep(0,Ndim)
  sig_lower_limit = 0 # note that when new value of sigma proposed is very low then likelihood is -Inf then exp(metropolis+hastings)==0 then new value is rejected
  print("Compute sigma upper limit(s)...")
  sig_upper_limit = unlist(lapply(lapply(mat_Dists,max),sigma_upperlim)) # limit on prior for sampling sigma during mcmc
  sigma_threshold = unlist(lapply(lapply(mat_Dists,mean),sigma_upperlim)) # threshold below which limited dispersal is inferred
  if (sum(est_sigma==2)>0){
    sig_lam = rep(0,Ndim)
    sig_toplim = rep(0,Ndim)
    alpha = .5 ## probability for drawing prior on sigma equal to Uniform
  }
  for (nd in 1:Ndim){
    if (est_sigma[nd]>0){
      sigma[nd] = sig_upper_limit[nd]/2
      print(paste("sigma upper limit",nd,":",sig_upper_limit[nd]))
    }else{
      sigma[nd] = sigma_simul[nd]
    }
    if(est_sigma[nd]==2){
      sig_lam[nd] = -log(0.01)/sig_upper_limit[nd]
      sig_toplim[nd] = 1-exp(-sig_lam[nd]*sig_upper_limit[nd])
    }
  }
  
  lam_lower_limit = 1e-3
  lam_upper_limit = 1e3
  if (est_lambda==1){ ## start at random
    #lambda = runif(1,min=lam_lower_limit,max=lam_upper_limit)
    lambda = 1 ## start at 1, i.e. no competition
  }else{
    lambda = lambda_simul
  }
  
  tau_lower_limit = 1e-3
  tau_upper_limit = 1e3
  if (est_tau==1){ ## start at random
    #tau = runif(1,min=tau_lower_limit,max=tau_upper_limit)
    # start at average rate according to tree height and number of internal nodes
    tau = max(internal_node_height(gtreel)) / (tnumi-1)
  }else{
    tau = tau_simul
  }
  
  if (sample_loc==1 | length(sample_geneal[[1]])>1){ ## sample locations of all internal nodes, if sample_loc==0 the true internal nodes locations are used
    if (length(sample_geneal[[1]])>1){
      ind = sample(1:length(sample_geneal),1)
      tree_phylo = sample_geneal[[ind]]
      gtreel = converttree(tree_phylo)
      #gtreelikelihood_curr = treelikelihood[[ind]]
    }else{
      tree_phylo = lconverttree(gtreel)
    }
    real_loc = locations_sim[,1]
    locations = internal_loc(gtreel,locations_sim[,1]) ## get list of possible locations for all nodes, Uniform
    ## draw migrations uniformly
    locations_sampled = sim_history(gtreel,space,sig_upper_limit,1,tau_upper_limit/2,locations,0,model,1)[[1]]
    check_treel(gtreel,locations_sampled)
    if (show_loc==1){
      sort_space_id = sort(space,decreasing=TRUE,index.return=TRUE)$ix
      par(mfrow=c(2,1));
      xbar = barplot(hist(locations_sampled[,1],breaks=0:ncol(space),plot=FALSE)$counts[sort_space_id])
      title('true internal locations')
      histloc = matrix(nrow=(tnumi-1),ncol=Nstep) ## store location of internal nodes
      histdraw = vector("list",Nstep) ## store which internal node to be updated
    }
  }else{
    locations_sampled = locations_sim
    tree_phylo = lconverttree(gtreel)
  }
  
  ### acceptance rates computed on 100 values, first line values, second line acceptance
  acc_lambda = vector('numeric',100)
  id_acc_lambda = 1
  acc_sigma = matrix(0,Ndim,100)
  id_acc_sigma = vector('numeric',Ndim) + 1
  acc_tau = vector('numeric',100)
  id_acc_tau = 1
  acc_geneal = vector('numeric',100)
  id_acc_geneal = 1
  
  #### trees_phylo sampled
  trees_sampled = rmtree(1,2)
  list_trees = rmtree(1,2)
  list_trees = list_trees[-1]
  list_nb = vector("numeric",0)
  list_ind = vector("list",Nstep/freq)
  
  ### tuning parameters
  target_rate = 0.3 # acceptance rate to achieve
  init_tuningp = 1e2 ## start high and lower it progressively
  tuningp_lambda = init_tuningp
  tune_tp_lambda = 1 ## need to tune tuning parameter for lambda?
  tuningp_tau = init_tuningp
  tune_tp_tau = 1 ## need to tune tuning parameter for tau?
  tuningp_sigma = rep(0,Ndim) + init_tuningp
  tune_tp_sigma = rep(1,Ndim) ## need to tune tuning parameter for sigma?
  
  freq_n = 1
  freq_n1 = 1
  output_val = matrix(nrow=Nstep/freq,ncol=(4+Ndim))
  colnames(output_val) = c("Sample", paste("Sigma", 1:dim(space)[1]), "Lambda", "Tau", "LogLikelihood")
  output_loc = matrix(nrow=Nstep/freq,ncol=(length(gtreel$nodes[,1])+1))
  
  ## compute likelihood of this particular phylogeny
  likeli_curr = sim_history(gtreel,space,sigma,lambda,tau,locations_sampled[,1],1,model)[[2]]
  likeli_new = likeli_curr # initialisation
  
  # debug
  rate_acc = c(0,0,0)
  
  for (n in 1:Nstep) {
    
    ## sample internal location and/or genealogy
    if (length(sample_geneal[[1]])>1){
      for (n2 in 1:Nstep_genealogy) {
        ind = sample(1:length(sample_geneal),1)
        tree_phylo_new = sample_geneal[[ind]]
        gtreel_new = converttree(tree_phylo_new)
        ##gtreelikelihood_new = treelikelihood[[ind]]
        locations = internal_loc(gtreel_new,locations_sim[,1]) ## get list of possible locations for all nodes, Uniform
        locations_sampled_new = sim_history(gtreel_new,space,sigma,lambda,tau,locations,0,model,1)[[1]]
#ordered = reorder_treel(gtreel_new,locations_sampled_new[,1])
#locations_sampled_new = cbind( ordered[[2]], matrix(0,nrow=length(gtreel$nodes[,1]),ncol=((dim(space)[1]))) )
#gtreel_new = ordered[[1]]
        check_treel(gtreel_new,locations_sampled_new)
        likeli_new = sim_history(gtreel_new,space,sigma,lambda,tau,locations_sampled_new[,1],1,model)[[2]]
        ##metropolis = likeli_new + gtreelikelihood_new - likeli_curr - gtreelikelihood_curr
        metropolis = likeli_new - likeli_curr
        #accepted = 0
        if (runif(1,min=0,max=1)<exp(metropolis)){
          tree_phylo = tree_phylo_new
          gtreel = gtreel_new
          locations_sampled = locations_sampled_new
          likeli_curr = likeli_new
          ##gtreelikelihood_curr = gtreelikelihood_new
          #accepted = 1
#if (likeli_new>0){
#print(paste(ind,likeli_new))
#ordered = reorder_treel(gtreel,locations_sampled[,1])
#locations_ordered = ordered[[2]]
#tree_phylo_ordered = reorder_treep(tree_phylo)
#print(locations_ordered)
#print(paste(sigma[1],sigma[2],lambda,tau))
#plot_trees_tips(tree_phylo_ordered,locations_ordered,space)
#readline()
#}
        }
        if (plot_phylo==1){
          gtreeape = lconverttree(gtreel)
          plot.phylo(gtreeape,show.node.label=TRUE)
          readline()
        }
      }
    } else if(sample_loc==1) { ## sample location of internal nodes
      for (n2 in 1:Nstep_loc) {
        locations_sampled_new = sim_history(gtreel,space,sigma,lambda,tau,locations,0,model,1)[[1]]
        check_treel(gtreel,locations_sampled_new)
        likeli_new = sim_history(gtreel,space,sigma,lambda,tau,locations_sampled_new[,1],1,model)[[2]]
        metropolis = likeli_new - likeli_curr
        if (is.na(metropolis)) {
          if (is.infinite(likeli_new)) {}
          else{print("metropolis error (locations)")}
        } else {
          if (runif(1,min=0,max=1)<exp(metropolis)){
            locations_sampled = locations_sampled_new
            likeli_curr = likeli_new
          }
        }
      }
    }
    
    ## sample sigma (one per space dimension)
    for (nd in 1:Ndim){##**##
#nd=1##**##
      if (est_sigma[nd]>0){
        for (n2 in 1:Nstep_sigma) {
          if (est_sigma[nd]==1){
            sigma_new = sigma
            accepted = 0
            #print(sigma_new)
            sigma_new[nd] = sigma[nd] * exp(tuningp_sigma[nd]*(runif(1,min=0,max=1)-.5))
#sigma_new[2]=sigma_new[nd]##**##
            #print(sigma_new)
            if (sigma_new[nd]>sig_lower_limit && sigma_new[nd]<sig_upper_limit[nd]) {
              likeli_new = sim_history(gtreel,space,sigma_new,lambda,tau,locations_sampled[,1],1,model)[[2]]
              metropolis = likeli_new - likeli_curr
              if (no_metropolis==1) {metropolis = 0}
              hastings = log(sigma_new[nd]) - log(sigma[nd])
              if (is.na(metropolis)) {
                if (is.infinite(likeli_new)) {}
                else{print("metropolis error (sigma)")}
              } else if (is.na(hastings)) {
                print("hastings error (sigma)")
              } else {
                if (runif(1,min=0,max=1)<(exp(metropolis+hastings))) {
                  sigma = sigma_new
                  likeli_curr = likeli_new
                  accepted = 1
                }
              }
            }
          } else if (est_sigma[nd]==2){
            sigma_new = sigma
            accepted = 0
            if (runif(1,min=0,max=1)<alpha){
              sigma_new[nd] = sig_upper_limit[nd]
            }else{
              sigma_new[nd] = -log(1-runif(1,min=0,max=sig_toplim[nd]))/sig_lam[nd]
            }
            likeli_new = sim_history(gtreel,space,sigma_new,lambda,tau,locations_sampled[,1],1,model)[[2]]
            if (runif(1,min=0,max=1)<exp(likeli_new-likeli_curr)){
              sigma = sigma_new
              likeli_curr = likeli_new
              accepted = 1
            }
          }
          acc_sigma[nd,id_acc_sigma[nd]] = accepted
          id_acc_sigma[nd] = id_acc_sigma[nd] + 1
          if (id_acc_sigma[nd]==100){
            if (est_sigma[nd]==1 && tune_tp_sigma[nd]==1){ ## adjust tuning parameter
              accept_rate = sum(acc_sigma[nd,])/100
              if ( accept_rate<(target_rate*.5) ){
                tuningp_sigma[nd] = tuningp_sigma[nd] * 0.75
              } else if ( accept_rate<(target_rate*.9) ){
                tuningp_sigma[nd] = tuningp_sigma[nd] * 0.5
              } else if ( accept_rate<(target_rate*1.1) && accept_rate>(target_rate*.9) ){
                tune_tp_sigma[nd] = 0 ## stop tuning when reached approximately target_rate for the first time
                print(paste("tuning parameter sigma",nd,":",tuningp_sigma[nd]))
              } else if ( accept_rate>(target_rate*1.5) ){
                tuningp_sigma[nd] = tuningp_sigma[nd] * 1.5
              } else if ( accept_rate>(target_rate*1.1) ){
                tuningp_sigma[nd] = tuningp_sigma[nd] * 1.25
              }
            }
            id_acc_sigma[nd] = 1
          }
        }
        if (nd==2) rate_acc = rbind(rate_acc,c(sigma_new[nd],accepted,likeli_new)) # debug
      }
    }##**##
    
    ## sample lambda
    if (est_lambda==1){
      for (n2 in 1:Nstep_lambda) {
        accepted = 0
        lambda_new = lambda * exp(tuningp_lambda*(runif(1,min=0,max=1)-.5))
        if (lambda_new>lam_lower_limit && lambda_new<lam_upper_limit) {
          likeli_new = sim_history(gtreel,space,sigma,lambda_new,tau,locations_sampled[,1],1,model)[[2]]
          metropolis = likeli_new - likeli_curr
          if (no_metropolis==1) {metropolis = 0}
          hastings = log(lambda_new) - log(lambda)
          if (is.na(metropolis)) {
            if (is.infinite(likeli_new)) {}
            else{print("metropolis error (lambda)")}
          } else if (is.na(hastings)) {
            print("hastings error (lambda)")
          } else {
            if (runif(1,min=0,max=1)<(exp(metropolis+hastings))) {
              lambda = lambda_new
              likeli_curr = likeli_new
              accepted = 1
            }
          }
        }
        acc_lambda[id_acc_lambda] = accepted
        id_acc_lambda = id_acc_lambda + 1
        if (id_acc_lambda==100){
          if (tune_tp_lambda==1){
            accept_rate = sum(acc_lambda)/100
            if ( accept_rate<(target_rate*.5) ){
              tuningp_lambda = tuningp_lambda * 0.75
            } else if ( accept_rate<(target_rate*.9) ){
              tuningp_lambda = tuningp_lambda * 0.5
            } else if ( accept_rate<(target_rate*1.1) && accept_rate>(target_rate*.9) ){
              tune_tp_lambda = 0 ## stop tuning when reached approximately 0.3 for the first time
              print(paste("tuning parameter lambda :",tuningp_lambda))
            } else if ( accept_rate>(target_rate*1.5) ){
              tuningp_lambda = tuningp_lambda * 1.5
            } else if ( accept_rate>(target_rate*1.1) ){
              tuningp_lambda = tuningp_lambda * 1.25
            }
            id_acc_lambda = 1
          }
        }
      }
    }
    
    ## sample tau
    if (est_tau==1){
      for (n2 in 1:Nstep_tau) {
        accepted = 0
        tau_new = tau * exp(tuningp_tau*(runif(1,min=0,max=1)-.5))
        if (tau_new>tau_lower_limit && tau_new<tau_upper_limit) {
          likeli_new = sim_history(gtreel,space,sigma,lambda,tau_new,locations_sampled[,1],1,model)[[2]]
          metropolis = likeli_new - likeli_curr
          if (no_metropolis==1) {metropolis = 0}
          hastings = log(tau_new) - log(tau)
          if (is.na(metropolis)) {
            if (is.infinite(likeli_new)) {}
            else{print("metropolis error (tau)")}
          } else if (is.na(hastings)) {
            print("hastings error (tau)")
          } else {
            if (runif(1,min=0,max=1)<(exp(metropolis+hastings))) {
              tau = tau_new
              likeli_curr = likeli_new
              accepted = 1
            }
          }
        }
        acc_tau[id_acc_tau] = accepted
        id_acc_tau = id_acc_tau + 1
        if (id_acc_tau==100){
          if (tune_tp_tau==1){
            accept_rate = sum(acc_tau)/100
            if ( accept_rate<(target_rate*.5) ){
              tuningp_tau = tuningp_tau * 0.5
            } else if ( accept_rate<(target_rate*.9) ){
              tuningp_tau = tuningp_tau * 0.75
            } else if ( accept_rate<(target_rate*1.1) && accept_rate>(target_rate*.9) ){
              tune_tp_tau = 0 ## stop tuning when reached approximately 0.3 for the first time
              print(paste("tuning parameter tau :",tuningp_tau))
            } else if ( accept_rate>(target_rate*1.5) ){
              tuningp_tau = tuningp_tau * 1.5
            } else if ( accept_rate>(target_rate*1.1) ){
              tuningp_tau = tuningp_tau * 1.25
            }
            id_acc_tau = 1
          }
        }
      }
    }
    
    ## record value
    if (n==(freq_n*freq)) {
      ## save tree, reorder before so that internal nodes are from oldest to most recent so that location_sampled is valid for both tree formats
      ordered = reorder_treel(gtreel,locations_sampled[,1])
      locations_ordered = ordered[[2]]
      tree_phylo_ordered = reorder_treep(tree_phylo)
      trees_sampled[[n/freq]] = tree_phylo # _ordered
      test <- test_tree(tree_phylo_ordered, list_trees, list_nb, list_ind, n ,freq)
      list_trees = test[[1]]
      list_nb = test[[2]]
      list_ind = test[[3]]
      if (!is.na(file_tree)){
        write(write.tree(tree_phylo_ordered,file=""), file = file_tree, append = TRUE)
      }
      freq_n = freq_n + 1
      output_val[n/freq,1] = n
      
      ## save parameters
      for (idn in 1:Ndim){
        output_val[n/freq,(idn+1)] = sigma[idn]
      }
      output_val[n/freq,(idn+2)] = lambda
      output_val[n/freq,(idn+3)] = tau
      output_val[n/freq,(idn+4)] = likeli_curr
      stringo = sprintf("%i    ",output_val[n/freq,1])
      for (idn in 2:length(output_val[1,])){
        stringo = paste(stringo,sprintf("%.4f    ",output_val[n/freq,idn]))
      }
      print(stringo)
      
      ## save locations
      output_loc[n/freq,1] = n
      #output_loc[n/freq,2:(length(gtreel[,1])+1)] = locations_sampled[,1]
      output_loc[n/freq,2:(length(gtreel$nodes[,1])+1)] = locations_ordered
      if (!is.na(filena_loc)){
        write.table(t(c(n,locations_ordered)), col.names = FALSE, row.names = FALSE, file = filena_loc, sep = "\t", eol = "\n", append = TRUE)
      }
      #check_treel(reorder_treel(converttree(trees_sampled[[n/freq]]))[[1]], output_loc[n/freq,2:(length(gtreel$nodes[,1])+1)])
      if (!is.na(file_tracer)){
        write.table(t(output_val[n/freq,]), file = file_tracer, row.names = FALSE, col.names = FALSE, sep = "\t", append = TRUE)
      }
      
      ## check effective sample size
      if (n==(10*freq_n1*freq)){
        freq_n1 = freq_n1 + 1
        ess_val = vector('numeric',(Ndim+2)) ## sigmas plus lambda and tau
        ness = 1
        for (idn in 2:(Ndim+3)){
          val = output_val[is.finite(output_val[,idn]),idn]
          ess_val[ness] = effSaSize(val)
          for (nd in 1:Ndim){
            if (est_sigma[nd] == 0){
              ess_val[nd] = ess_lim	
            }
          }
          if (est_lambda == 0){
            ess_val[3] = ess_lim	
          }
          if (est_tau == 0){
            ess_val[4] = ess_lim	
          }
          ness = ness+1
        }
        stringo1 = sprintf("            ")
        stringo2 = sprintf("ess:    ")
        for (idn in 1:Ndim){
          stringo1 = paste(stringo1, sprintf("Sigma       "))
          stringo2 = paste(stringo2, sprintf("%.4f    ", ess_val[idn]))
        }
        stringo1 = paste(stringo1, sprintf("Lambda      Tau"))
        stringo2 = paste(stringo2, sprintf("%.4f    %.4f", ess_val[idn+1], ess_val[idn+2]))
        print(stringo1)
        print(stringo2)
        #stringo2 = paste(stringo2,sprintf("Lambda    Tau    LogLikelihood"))
        #print(paste(rep('Sigma',Ndim),'Lambda','Tau','LogLikelihood'))
        #print(paste('ess:',toString(ess_val)))
        if (sum(ess_val<ess_lim)==0){
          L = sum(is.finite(output_val[,2]))
          output_val = output_val[1:L,]
          output_loc = output_loc[1:L,]
          break
        }
      }
    }
  }
  if (sample_loc==1 && show_loc==1){
    return(list(output_val,histloc,histdraw,tnumi,real_loc,locations_sampled,space))
  }else{
    return(list(output_val, output_loc, trees_sampled, list_trees, list_nb, list_ind[1:length(list_trees)], rate_acc, sigma_threshold))
  }
}

## migration in C
new_migration <- function(space, occupied, sigma, lambda, tau, departure, destination, mig){
  migr <- function(space, occupied, sigma, lambda, tau, mig, mig_event, space_dim, space_size, length_mig){
    .C("migC", space, occupied, sigma, lambda, tau, mig, mig_event, space_dim, space_size, length_mig, PACKAGE = "phyloland")}
  
  space_dimC = as.integer(nrow(space))
  space_sizeC = as.integer(ncol(space))
  
  spaceC = as.double(as.vector(t(space)))
  occupiedC = as.integer(occupied)
  sigmaC  = as.double(sigma)
  lambdaC = as.double(lambda)
  tauC = as.double(tau)
  migC = as.double(mig)
  
  mig_eventC = as.double(vector("numeric", (4+space_dimC)))
  length_migC = as.integer(length(migC))
  
  if (length(sigma) != space_dimC) stop('new_migration error 1')
  
  mig_event = migr(spaceC, occupiedC, sigmaC, lambdaC, tauC, migC, mig_eventC, space_dimC, space_sizeC, length_migC)[[7]]
  
  return(mig_event)
}

## compute effective sample size as N/(1+2*sum(autocorr))
effSaSize <- function(x){
  if (sum(x-x[1]) == 0) return(0) ## return 0 if all elements of x are equals
  autocor = acf(x,type="correlation",plot=FALSE,lag.max=length(x))[[1]]
  cutoff = which(autocor<(2*sd(x)))[1]
  if (is.na(cutoff)){
    return( NA )
  }else{
    return( length(x)/(1+2*sum(autocor[2:cutoff])) )
  }
}

## 
stat_phyloland <- function(mcmco, plot_MCMC = 0, save_MCMC = 0){
    output_val = mcmco
    smod = vector('numeric',dim(mcmco)[2]-4)
    smean = vector('numeric',dim(mcmco)[2]-4)
    squant = matrix(0,(dim(mcmco)[2]-4),7)
    slast = vector('numeric',dim(mcmco)[2]-4)
    for (n in 1:(dim(mcmco)[2]-4)){
      id = n+1
      samplevals = mcmco[,id]
      if (length(samplevals)>1){ 
        den = density(samplevals)
        smod[n] = den$x[which(den$y==max(den$y))[1]]
        squant[n,] = quantile(samplevals,probs=c(0,.05,.25,.5,.75,.95,1),names=FALSE)
        smean[n] = mean(samplevals)
        slast[n] = samplevals[length(samplevals)]
      }else if (length(samplevals)==1){
        smod[n] = samplevals
        squant[n,] = quantile(samplevals,probs=c(0,.05,.25,.5,.75,.95,1),names=FALSE)
        smean[n] = samplevals
        slast[n] = samplevals
      }
    }
    den = density(mcmco[,(dim(mcmco)[2]-2)])
    lmod = den$x[which(den$y==max(den$y))[1]]
    lquant = quantile(mcmco[,(dim(mcmco)[2]-2)],probs=c(0,.05,.25,.5,.75,.95,1),names=FALSE)
    lmean = mean(mcmco[,(dim(mcmco)[2]-2)])
    llast = mcmco[nrow(mcmco),(dim(mcmco)[2]-2)]
    den = density(mcmco[,(dim(mcmco)[2]-1)])
    tmod = den$x[which(den$y==max(den$y))[1]]
    tquant = quantile(mcmco[,(dim(mcmco)[2]-1)],probs=c(0,.05,.25,.5,.75,.95,1),names=FALSE)
    tmean = mean(mcmco[,(dim(mcmco)[2]-1)])
    tlast = mcmco[nrow(mcmco),(dim(mcmco)[2]-1)]
    
    for (n in 1:(dim(mcmco)[2]-4)){
      if (n==1){ strout = paste(toString(round(squant[n,],3)),round(smod[n],3),round(smean[n],3),round(slast[n],3),sep=",") }
      else if (n>1){ strout = paste(strout,toString(round(squant[n,],3)),round(smod[n],3),round(smean[n],3),round(slast[n],3),sep=",") }
    }
    strout = paste(strout,toString(round(lquant,3)),round(lmod,3),round(lmean,3),round(llast,3),
                   toString(round(tquant,3)),round(tmod,3),round(tmean,3),round(tlast,3),sep=",")
    return(strout)
  }

# read space file : a text file with 3 columns (separated by tabs). The first one with the tips names (the same as in \code{fileTrees}), the second one with the location latitudes and the last one with the location longitudes (both in decimal degrees). No header.
# return the tips locations, the space with latitudes and longitudes, the space [0,1]
read_space <- function(fileDATA, tips){
	sampled_loc <- as.matrix(read.table(fileDATA, header = FALSE, sep="\t"))
  if (length(sampled_loc)<(3*length(readLines(fileDATA)))){ # parsing failed, try sep=" " 
    sampled_loc <- as.matrix(read.table(fileDATA, header = FALSE, sep=" "))
  }
	tips_file = as.character(sampled_loc[,1])	
	if (sum(sort(tips_file) == sort(tips)) != length(tips)){
		print(sort(tips_file)[which(sort(tips_file) != sort(tips))])
		stop("in interface : tips names in fileDATA not in fileTREES",  call.=FALSE)	
	}
	coord_tips = matrix(c(as.double(sampled_loc[,2]),as.double(sampled_loc[,3])),ncol = 2)
	loc_tips <- vector("numeric",dim(coord_tips)[1])
	space = t(coord_tips[-which(duplicated(coord_tips[,1]) & duplicated(coord_tips[,2])),])
	for (i in 1:length(loc_tips)){
		loc_tips[which(tips==tips_file[i])]  =  which(space[1, ] == coord_tips[i, 1] & space[2, ] == coord_tips[i, 2])
	}
	space_dec = space
  # normalise space in 0-1 in each dimension
	for ( i in 1:dim(space)[1]){
		L = space[i,]
		if (sum(L < 0) != 0){
			abso = abs(min(L))
			for (j in 1:length(L)){
				L[j] = L[j] + abso	
			}
		}
		mL = max(L) - min(L)
		space[i,] = (L - min(L))/mL	
	}	
	return(list(loc_tips,space,space_dec))	
}

# read locations file (output interface)
read_loc <- function(fileLOC){
	sampled_loc <- read.table(file=fileLOC, header=FALSE, sep="\t", quote='"', comment.char="", check.names=FALSE)
	sampled_loc = sampled_loc[,-1]
	return(sampled_loc)
}

# verif function arguments : type, length...
verif_arg <- function(x, type = NA, len = NA, name_function, positive = 0){
	if (sum(is.na(type)) == 0){
		if(!(is.element(class(x),type))){
			stop("in ",name_function ," : invalid 'class' (", class(x), ") for '", deparse(substitute(x)),"'", call. = FALSE)
		} 
	}
	if (sum(is.na(len)) == 0){
		if(!(is.element(length(x),len))){
			stop("in ",name_function ," : parameter '",deparse(substitute(x)),"' is of incorrect size", call. = FALSE)
		}		
	}
	if (positive == 1){
		if(sum(x<0) >0){
			stop("in ",name_function ," : invalid negative value specified for parameter '",deparse(substitute(x)),"'. ","'",deparse(substitute(x)),"' must be positive", call. = FALSE)
		}
	}else if (positive == 2){
		if(sum(x<=0) >0){
			stop("in ",name_function ," : invalid value specified for parameter '",deparse(substitute(x)),"'. ","'",deparse(substitute(x)),"' must be strictly positive", call. = FALSE)
		}	
	}
}

# plot tree with location (print tip names)
# plot_trees_tips0<-function(tree_phylo, location){
# 	tips = tree_phylo$tip.label
# 	gtreel = converttree(tree_phylo)
# 	ind_tips = vector("numeric",length(tips))
# 	for (i in 1:length(tips)){
# 		ind_tips[i] = paste(i, '-', tips[i], sep = "")
# 	}
# 	tree_phylo = lconverttree(gtreel)
# 	tree_phylo$tip.label = ind_tips
# 	tree_par = write.tree(tree_phylo)
# 	tree_par = gsubfn('t([0-9]*):', function(x) paste(as.numeric(x), '-L',location[as.numeric(x)], ':', sep=''), tree_par, engine = "R")
# 	tree_par = gsubfn('t([0-9]*);', function(x) paste(as.numeric(x), '-L',location[as.numeric(x)], ';', sep=''), tree_par, engine = "R")
# 	for (i in tips){
# 		tree_par = replace_tips(i, tree_par, tips, location)
# 	}
# 	plot.phylo(read.tree(file="", text=tree_par), show.tip.label=TRUE, show.node.label=TRUE)
# 	print(tree_par)	
# }

# replace_tips <- function(tip, tree_par, tips, location){
# 	tree_par = gsubfn(paste(which(tips==tip),'-',tip,sep=""), paste(which(tips==tip),'-',tip,'-L',location[which(tips == tip)],sep =''), tree_par, engine = "R")
# 	return(tree_par)
# }

# plot tree with location names for internal nodes and tips
# "location" contains the index of location from the oldest up to the most recent internal nodes, therefore needs to make sure tree_phylo is ordered the same way
plot_trees_tips<-function(tree_phylo, location, space, show_index=0){
  # make sure the internal nodes indexes are in correct order with location
  ordered = reorder_treel(converttree(tree_phylo),location)
  #location = ordered[[2]]
  tree_phylo = reorder_treep(tree_phylo)
  if (length(colnames(space))>0){ # locations of space are named
    nodes = colnames(space)
  }
#   if (show_index==1){
#     tips = paste(tree_phylo$tip.label, colnames(location)[as.numeric(sub("t","",tree_phylo$tip.label))], sep="_")
#     tree_phylo$node.label = paste(as.numeric(location[(length(tree_phylo$tip.label)+1):(length(location))]), nodes[ as.numeric(location[(length(tree_phylo$tip.label)+1):(length(location))]) ], sep="_")
#   }else{
#     tips = colnames(location)[as.numeric(sub("t","",tree_phylo$tip.label))]
#     tree_phylo$node.label = nodes[ as.numeric(location[(length(tree_phylo$tip.label)+1):(length(location))]) ]
#   }
  tree_phylo$node.label = colnames(space)[ as.numeric(location[(Ntip(tree_phylo)+1):(2*Ntip(tree_phylo)-1)]) ]
#  tree_phylo$tip.label = tips
  tree_par = write.tree(tree_phylo)
  plot.phylo(read.tree(file="", text=tree_par), show.tip.label=TRUE, show.node.label=TRUE, cex=.7, label.offset=.001)
  #plot.phylo(read.tree(file="", text=tree_par), show.tip.label=TRUE, show.node.label=TRUE, cex=.7) # Ape's label.offset option seems unstable??
  print(tree_par)	
}

# Return list of all the ancestors of a node
ancestors<-function(node,gtreel){
	res=vector("numeric",0)
	ind<-node
	is_root=is.na(gtreel$nodes[ind,1])
	if (is.na(gtreel$nodes[node,1])){
		is_root=TRUE
	}
	while (is_root==FALSE){
		ind=gtreel$nodes[ind,1]
		res[length(res)+1]=ind
		is_root=is.na(gtreel$nodes[ind,1])
	}
	return(res)
}
	
# Return first common element to lists
ancestor<-function(lists){
	res = lists[[1]]
	for (i in lists[2:length(lists)]) {
			res = intersect(res, i)
	}
	return(res[1])
}

# Return location of a common ancestor
loc_ancestor<-function(gtreel,tips,sampled_loc){
	anc=ancestor(lapply(X=tips,FUN=ancestors,gtreel))
	return(as.integer(sampled_loc[anc]))
}

# Checks if the tree is contained in the list of trees already sampled
test_tree <- function(tree_phylo, list_trees, list_nb, list_ind, n, freq){
	equal <- lapply(list_trees, "all.equal.phylo", current = tree_phylo, use.edge.length = FALSE)
	if (length(which(equal == TRUE)) > 0){
		list_nb[which(equal == TRUE)] = list_nb[which(equal == TRUE)] + 1
		list_ind[[which(equal == TRUE)]][length(list_ind[[which(equal == TRUE)]]) + 1] = n/freq
	}else{
		list_trees[[length(list_trees) + 1]] <- tree_phylo
		list_nb[length(list_nb) + 1] = 1	
		for ( i in 1:length(list_ind) ){
			if (length(list_ind[[i]]) == 0){
				list_ind[[i]][1] = n/freq
				break	
			}	
		}
	}
	return(list(list_trees, list_nb, list_ind))
}

# likelihood extraction
treeLikeli <- function(line, pattern){
  if(regexpr(pattern,line)[1] < 0 ){
    stop(" in interface : wrong 'pattern_trees_likelihood'",call.=FALSE)
  }
  start = regexpr(pattern,line)[1] + nchar(pattern)
  stop = regexpr("]",line)[1] -1
  Lchar = strsplit(line, split = "")[[1]][start:stop]
  return(as.numeric(paste(Lchar, collapse="")))
}

# compute pairwise distances of space
space_dist <- function(space, dmethod="euclidean"){
  mat_Dists = list()
  for (n in 1:dim(space)[1]){
    mat_Dists[[n]] = dist(space[n,], diag=T, upper=T, method=dmethod)
  }
  return(mat_Dists)
}

# find upper limit for sigma so that area of normal density represent 95% of the uniform U(0,dmax)
#sigma_upperlim <- function(dmax,aire=.95,n_itera=100,sigma_min=1e-3,sigma_max=1e3,tries=1e4,plotting=0){
sigma_upperlim <- function(dmax,aire=.95,n_itera=10,sigma_min=1e-3,sigma_max=1e3,tries=10,plotting=0){
  for (i in 1:n_itera){
    testsigma = seq(sigma_min,sigma_max,length.out=tries)
    proba1 = matrix(0,nrow=length(testsigma),ncol=2)
    for(i in 1:length(testsigma)){
      proba1[i,1] = testsigma[i]
      norma_factor = pnorm(dmax,0,sqrt(testsigma[i]))-pnorm(0,0,sqrt(testsigma[i]))
      x = sqrt( -2 * testsigma[i] * log(norma_factor * (1/dmax) * sqrt(2*pi) * sqrt(testsigma[i])) )
      proba1[i,2] = x*(1/dmax) + (pnorm(dmax,0,sqrt(testsigma[i]))-pnorm(x,0,sqrt(testsigma[i])))/norma_factor
    }
    if (proba1[length(testsigma),2]<aire) {
      sigma_min = sigma_max
      sigma_max = sigma_max^2
    }else if (proba1[1,2]>aire) {
      sigma_max = sigma_min
      sigma_min = sigma_min/2
    }else{
      # upper limit
      n = which(proba1[,2]>=aire)[1]
      sigma_min = proba1[n-1,1]
      sigma_max = proba1[n,1]
    }
  }
  # upper limit
  sig_upper_limit = testsigma[which(proba1[,2]>=aire)][1]
  return(sig_upper_limit)
}

##check the structure of a gtreel phylogenetic tree
check_treel <- function(gtreel,locations=NA){
  error_found = 0
  if ( length(which(is.na(gtreel$nodes[,1])))!=1 ){
    print("error root")
    error_found = 1
  }
  for (n in 1:length(gtreel$nodes[,1])){
    pere = gtreel$nodes[n,1]
    filsg = gtreel$nodes[n,2]
    filsd = gtreel$nodes[n,3]
    if ( !is.na(pere) ){
      if ( gtreel$nodes[pere,2]!=n && gtreel$nodes[pere,3]!=n ){
        print("error 1")
        error_found = 1
      }
    }
    if (  (is.na(filsg) & !is.na(filsd))  |  (!is.na(filsg) & is.na(filsd))  ){
      print("error 2")
      error_found = 1
    }
    if ( !is.na(filsg) ){
      if ( gtreel$nodes[filsg,1]!=n ){
        print("error filsg")
        error_found = 1
      }
    }
    if ( !is.na(filsd) ){
      if ( gtreel$nodes[filsd,1]!=n ){
        print("error filsd")
        error_found = 1
      }
    }
    if (!is.na(locations[1])){
      if ( !is.na(filsg) & !is.na(filsd) ){
        if ( locations[filsg]!=locations[n] & locations[filsd]!=locations[n] ){
          print(paste("error locations, node",n))
          error_found = 1}
      }
    }
  }
  return(error_found)
}

# reOrder tree format louis, so that oldest internal node is first after the root
# tips must appear first, then root and finally internal nodes (case when gtreel is output by converttree())
reorder_treel <- function(gtreel,locations=0){
  #gtreelO = matrix(0,nrow=nrow(gtreel$nodes),ncol=ncol(gtreel$nodes))
  gtreelO <- list( nodes=matrix(0,nrow=nrow(gtreel$nodes),ncol=ncol(gtreel$nodes)), nodes.label=seq(1,length(gtreel$nodes.label)) )
  locationsO = locations*0 # intialise to 0
  ntips = ( length(gtreel$nodes[,1])+1 )/2
  nnodes = ntips-1
  if (!is.na(gtreel$nodes[ntips+1,1])){stop('reorder_treel() error: root not following tips')}
  gtreelO$nodes[1:ntips,2:3] = cbind(rep(NA,ntips),rep(NA,ntips)) # tips first (no children)
  gtreelO$nodes[ntips+1,c(1,4)] = NA # then root (no parent)
  gtreelO$nodes[1:ntips,4] = gtreel$nodes[1:ntips,4] # height for tips stay the same
  ## get the absolute height of each internal node
  inode_height = internal_node_height(gtreel)
  ## need to treat the nodes from the most recent to the oldest 
  node_id = sort(inode_height,decreasing=FALSE,index.return=TRUE)$ix
  rootn = which(is.na(gtreel$nodes[,1]))
  if (node_id[length(node_id)]!=rootn){stop('reorder_treel() error: root')}
  nh = cbind(node_id,c(node_id[1:ntips],nrow(gtreel$nodes):(rootn+1),rootn))
  locationsO[nh[,2]] = locations[nh[,1]]
  for (n in (ntips+1):nrow(gtreel$nodes)) {
    if(gtreel$nodes[nh[n,1],2]<=ntips){
      filsg = gtreel$nodes[nh[n,1],2]
    }else{
      filsg = nh[which(nh[,1]==gtreel$nodes[nh[n,1],2]),2]
    }
    if(gtreel$nodes[nh[n,1],3]<=ntips){
      filsd = gtreel$nodes[nh[n,1],3]
    }else{
      filsd = nh[which(nh[,1]==gtreel$nodes[nh[n,1],3]),2]
    }
    gtreelO$nodes[nh[n,2],2:3] = c(filsg,filsd)
    gtreelO$nodes[filsg,1] = nh[n,2]
    gtreelO$nodes[filsd,1] = nh[n,2]
    gtreelO$nodes[filsg,4] = gtreel$nodes[gtreel$nodes[nh[n,1],2],4]
    gtreelO$nodes[filsd,4] = gtreel$nodes[gtreel$nodes[nh[n,1],3],4]
    gtreelO$nodes.label[filsg] = gtreel$nodes.label[gtreel$nodes[nh[n,1],2]]
    gtreelO$nodes.label[filsd] = gtreel$nodes.label[gtreel$nodes[nh[n,1],3]]
  }
  return(list(gtreelO,locationsO))
}

# reOrder tree format Ape, so that oldest internal node appear first
# in Ape tree format, tips of the tree are always from 1 to number of tips of tree_phylo and internal nodes are indexed from Ntip(tree_phylo)+1
reorder_treep <- function(tree_phylo){
  tree_phyloO = tree_phylo
  ntips = Ntip(tree_phylo)
  edges = tree_phylo$edge
  edgesO = edges #rep(0,2*ntips-1)
  ## get the depth of each node
  nod_depth = node.depth.edgelength(tree_phylo)
  ## need to treat the internal nodes from the most recent to the oldest 
  nod_depth = nod_depth[(ntips+1):length(nod_depth)]
  node_id = sort(nod_depth,decreasing=TRUE,index.return=TRUE)$ix
  node_id = c(1:ntips,node_id+ntips)
  nh = cbind(node_id,c(node_id[1:ntips],(2*ntips-1):(ntips+1)))
  for (n in 1:length(nh[,1])){
    edgesO[which(edges==nh[n,1])] = nh[n,2]
  }
  tree_phyloO$edge = edgesO
  return(tree_phyloO)
}

# look for node of equal heights in the tree, if found, break the try by randomly changing the corresponding node heights
# write a set of trees if numrep>1
break_node_height_tie <- function(gtreel, tiplabel, numrep=100){
  inode_height = internal_node_height(gtreel)
  if ( sum( duplicated(inode_height[which(inode_height>0)]) )==0 ) {stop('no pair of internal nodes with equal heights')} # no internal nodes have identical heights
  node_id = sort(inode_height,decreasing=TRUE,index.return=TRUE)$ix # treat nodes from oldest to most recent
  treeset <- vector("list", numrep)
  class(treeset) <- "multiPhylo"
  for (nrep in 1:numrep){
#     print(nrep)
    gtreelN = gtreel
    n=1
    while (inode_height[node_id[n]]>0) { # stop when reached the first tip
      #     x11(); plotree(gtreel,seq(1,15))
      m=n+1
      if (inode_height[node_id[n]]==inode_height[node_id[m]]) { # tie found need to break it
#         print(paste("node:",node_id[n]))
        L_to_break=inode_height[node_id[n]]
        maxL=inode_height[node_id[n-1]]
        while (inode_height[node_id[m]]==L_to_break) {
          m=m+1
          minL=inode_height[node_id[m]]
        }
        #       print(c(node_id[n],maxL,minL))
        for (p in n:(m-1)){
          #         print(node_id[p])
          # sample a truncated normal between the height of previous and following nodes: Compute the CDF for the upper and lower limits of the interval and 
          # generate a uniform random numbers within that range of CDF values, then compute the inverse CDF
          # standard deviation is set as 5% of interval of height
          sdev=0.05*(maxL-minL)
          new_height = qnorm(runif(1, pnorm(minL, mean=L_to_break, sd=sdev), pnorm(maxL, mean=L_to_break, sd=sdev)), mean=L_to_break, sd=sdev)
          height_correction = new_height - inode_height[node_id[p]]
#           print(height_correction)
          gtreelN$nodes[node_id[p],4] = gtreelN$nodes[node_id[p],4] + (height_correction)
          gtreelN$nodes[gtreelN$nodes[node_id[p],2:3],4] = gtreelN$nodes[gtreelN$nodes[node_id[p],2:3],4] - (height_correction)
        }
      }
      n=m
    }
#     x11(); plotree(gtreelN,seq(1,15))
    treeset[[nrep]] = lconverttree(gtreelN)
  }
  return(treeset)
}
