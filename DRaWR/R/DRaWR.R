#' threeCol2MaxMat
#'
#' This function takes a three vectors of equal length (source nodes, target nodes, and edge weights) and return the adjacency matrix as a sparse Matrix.
#' @param a (vector): vector of source node names.
#' @param b (vector): vector of target node names.
#' @param v (vector): vector of edge weights names.
#' @return sparce Matrix.
#' @keywords sparce matrix
#' @import Matrix
#' @export
#' @examples
#' threeCol2MaxMat(a = c("a","b","c","c"), b = c("a","b","b","b"), v = c(1,2,3,4))
threeCol2MaxMat<- function(a = c("a","b","c","c"), b = c("a","b","b","b"), v = c(1,2,3,4)){
#	library(Matrix)
	if(length(a)!=length(b) | length(a)!=length(v)){
		return(-1)
	}
	a = as.character(a)
	b = as.character(b)
	v = as.numeric(v)
	avals = sort(unique(a))
	bvals = sort(unique(b))

	nrows = length(avals)
	ncols = length(bvals)

	aidxs = match(a, avals)
	bidxs = match(b, bvals)

	vidxs = (bidxs-1)*nrows+aidxs

	retMat = Matrix(0, nrows, ncols, dimnames = list(avals,bvals), sparse=T)
	retMat[vidxs] = v
# arbitrarily selects one value if multiple values for same edge
	return(retMat)
}

#' threeCol2listMat
#'
#' This function takes a three vectors of equal length (source nodes, target nodes, and edge weights) and return the adjacency matrix as a list of vectors.
#' @param a (vector): vector of source node names.
#' @param b (vector): vector of target node names.
#' @param v (vector): vector of edge weights names.
#' @return list of vectors matrix representation.
#' @keywords sparce matrix
#' @export
#' @examples
#' threeCol2listMat(a = c("a","b","c","c"), b = c("a","b","b","b"), v = c(1,2,3,4))
threeCol2listMat<- function(a = c("a","b","c","c"), b = c("a","b","b","b"), v = c(1,2,3,4)){
	if(length(a)!=length(b) | length(a)!=length(v)){
		return(-1)
	}
	a1 = as.character(a)
	b1 = as.character(b)
	v1 = as.numeric(v)

	show("sorting...")
	ss = sort(a1, index.return=T)
	a2 = a1[ss$ix]
	b2 = b1[ss$ix]
	v2 = v1[ss$ix]

	avals = unique(a2)
	bvals = sort(unique(b2))
#	nvals = 10000  # for testing limit
	nvals = length(a2)
	nrows = length(avals)
	ncols = length(bvals)

	listm = list()
	colsums = 1:length(bvals)*0
	show("creating list...")
	curval = a2[1]
	idxstart = 1
	for (i in 2:nvals){
		if(i%%100000==1){show(i)}

		if(a2[i] != curval){
			# process values for row
			cols = b2[idxstart:(i-1)]
			vals = v2[idxstart:(i-1)]
# arbitrarily selects one value if multiple values for same edge
			listidxs = match(unique(cols), cols)
			colidxs = match(cols[listidxs],bvals)
			listm[[as.character(curval)]] = list("colidxs" = colidxs,"vals" = vals[listidxs])
			colsums[colidxs] = colsums[colidxs] + vals[listidxs]
			# update pointers
			curval = a2[i]
			idxstart = i
		}
	}
	# process values for end
	if(a2[i]==a2[i-1]){
		cols = b2[idxstart:i]
		vals = v2[idxstart:i]
# arbitrarily select one value if multiple values for same edge
		listidxs = match(unique(cols), cols)
		colidxs = match(cols[listidxs],bvals)
		listm[[as.character(curval)]] = list("colidxs" = colidxs,"vals" = vals[listidxs])
		colsums[colidxs] = colsums[colidxs] + vals[listidxs]
	}
	return(list("listm" = listm, "avals" = avals, "colsums"=colsums))
}


#' RWR
#'
#' This function runs a random walk with restart using two supported matrix representations.
#' @param boolSparceMat (bool): Boolean to indicate sparce Matrix or list matrix.
#' @param transmat (sparce Matrix / list matrix): transition probabilites.
#' @param restart (float): probability of restart.
#' @param query (vector): probability of restarting at all nodes.
#' @param startvec (vector): initial probability of being at any node.
#' @param maxiters (int): maximum number of allowable iterations.
#' @param thresh (float): threshold for L1 norm convergence.
#' @return list of 'iter':number of iterations, 'diff': L1 norm of difference, 'vec': converged probability distribution vector.
#' @keywords random walk with restart
#' @import Matrix
#' @export
#' @examples
#' RWR(boolSparceMat=TRUE, transmat=transmat, restart=.3, query=c(rep(0.1,10),rep(0,5)),
#'		startvec=rep(1/15,15), maxiters=10, thresh=0.001)
RWR<- function(boolSparceMat, transmat, restart, query, startvec, maxiters, thresh){
#	library(Matrix)
	damping = 1-restart
	query = query / sum(query)
	vec = startvec
	mycolsum = NULL
	mylistm = NULL
	mismatch = 0

	if(!boolSparceMat) {
		mycolsums = transmat$colsums
		mylistm = transmat$listm
		mismatch = sum(names(query)!=names(mylistm))
	} else {
		mismatch = sum(names(query)!=colnames(transmat))
	}

	if(mismatch>0){
		show("data names do not match")
		return(-1)
	}

	for (i in 1:maxiters){
		newvec = NULL
		if(boolSparceMat) {
			newvec = damping*transmat%*%vec+restart*query
		} else {
			rr = sapply(mylistm,function(x){vec[x$colidxs]%*%(x$vals/mycolsums[x$colidxs])})
			newvec = damping*rr+restart*query
		}
		diff = sum(abs(newvec-vec))
		vec = newvec
#		show(c(i,signif(diff,3)))
		if(diff < thresh) { break }
	}
	vec = as.vector(vec)
	names(vec) = names(query)
	return(list("iter" = i, "diff" = diff, "vec"=vec))
}


#' DRaWR
#'
#' This function runs the DRaWR two stage random walk with restart method.
#' @param possetfile (string): location of file containing location of gene sets to test.
#' @param unifile (string): location of file listing gene universe.
#' @param networkfile (string): location of file containing network contents.
#' @param outdir (string): prefix of location of file to write performance results (optionally prediction results).
#' @param restarts (vector): vector of restart values to test. Default is c(0.7).
#' @param nfolds (int): number of folds for cross validation, 1 is no cross-validation. Default is 4.
#' @param st2keep (int): number of property nodes to keep in second stage for each property type. Default is 50.
#' @param undirected (bool): boolean to make network undirected.
#' @param unweighted (bool): boolean to make network unweighted.
#' @param normalize (string): "type" or "none". Default is 'type'.
#' @param maxiters (int): maximum number of allowable iterations. Default is 50.
#' @param thresh (float): threshold for L1 norm convergence. Default is 0.001.
#' @param property_types (vector): list of possible property types. Default is c("go_curated_evidence", "go_inferred_evidence", "pfam_domain").
#' @param writepreds (boolean): write predictions out to a file. Default is FALSE
#' @keywords random walk with restart
#' @import Matrix
#' @importFrom Matrix colSums
#' @importFrom ROCR prediction
#' @importFrom ROCR performance
#' @export
#' @examples
#' DRaWR(possetfile = system.file("extdata", "sample.setlist", package="DRaWR"),
#'		unifile = system.file("extdata", "sample.uni", package="DRaWR"),
#'		networkfile = system.file("extdata", "sample.edge", package="DRaWR"),
#'		outdir = "output_", restarts = c(.7), nfolds = 1, st2keep = 1,
#'		undirected = TRUE, unweighted = FALSE, normalize = "type", maxiters = 50,
#'		thresh = 0.0001, property_types = c("T1", "T2"), writepreds = 0)
DRaWR<- function(possetfile = "extdata/sample.setlist", unifile = "extdata/sample.uni", networkfile = "extdata/sample.edge", outdir = "output_", restarts = c(.7), nfolds = 1, st2keep = 1, undirected = TRUE, unweighted = FALSE, normalize = "type", maxiters = 50, thresh = 0.0001, property_types = c("allen_brain_atlas", "chip_binding", "gene_ontology", "motif_u5", "pfam_domain", "T1", "T2"), writepreds = 0){

#	library(Matrix)
#	library(ROCR)

	# set up results values
	uni = tail(unlist(strsplit(unifile, "/")),1)
	uni = gsub(".txt","",uni)
	network = tail(unlist(strsplit(networkfile, "/")),1)
	network = gsub(".edge","",network)
	network = gsub(".txt","",network)
	possetname = tail(unlist(strsplit(possetfile, "/")),1)
	possetname = gsub(".txt","",possetname)
	weight = 'weight'
	if(unweighted){ weight = 'unweight' }
	directed = 'dir'
	if(undirected){ directed = 'undir' }

	# set up results output
	restable = NULL
	resfile = paste(sep="", outdir, uni, ".", network, ".", directed, ".", weight, ".", normalize, ".", maxiters, ".", thresh, ".", st2keep, ".", nfolds, ".", paste(collapse="_", restarts), ".", possetname, ".stats")
	if(file.exists(resfile)){
		show(c("WARNING: File already exists ", resfile));
		#return(1)
	}
	write.table(restable, resfile, quote=F, sep="\t", row.names=F)
	show(resfile)

	# read in edge file
	show('Reading Original File')
	edges = read.table(networkfile)
	colnames(edges) = c("src", "target", "weight","type")
	edges$weight <- as.numeric( as.character( edges$weight ) )

	# check for strange values, 0, negs, max
	wmin = min(edges$weight)
	wmins = edges[which(edges$weight==wmin),]
	show("min edges")
	show(wmins[1:min(3,length(wmins[,1])),])
	if(wmin <= 0){
		show("Invalid Edge Weights")
		return(-1)
	}
	wmax = max(edges$weight)
	wmaxs = edges[which(edges$weight==wmax),]
	show("max edges")
	show(wmaxs[1:min(3,length(wmaxs[,1])),])

	# collect info about feature nodes
	all_etypes = as.character(unique(edges$type))
	prop_etypes = intersect(all_etypes, property_types)
	ntypes = length(prop_etypes)
	features = unique(edges[which(edges$type %in% prop_etypes),c("src","type")])
	featnodes = as.character(features[,"src"])
	nfeats = length(featnodes)
	nkeep = min(st2keep*ntypes,nfeats)
	show(c("num_prop_types:", ntypes, 'num_prop_nodes', nfeats))

	# make undirected
	if(undirected){
		show("Making Network Undirected")
		tmpsrc = c(as.character(edges$src), as.character(edges$target))
		tmptarget = c(as.character(edges$target), as.character(edges$src))
		tmpweight = c(edges$weight, edges$weight)
		tmptype = c(as.character(edges$type), as.character(edges$type))
		edges = data.frame(tmpsrc, tmptarget, tmpweight, tmptype)
		colnames(edges) = c("src", "target", "weight","type")
	}

	# make unweighted
	if(unweighted){
		show("Removing Edge Weights")
		edges[,"weight"] = 1
	}

	# resolve node duplication within edge_types
	show("Removing (n1,n2,et) Duplicates")
	edges = aggregate(x = edges$weight, by = list(edges$src,edges$target,edges$type), FUN=max)[,c(1,2,4,3)]
	colnames(edges) = c("src", "target", "weight","type")

	# edge type normalization
	if(normalize == "type"){
		show("Edge Type Normalization")
		typetable = aggregate(weight~type, data = edges, FUN=sum)
 		show(typetable)
		rownames(typetable) = as.character(typetable[,1])
		edges$weight = edges$weight / typetable[as.character(edges$type),"weight"]
	}

	# resolve node duplication across edge_types

	#output normalized matrix
	if(FALSE){
		write.table(edges, paste(networkfile,".clean"), quote=F, sep="\t", col.names=F, row.names=F)
	}

	# calculate what matrix format is need (sparce matrix or personal list format)
	forcelarge = 0
	node_estimate = length(unique(c(as.character(edges[,1]),as.character(edges[,2]))))
	boolSparceMat = (node_estimate^2 < .Machine$integer.max) * (1-forcelarge)   # will fit in sparce mat

	# create adjacency matrix
	transmat = NULL
	nodenames = NULL
	colsum = NULL

	if(boolSparceMat) {
		# convert edgelist to matrix
		n1Matrix = threeCol2MaxMat( as.character(edges[,"src"]),  as.character(edges[,"target"]), as.numeric(edges[,"weight"]))
		nodenames = rownames(n1Matrix)
		# column normalize
		colsum = colSums(n1Matrix)
		transmat = t(t(n1Matrix)/colsum)
		rm(n1Matrix)
	} else {  # must use personal list format
		ll = threeCol2listMat( as.character(edges[,"src"]),  as.character(edges[,"target"]), as.numeric(edges[,"weight"]))
		transmat = ll
		nodenames = ll$avals
		colsum = ll$colsums
		rm(ll)
	}
	nnodes = length(nodenames)

	# Stage 0
	# read gene universe file
	universe = read.table(unifile)
	rownames(universe) = as.character(universe[,1])
	if(dim(universe)[2]<2){
		universe = cbind(as.character(universe[,1]), rep(1, length(universe[,1])) )
	}
	uniIDs = sort(intersect(nodenames, as.character(unique(universe[,1]))))

	if(length(uniIDs)<1){ return(-1)}

	for (restart in restarts){
		# restart = 0.3
		# store results in table
		evaltabstart = cbind(nodenames,-1,0)
		colnames(evaltabstart) = c("node", "type", "universe")
		rownames(evaltabstart) = nodenames
		evaltabstart[uniIDs,"universe"] = 1

		# store node type info in table
		evaltabstart[as.character(features[,1]),"type"] = as.character(features[,2])

		# run baseline rwr network
		basefile = paste(sep="", outdir, uni, ".", network, ".", directed, ".", weight, ".", normalize, ".", maxiters, ".", thresh, ".", restart, ".base")
		show(basefile)

		blankvec = structure(rep(0,nnodes), names = nodenames)
		startvec = blankvec + 1 / nnodes
		biter=0
		if(!file.exists(basefile)){
			query = blankvec
			midxs = match(uniIDs, nodenames)
			query[midxs] = as.numeric(universe[uniIDs,2])

			rwr_res = RWR(boolSparceMat, transmat, restart, query, startvec, maxiters, thresh)

			biter = rwr_res$iter
			evaltabstart = cbind(evaltabstart,as.numeric(rwr_res$vec))
			colnames(evaltabstart)[4] = "baseline"

			if(!file.exists(basefile)){
				write.table(evaltabstart, basefile, quote=F, sep="\t", row.names=T, col.names=NA)
			}
		} else { # read in base results into startvec
			evaltabstart = read.table(basefile)
		}

		startvec[nodenames] = as.numeric(evaltabstart[nodenames,"baseline"])  #starting distribution at baseline may speed convergence

		# read positive set names
		tmpname = tail(unlist(strsplit(possetfile, "/")),1)
		possetpref = gsub(tmpname,"",possetfile)
		possets = as.character(read.table(possetfile)[,1])

		for (qfile in possets){
			#qfile = possets[1]
			queryfile = paste(sep="/", possetpref, qfile)

			posset = tail(unlist(strsplit(queryfile, "/")),1)
			posset = gsub(".txt","",posset)

			show(queryfile)
			if(file.info(queryfile)$size == 0){show(c("Empty file ", queryfile)); next}

			# Stage 1
			# read query file
			query_gs = read.table(queryfile)
			if(dim(query_gs)[2]<2){
				query_gs = cbind(as.character(query_gs[,1]), rep(1, length(query_gs[,1])) )
			}
			rownames(query_gs) = as.character(query_gs[,1])
			queryIDs = sort(intersect(uniIDs, as.character(unique(query_gs[,1]))))
			nquery = length(queryIDs)

			set.seed(041185)
			folds = rep(0, nquery)
			if(nfolds>1){
				folds = sample(cut(seq(1,nquery),breaks=nfolds,labels=FALSE))
			}

			for(iter in 1:nfolds){
				# iter = 1
				## read query set
				outfile = paste(sep="", outdir, uni, ".", network, ".", directed, ".", weight, ".", normalize, ".",  maxiters, ".", thresh, ".", st2keep, ".", nfolds, ".", restart, ".", posset, ".", iter, ".rwr")
				show(outfile)
				if(file.exists(outfile)){
					show(c("WARNING: File already exists ", outfile));
					#next;
				}

				## separate training and testing
				train_idxs = which(folds!=iter)
				train_nidxs = queryIDs[train_idxs]
				ntrain = length(train_nidxs)

				test_idxs = which(folds==iter)
				testuni = as.character(setdiff(uniIDs,train_nidxs))
				if(nfolds==1){
					test_idxs = train_idxs
					testuni = as.character(uniIDs)
				}
				test_nidxs = queryIDs[test_idxs]
				ntest = length(test_nidxs)
				ntestuni = length(testuni)

				if(ntrain*ntest*1.0*ntestuni==0){ ## either no training or testing examples
					row1 = c(network, directed, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "baseline", -1000, -1000, -1000)
					row2 = c(network, directed, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage1", -1000, -1000, -1000)
					row3 = c(network, directed, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "diff", -1000, -1000, -1000)
					row4 = c(network, directed, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage2", -1000, -1000, -1000)
					restable = rbind(restable, row1, row2, row3, row4)
					next
				}

				# run RWR on training set
				query = blankvec
				midxs = match(train_nidxs, nodenames)
				query[midxs] = as.numeric(query_gs[train_nidxs,2])
				rwr_res = RWR(boolSparceMat, transmat, restart, query, startvec, maxiters, thresh)

				# store results
				qiter = rwr_res$iter
				evaltab = cbind(evaltabstart,0,0,as.numeric(rwr_res$vec[nodenames]))
				colnames(evaltab) = c("node", "type", "universe", "baseline", "train", "test", "stage1")
				diff = as.numeric(evaltab[,"stage1"]) - as.numeric(evaltab[,"baseline"])
				evaltab = cbind(evaltab,diff)
				evaltab[train_nidxs,"train"] = 1
				evaltab[test_nidxs,"test"] = 1

				# eval of baseline rwr on pos set
				model = prediction(as.numeric(evaltab[testuni,"baseline"]), evaltab[testuni,"test"])
				auc  = performance(model, "auc")
				perf = performance(model,"tpr","fpr")
				aucval = round(as.numeric(slot(auc, "y.values")),3)
				row	= c(network, directed, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "baseline", aucval, 0, length(uniIDs))
				show(paste(collapse=' ', row[12:15]))
				restable = rbind(restable, row)

				# eval of query rwr on pos set
				model = prediction(as.numeric(evaltab[testuni,"stage1"]), evaltab[testuni,"test"])
				auc  = performance(model, "auc")
				perf = performance(model,"tpr","fpr")
				aucval = round(as.numeric(slot(auc, "y.values")),3)
				row = c(network, directed, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage1", aucval, qiter, ntrain)
				show(paste(collapse=' ', row[12:15]))
				restable = rbind(restable, row)
				#plot(perf, lwd = 5)

				# eval of diff rwr on pos set
				model = prediction(as.numeric(evaltab[testuni,"diff"]), evaltab[testuni,"test"])
				auc  = performance(model, "auc")
				perf = performance(model,"tpr","fpr")
				aucval = round(as.numeric(slot(auc, "y.values")),3)
				row = c(network, directed, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "diff", aucval, 0, ntrain)
				show(paste(collapse=' ', row[12:15]))
				restable = rbind(restable, row)
				#plot(perf, lwd = 5)

				colnames(restable) = c("network", "direct", "weight", "normalize", "uni", "restart", "maxiters", "thresh", "st2keep", "posset", "nfolds", "iter", "stage", "aucval", "rwr_iters", "ntrain")
				write.table(restable, resfile, quote=F, sep="\t", row.names=F)

				# stage 2
				# can skip in no feature nodes
				if(nfeats==0){
					row4 = c(network, directed, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage2", -1000, -1000, -1000)
					restable = rbind(restable, row4)
					write.table(restable, resfile, quote=F, sep="\t", col.names=F, row.names=F)
					if(iter == 1 && writepreds){
						write.table(evaltab, outfile, quote=F, sep="\t", row.names=F, col.names=T)
					}
					next;
				}

				# extract best features nodes
				keep = rep(-1, nnodes)
				evaltab = cbind(evaltab, keep)
				ss = sort(as.numeric(evaltab[featnodes,"diff"]), decreasing=T, index.return = T)
				sortedfeats = featnodes[ss$ix]
				keepfeats = sortedfeats[1:nkeep]
				evaltab[featnodes,"keep"] = 0
				evaltab[keepfeats,"keep"] = 1

				# keep all non-property edges
				newedges = edges[which(edges$type %in% setdiff(all_etypes, prop_etypes)),]

				# add in edges connected to kept features
				keepidxs = which(edges[,"src"] %in% keepfeats | edges[,"target"] %in% keepfeats)
				newedges = rbind(newedges, edges[keepidxs,])

				tmpnames = unique(c(as.character(newedges[,1]),as.character(newedges[,2])))
				train_nidxs2 = intersect(train_nidxs,tmpnames)
				test_nidxs2 = intersect(test_nidxs,tmpnames)
				testuni2 = intersect(testuni,tmpnames)

				# can skip if no overlap with the train or test set
				if(length(train_nidxs2)*length(test_nidxs2)*1.0*length(testuni2)==0){
					row4 = c(network, directed, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage2", -1000, -1000, -1000)
					restable = rbind(restable, row4)
					write.table(restable, resfile, quote=F, sep="\t", col.names=F, row.names=F)
					if(iter == 1 && writepreds){
						write.table(evaltab, outfile, quote=F, sep="\t", row.names=F, col.names=T)
					}
					next;
				}

				# edge type renormalized
				typetable2 = NULL
				if(normalize=="type"){
					typetable2 = aggregate(weight~type, data = newedges, FUN=sum)
					#show(typetable2)
					rownames(typetable2) = as.character(typetable2[,1]	)
					typetable2[as.character(newedges$type),"weight"]
					newedges$weight = newedges$weight / typetable2[as.character(newedges$type),"weight"]
					#show(newedges[1:5,])
				}

				transmat2 = NULL
				nodenames2 = NULL
				colsum2 = NULL
				# check if sparce matrix will work
				boolSparceMat2 = (length(tmpnames)^2 < .Machine$integer.max) * (1-forcelarge)  # will fit in sparce mat
				if(boolSparceMat2) {
					# convert edgelist to matrix
					n2Matrix = threeCol2MaxMat( as.character(newedges[,"src"]),  as.character(newedges[,"target"]), as.numeric(newedges[,"weight"]))
					nodenames2 = rownames(n2Matrix)
					# column normalize
					colsum2 = colSums(n2Matrix)
					transmat2 = t(t(n2Matrix)/colsum2)
					rm(n2Matrix)
				} else {  # must use personal list format
					ll2 = threeCol2listMat( as.character(newedges[,"src"]),  as.character(newedges[,"target"]), as.numeric(newedges[,"weight"]))
					transmat2 = ll2
					nodenames2 = ll2$avals
					colsum2 = ll2$colsums
					rm(ll2)
				}
				nnodes2 = length(nodenames2)
				rm(newedges)

				# run RWR on s2 network on training set
				blankvec2 = structure(rep(0,nnodes2), names = nodenames2)
				query2 = blankvec2
				midxs = match(train_nidxs2, nodenames2)
				query2[midxs] = as.numeric(query_gs[train_nidxs2,2])
				startvec2 = blankvec2 + 1 / nnodes2
				rwr_res = RWR(boolSparceMat2, transmat2, restart, query2, startvec2, maxiters, thresh)

				# store results
				q2iter = rwr_res$iter
				stage2 = rep(0,nnodes)
				evaltab = cbind(evaltab, stage2)
				evaltab[nodenames2,"stage2"] =  rwr_res$vec[nodenames2]

				# eval of stage2 rwr on pos set
				model = prediction(as.numeric(evaltab[testuni2,"stage2"]), evaltab[testuni2,"test"])
				auc  = performance(model, "auc")
				perf = performance(model,"tpr","fpr")
				aucval = round(as.numeric(slot(auc, "y.values")),3)
				row	= c(network, directed, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage2", aucval, q2iter, length(train_nidxs2))
				show(paste(collapse=' ', row[12:15]))
				restable = rbind(restable, row)
				write.table(restable, resfile, quote=F, sep="\t", row.names=F)

				if(iter == 1 && writepreds){
					write.table(evaltab, outfile, quote=F, sep="\t", row.names=F, col.names=T)
				}

			} #end iter
		} #end queryset
	} #end restart
} #end function


#' Sample transition matrix.
#'
#' A Matrix containing tthe normalized transition matrix from the test network
#'
#' @format a Matrix containing the normalized transition matrix from the test network.
"transmat"