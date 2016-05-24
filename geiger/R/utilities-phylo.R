
## GENERIC FUNCTION
heights <- function(x)
UseMethod("heights")


heights.phylo=function(x){
    phy=x
	phy <- reorder(phy, "postorder")
	n <- length(phy$tip.label)
	n.node <- phy$Nnode
	xx <- numeric(n + n.node)
	for (i in nrow(phy$edge):1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	root = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
	labs = c(phy$tip.label, phy$node.label)
	depth = max(xx)
	tt = depth - xx
	idx = 1:length(tt)
	dd = phy$edge.length[idx]
	mm = match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
	dd = c(phy$edge.length, root)[mm]
	ss = tt + dd
	res = cbind(ss, tt)
	rownames(res) = idx
	colnames(res) = c("start", "end")
	res = data.frame(res)
	res
}


heights.multiPhylo=function(x){
    phy=x
	phy=hashes.phylo(phy)
	hh=unique(unlist(lapply(phy, function(x) x$hash)))
	times=lapply(phy, function(x) {tmp=heights.phylo(x); tmp$hash=x$hash; tmp})
	if(is.null(names(times))) names(times)=1:length(times)
	out=lapply(hh, function(hash){
        tmp=sapply(times, function(dat){
            if(hash%in%dat$hash) dat[which(dat$hash==hash),c("start","end")]  else return(c(0,0))
        })
        res=data.frame(t(tmp))
        rownames(res)=names(times)
        names(res)=c("start","end")
        gg=apply(res, 1, function(x) all(x==0))
        if(any(gg)) res=res[-which(gg),]
        return(res)
    })
	attr(out,"hash")=hh
	out
}

## as.phylo.hyplo
as.phylo.hphylo=function(x, ...){
     x$desc=NULL
     x$hash=NULL
     cl=class(x)
     class(x)=cl[!cl%in%"hphylo"]
     x
}

.unique.phylo=function(phy){
    # phy is assumed to have tips that are redundant (and whose exemplars are monophyletic)
    # prunes tree to leave one member of each unique label

	mon=table(phy$tip.label)
	if(all(mon==1)) return(phy)
	mon=mon[mon>1]
	for(i in 1:length(mon)){
		m=names(mon[i])
		drop=which(phy$tip.label==m)
		drop=drop[2:length(drop)]
		phy$tip.label[drop]="null"
	}
	phy=.drop.tip(phy, phy$tip.label[phy$tip.label=="null"])
	return(phy)
}


unique.phylo=function(x, incomparables=FALSE, ...){
    # phy: has tip.labels that are not unique
    # returns phylogeny with unique tip labels
    # if a taxon (indicated by multiple tip labels of same string) is non-monophyletic, all tips of that taxon are removed with warning
    # likely used where an exemplar phylogeny is needed and phy$tip.label is 'tricked' into a non-unique vector of names

	phy=x
	if(incomparables) warning("'incomparables' exerts no effect in this setting")
	.exemplar.monophyly=function(tip, phy) {
        # tip: occurs multiply in phy$tip.label
		nn=which(phy$tip.label==tip)
		if(length(nn)==1) return(TRUE)
		anc=getMRCA(phy, nn)
		dd=unique(phy$tip.label[.get.descendants.of.node(anc, phy, tips=TRUE)])
		if(length(dd)==1) return(TRUE) else return(FALSE)
	}


	tt=table(phy$tip.label)
	if(!any(tt>1)) {
		return(phy)
	} else {
		todo=names(tt[tt>1])
        #		cat("checking monophyly of groups...\n")
		mon=sapply(todo, function(t) {
            #				   cat(paste("\n\t",t, sep=""))
            .exemplar.monophyly(t, phy)
        })
        #		cat("\n")

		if(any(!mon)) {
			warning(paste("non-monophyletic lineages encountered:\n\t", paste(todo[!mon], collapse="\n\t"), "\n", sep=""))
			phy=.drop.tip(phy, names(!mon))
		}
		return(.unique.phylo(phy))
	}
}

unique.multiPhylo=function(x, incomparables=FALSE, ...){
	phy=x
	if(incomparables) warning("'incomparables' exerts no effect in this setting")
	ss=sapply(phy, digest)
	if(any(dd<-duplicated(ss))){
		sub=phy[-which(dd)]
		class(sub)="multiPhylo"
		return(sub)
	} else {
		return(phy)
	}
}

# wrapper for plot.phylo: plot tree with internal nodes labelled with ape numbering
plotNN <- function (phy, time = TRUE, margin = TRUE, ...) {
    phy$node.label <- (length(phy$tip.label) + 1):max(phy$edge);
    plot.phylo(phy, show.node.label = TRUE, no.margin = !margin, cex = 0.5, ...);
    if (time && margin) {
        axisPhylo(cex.axis = 0.75);
    }
}

# Get b and d values from r (b-d) and epsilson (d/b)
# Used in previous version of program; now in terms of r and epsilon
# Possibly of use to users wishing to translate results
get.bd <- function (r, epsilon) {
      b <- r/(1 - epsilon);
      d <- b - r;   # Alternatively: d <- epsilon * r/(1 - epsilon)
      return(list(b = b, d = d));
}


#general phylogenetic utility for determining whether a node is the root of the phylogeny
#author: JM EASTMAN 2011

is.root <- function (node,phy) {
	if (node == (Ntip(phy) + 1)) return(TRUE) else return(FALSE);
}


.nodefind.phylo=function(phy, label){
	N=Ntip(phy)
	ww=match(label, phy$node.label)
	if(!all(is.na(ww))) {
		return(N+ww)
	} else {
		warning(paste("Encountered no node.label for ",label,sep=""))
		return(NULL)

	}
}

.cache.descendants=function(phy){
    # fetches all tips subtended by each internal node

	N=as.integer(Ntip(phy))
	n=as.integer(Nnode(phy))

	phy=reorder(phy, "postorder")

	zz=list( N=N,
    MAXNODE=N+n,
    ANC=as.integer(phy$edge[,1]),
    DES=as.integer(phy$edge[,2])
    )

	res=.Call("cache_descendants", phy=zz, package="geiger")
	return(res)
}


.polytomy.phylo=function(tips, age=1){
	N=length(tips)+1
	edge=cbind(N,1:(N-1))
	length=rep(age, N-1)
	phy=list(edge=edge, edge.length=length, tip.label=tips, Nnode=1)
	class(phy)="phylo"
	return(phy)
}


nodelabel.phylo=function(phy, taxonomy, strict=TRUE, ncores=NULL){
    # all phy$tip.label must be in taxonomy
    # taxonomy: exclusivity highest on left, lowest on right (species, genus, family, etc., as columns)
    # columns in 'taxonomy' should ONLY be taxonomic ranks


    #	tt=as.matrix(taxonomy)
    #	tt[!is.na(tt)]=""
    #	drp=apply(tt, 1, function(x) all(x==""))
    #	if(any(drp)) taxonomy=taxonomy[-which(drp),]


	taxonomy=cbind(rownames=rownames(taxonomy),taxonomy)
	rank="rownames"

	taxonomy=as.data.frame(as.matrix(taxonomy),stringsAsFactors=FALSE)

	if(!all(xx<-phy$tip.label%in%taxonomy[,rank])) {
		warning(paste("taxa not found in 'taxonomy':\n\t", paste(phy$tip.label[!xx], collapse="\n\t"), sep=""))
	}
	taxonomy[taxonomy==""]=NA

    op<-orig<-options()
    op$expressions=max(op$expressions, 500000)
	options(op)

	unmatched=phy$tip.label[!xx]
	idx=match(phy$tip.label, taxonomy[,rank])
	tax=taxonomy[idx,]

	labels=unique(unlist(tax[,-which(names(tax)==rank)]))
	labels=labels[!is.na(labels)]

	dat=tax[,-which(names(tax)==rank)]
	hashes_labels=character(length(labels))
	zz=tax[,rank]
	tips=phy$tip.label[xx]
	for(i in 1:ncol(dat)){
		uu=unique(dat[,i])
		uu=uu[!is.na(uu)]
		for(j in uu){
			cur=zz[which(dat[,i]==j)]
			if(length(cur)>1){
				hashes_labels[which(labels==j)]=.hash.tip(cur, tips)
			} else {
				hashes_labels[which(labels==j)]=NA
			}

		}
	}
	names(hashes_labels)=labels
	hashes_labels=hashes_labels[!is.na(hashes_labels)]

	# redundancies for labels
	tmp=table(hashes_labels)
	if(any(tmp[tmp>1]->redund)){
		root=.hash.tip(tips, tips)
		for(r in names(redund)){
			if(r==root){
				rdx=which(hashes_labels==r)
				warning(paste("redundant labels encountered at root:\n\t", paste(names(hashes_labels[rdx]), collapse="\n\t"), sep=""))
				hashes_labels=hashes_labels[-rdx[!rdx==min(rdx)]]
			} else {
				rdx=which(hashes_labels==r)
				warning(paste("redundant labels encountered:\n\t", paste(names(hashes_labels[rdx]), collapse="\n\t"), sep=""))
				hashes_labels=hashes_labels[-rdx[!rdx==max(rdx)]]
			}
		}
	}

    #	cat("resolving descendants for splits in tree...\n")
	tmp=hashes.phylo(phy, tips, ncores=ncores)
	hashes_tree=tmp$hash
	phy$node.label=rep("",max(phy$edge))
	mm=match(hashes_tree, hashes_labels)
	nodelabels=ifelse(is.na(mm), "", names(hashes_labels[mm]))
	nodelabels[is.na(nodelabels)]=""
	nodelabels=nodelabels[(Ntip(phy)+1):max(phy$edge)]
	tmp=table(nodelabels[nodelabels!=""])
	if(any(tmp[tmp>1]->redund)){
		for(r in names(redund)){
			rdx=which(nodelabels==r)
			nodelabels[rdx[rdx!=min(rdx)]]=""
		}
	}

	phy$node.label=nodelabels

	edges=NULL

	desc=.cache.descendants(phy)$tips[-c(1:Ntip(phy))]
    tidx=match(tips, phy$tip.label)

	FUN=function(taxon){
		nm=rownames(dat)
		dat=as.matrix(dat, ncol=ncol(dat))
		rownames(dat)=nm
		if(!taxon%in%dat) {
			try=agrep(taxon, unique(c(dat)), value=TRUE)
			if(length(try)){
				warning(paste(sQuote(taxon), " not encountered in 'taxonomy'\n\nIntended search may have been:\n\t", paste(try, collapse="\n\t", sep=""), sep=""))
			} else {
				warning(paste(sQuote(taxon), " not encountered in 'taxonomy'", sep=""))

			}

			return(NULL)
		}
		expected=rownames(which(dat==taxon, arr.ind=TRUE))
		if(length(expected)==1){ # single tip
			x=which(phy$tip.label==expected)
			hs=.hash.tip(expected, tips)
			res=list(unexpected=c(), missing=c())
			attr(res, "node")=x
			attr(res, "hash")=hs
			attr(res, "expected")=sort(expected)
			return(res)

		}

		bin=as.integer(tips%in%expected)

        #		rownames(edges)=1:nrow(edges)
		N=Ntip(phy)

        ff=.get.parallel(ncores)

		dist=unlist(ff(desc, function(x) {
            y=tidx%in%x
            sum(abs(y-bin))
        }))
		nearest=N+which(dist==min(dist))
		hs=character(length(nearest))
		for(i in 1:length(nearest))hs[i]=.hash.tip(edges[nearest[i],], tips)
		res=lapply(nearest, function(x) {
            tt=tips[which(tidx%in%desc[[x-N]])]
            unexpected=sort(setdiff(tt, expected)) ## in tree but unexpected
            missing=sort(setdiff(expected, tt)) ## expected but not in tree
            tmp=list(unexpected=unexpected, missing=missing)
            hs=.hash.tip(edges[x,], tips)
            attr(tmp, "node")=x
            attr(tmp, "hash")=hs
            attr(tmp, "expected")=sort(expected)
            return(tmp)
        })
		res
	}

	# missed labels
	mm=match(hashes_labels, hashes_tree)
	if(any(is.na(mm))){
		mss=hashes_labels[is.na(mm)]
		if(!strict){
			N=Ntip(phy)

			ll=list()
			nm=rev(names(mss))
			for(x in 1:length(nm)){
				ll[[x]]=FUN(nm[x])
			}
			near_tmp=ll
			names(near_tmp)=nm
			near=lapply(1:length(near_tmp), function(idx){
                tmp=near_tmp[[idx]]
                x=nm[idx]
                if(is.null(tmp)) return(c())
                if(length(tmp)==1) {
                    nd=attributes(tmp[[1]])$node
                    names(nd)=paste("\"", x, "\"", sep="")
                    if(phy$node.label[nd-N]=="") return(nd) else return(c())
                } else {
                    return(c())
                }
			})
			near=unlist(near)
			dd=duplicated(near)
			true_missed=nm[!nm%in%gsub("\"","",names(near))]
			if(any(dd)) near=near[-which(dd)]
			phy$node.label[near-N]=names(near)
		}
	}

	nn=gsub("\"", "", unique(phy$node.label[phy$node.label!=""]))
	hl=names(hashes_labels)
	ms=sort(hl[!hl%in%nn])
	phy$missed=ms
	if(length(ms)){
		warning(paste("labels missing from 'phy':\n\t", paste(ms, collapse="\n\t"), sep=""))
	}

	phy$FUN=FUN

	phy
}




#general phylogenetic utility for returning first ancestor (as a numeric referent) of the supplied node
#author: JM EASTMAN 2010

.get.ancestor.of.node <-
function(node, phy) {
	return(phy$edge[which(phy$edge[,2]==node),1])
}




#general phylogenetic utility for returning all ancestors (listed as given in phy$edge[,2]) of a node
#author: JM EASTMAN 2010

.get.ancestors.of.node <-
function(node, phy) {
	a=c()
	if(node==(Ntip(phy)+1->root)) return(NULL)
	f=.get.ancestor.of.node(node, phy)
	a=c(a,f)
	if(f>root) a=c(a, .get.ancestors.of.node(f, phy))
	return(a)
}


#general phylogenetic utility for returning the first (usually, unless a polytomy exists) two descendants of the supplied node
#author: JM EASTMAN 2010

.get.desc.of.node <-
function(node, phy) {
	return(phy$edge[which(phy$edge[,1]==node),2])
}



#general phylogenetic utility for returning all descendants (listed as given in phy$edge[,2]) of a node (excluding the supplied node)
#author: JM EASTMAN 2011

.get.descendants.of.node <-
function(node, phy, tips=FALSE){
	n=Ntip(phy)
	all=ifelse(tips, FALSE, TRUE)
	out <- .Call("get_descendants", tree=list(
    NODE = as.integer(node),
    ROOT = as.integer(n+1),
    ALL = as.integer(all),
    ENDOFCLADE = as.integer(dim(phy$edge)[1]),
    ANC = as.integer(phy$edge[,1]),
    DES = as.integer(phy$edge[,2])),
    PACKAGE = "geiger")
	res=out$TIPS
	if(!length(res)) res=NULL
	return(res)
}

.mrca=function(labels, phy){
	mm=labels

	if(all(is.character(labels))){

		ll=c(phy$tip.label, phy$node.label)
		mm=match(labels, ll)
		if(any(is.na(mm))) stop("Some 'labels' not encountered in 'phy'")


	}

	if(!all(is.numeric(mm))) stop("Supply 'labels' as a character or integer vector")

    if(length(u<-unique(mm))==1) return(u)

	aa=unlist(lapply(mm, function(x) .get.ancestors.of.node(x, phy)))
	tt=table(aa)
	max(as.integer(names(tt[tt==length(labels)])))
}

is.phylo=function(x) "phylo"%in%class(x)

## FUNCTIONS
## grabs most exclusive tips from clade definitions (and whose tips can then be reconciled with a taxonomic database -- a lookup table)
# allows for recursion in clade definitions (clade defined in part using another clade definition)
# returns trees representing each clade definition
# 'nested': an important variable -- this is the subset of clades that are defined (recursively) within 'clades'
phylo.clades=function(clades, phy=NULL, unplaced=TRUE, ncores=NULL){

    ## give 'phy' as a multiPhylo, named list of trees (whose labels appear in the clade defs)
    ## clades:
    #	clades=list(
    #			 Sirenoidea=c("Siren", "Pseudobranchus"),
    #			 Ambystomatidae=c("Ambystomatidae", "Plethodontidae"),
    #			 Cryptobranchoidea=c("Hynobiidae", "Cryptobranchus_alleganiensis"),
    #			 CAUDATA=c("Sirenoidea","Salamandroidea","Cryptobranchoidea")
    #	)

    if (!is.null(phy)) {
        if (!"multiPhylo" %in% class(phy) | is.null(phynm <- names(phy)))
        stop("Supply 'phy' as a 'multiPhylo' object with names")
        tmp = character(length(phy))
        for (i in 1:length(phy)) {
            cur = phy[[i]]
            cur$edge.length = NULL
            tre = write.tree(cur)
            tmp[i] = gsub(";", "", tre)
        }
        phy = tmp
        names(phy) = phynm
        for (i in 1:length(clades)) {
            cur = clades[[i]]
            if (any(cur == names(clades)[i]))
            stop("Encountered self-referential clade")
            mm = match(cur, phynm)
            if (any(!is.na(mm))) {
                ww = which(!is.na(mm))
                for (j in 1:length(ww)) {
                    clades[[i]][ww[j]] = phy[[mm[ww[j]]]]
                }
            }
        }
    }

    ll = sapply(clades, length)
    if (any(ll == 1)) {
    	unis=clades[ll==1]
    	clades=clades[ll!=1]
        warning("Non-splitting lineages found in 'clades'")
    } else {
        unis=NULL
    }
    cc = unique(c(unlist(clades)))
    tt = cc %in% names(clades)
    nested = cc[tt]
    is.nestedclade = function(clade, nested) {
        if (clade %in% nested)
        return(TRUE)
        else return(FALSE)
    }

    fetch_nestedclades = function(clade, clades, nested) {
        if (is.nestedclade(clade, nested)) {
            desc = clades[[clade]]
            new = as.list(desc)
            names(new) = desc
            tt = sapply(new, is.nestedclade, nested)
            for (i in 1:length(new)) {
                if (tt[i])
                new[[i]] = fetch_nestedclades(new[[i]], clades, nested)
                else new[[i]] = NA
            }
        }
        else {
            new = NA
        }
        return(new)
    }

    paths_through_clades = function(clades, nested) {
        nn = names(clades)
        res = lapply(nn, function(x) {
            dd = clades[[x]]
            y = lapply(dd, fetch_nestedclades, clades, nested)
            names(y)[1:length(dd)] = dd
            y
        })
        names(res) = nn
        res
    }

    unplaced_phy = function(phy, cladepath) {
        if (any(names(cladepath) == "unplaced")) {
            tips = names(cladepath$unplaced)
            y = .polytomy.phylo(tips)
            y$edge.length = NULL
            new = bind.tree(phy, y)
            new = .drop.tip(new, "unplaced")
            return(new)
        }
        else {
            return(phy)
        }
    }

    tree_cladepath = function(cladepath, nested) {
        print_group = function(cladepath) {
            xx = sapply(names(cladepath), is.nestedclade, nested)
            middle = character(length(xx))
            for (i in 1:length(xx)) {
                if (xx[i]) {
                    new = cladepath[[names(xx)[i]]]
                    middle[i] = print_group(new)
                }
                else {
                    middle[i] = names(cladepath)[i]
                }
            }
            paste("(", paste(middle, collapse = ", "), ")", sep = "")
        }
        tmp = paste(print_group(cladepath), ";", sep = "")
        phy = read.tree(text = tmp)
        return(unplaced_phy(phy, cladepath))
    }

    cladepaths = paths_through_clades(clades, nested)

    if (unplaced) {
        for (i in 1:length(cladepaths)) {
            cur = cladepaths[[i]]
            if (!is.null(unplc <- attributes(clades[[i]])$unplaced)) {
                unplaced_taxa = lapply(unplc, function(x) return(NA))
                names(unplaced_taxa) = unplc
                cladepaths[[i]]$unplaced = unplaced_taxa
            }
        }
    }

    alt_tip_label=function(phy, unis){
    	phy$alt.tip.label=character(Ntip(phy))
    	pt=phy$tip.label
    	for(i in 1:length(unis)){
    		if(length(cur<-unis[[i]])!=1) stop("'unis' should have a single member in each element")
    		if(!is.na(mm<-match(cur, pt))) phy$alt.tip.label[mm]=names(unis)[i]
    	}
    	phy
    }

    phy = lapply(cladepaths, tree_cladepath, nested)
    if(length(unis)) phy = lapply(phy, alt_tip_label, unis)
    nn=sapply(phy, Ntip)
    midx<-which(nn==max(nn))
   	master=phy[midx]

    ## RESOLVE NODE LABELS for 'master' tree (most comprehensive)
   	if(length(master)==1){
   		master=master[[1]]
   		tt=master$tip.label
   		null=.hash.tip(c(), tt)

   		mm=hashes.phylo(master, tt, ncores=ncores)
   		ss=sapply(phy, function(x) .hash.tip(x$tip.label, tt))
   		if(any(ss==null)){
			warning(paste("The following not encountered:\n\t", paste(names(ss)[which(ss==null)], collapse="\n\t")))
		}
   		ts=table(ss)
   		if(any(ts>1)){
   			drp=numeric()
   			ts=ts[ts>1]
   			for(i in 1:length(ts)){
   				ww=which(ss==names(ts[i]))
   				drp=c(drp, ww[which(ww!=min(ww))])
   			}
   			warning("The following node labels appear redundant:\n\t", paste(names(ss)[drp], collapse="\n\t"))
   			ss=ss[-drp]
   		}
   		xx=match(ss, mm$hash)
   		if(any(is.na(xx))){
   			warning(paste("The following not encountered:\n\t", paste(names(ss)[is.na(xx)], collapse="\n\t")))
   			ss=ss[-which(is.na(xx))]
   		}
   		N=Ntip(master)
   		nd=character(Nnode(master))
   		nd[xx-N]=names(ss)
   		master$node.label=nd
   		phy[[midx]]=master
   	}
    return(phy)
}


lookup.phylo=function(phy, taxonomy=NULL, clades=NULL, ncores=NULL){
    ## taxonomy expected to have first column at same level as tip labels in phy
    ## first row in taxonomy is most exclusive
    ## clade_defs are phylogenetic trees of a clade representation

	if(!is.null(taxonomy)){
		if(!any(taxonomy[,1]%in%phy$tip.label) & !is.null(rownames(taxonomy))) {
			taxonomy=as.data.frame(as.matrix(cbind(species=rownames(taxonomy), taxonomy)),stringsAsFactors=FALSE)
		}
	}

	related_tips=function(tips, phy){
		tips=tips[tips%in%phy$tip.label]
		if(length(tips)<2) stop("related_tips(): 'tips' to be found in 'phy' are too few")
		nd=unique(sapply(tips, function(x) match(x, phy$tip.label)))
		if(length(nd)==1) return(nd)
		anc=getMRCA(phy, nd)
		dd=.get.descendants.of.node(anc, phy, tips=TRUE)
		return(phy$tip.label[dd])
	}

	tips_in_group=function(phy, taxonomy=NULL, clade_def){

		lengths=sapply(clade_def$tip.label, function(grp){
            if(!is.null(taxonomy)){
                ww=which(taxonomy==grp, arr.ind=TRUE)[,"row"]
                if(length(ww)>0){
                    return(TRUE)
                } else {
                    return(FALSE)
                }
            } else {
                return(length(phy$tip.label[which(phy$tip.label%in%grp)])>0)
            }
		})

		if(sum(lengths)<2) return(NA)

		## FIXME: use 'lengths' to determine if 'clade_def' is satisfied (at least two members needed)


		tmp=sapply(clade_def$tip.label, function(grp){
            if(!is.null(taxonomy)){
                ww=which(taxonomy==grp, arr.ind=TRUE)[,"row"]
                if(length(ww)>0){
                    return(unique(taxonomy[ww,1]))
                } else {
                    return(c())
                }
            } else {
                return(phy$tip.label[which(phy$tip.label%in%grp)])
            }
		})

		## most exclusive 'spp' (recognized from tree)
		spanning_taxa=unique(unlist(tmp))
		recovered_taxa=unique(related_tips(spanning_taxa, phy))

		return(recovered_taxa)
	}

	orig_phy=phy

	if(!is.null(taxonomy)){
		## check ordering of taxonomy
		oo=order(apply(taxonomy, 2, function(x) length(unique(x))),decreasing=TRUE)
		if(!all(oo==c(1:ncol(taxonomy)))){
			warning("Assuming 'taxonomy' is not from most to least exclusive")
			taxonomy=taxonomy[,ncol(taxonomy):1] #reverse ordering of columns
		}
		original_taxonomy=taxonomy

		## check rank of tree
		phy$node.label=NULL
		gtips=sapply(phy$tip.label, function(x) unlist(strsplit(gsub(" ", "_", x), "_"))[1])
		check=sapply(list(phy$tip.label, gtips), function(x) sum(x%in%taxonomy[,1]))
		if(max(check)==check[2]) {
			genuslevel_tree=TRUE
			phy$tip.label=unname(gtips)
		} else {
			genuslevel_tree=FALSE
		}
		tips=phy$tip.label

		## prune taxonomy down to members in the phylogeny
		taxonomy=taxonomy[taxonomy[,1]%in%tips,]
		matching_tips=taxonomy[,1]

		## BUILD taxonomic mapping to each row in 'phy'
		mm=match(tips, matching_tips)
		tt=taxonomy[mm,]
		colnames(tt)=names(taxonomy)
	} else {
		tt=c()
	}

	if(!is.null(clades)){
		tips=phy$tip.label
		clade_defs=phylo.clades(clades, ncores=ncores)
        #		cat("resolving clades...\n\t")
		res=lapply(1:length(clade_defs), function(idx) {
            def=clade_defs[[idx]]
            #				   cat(paste(names(clade_defs)[idx], "\n\t", sep=""))
            tt=try(tips_in_group(phy, taxonomy, def),silent=TRUE)
            if(inherits(tt, "try-error")) return(NA) else return(tt)
        })
        #		cat("\n")
		names(res)=names(clade_defs)

		zz=sapply(res, function(x) all(is.na(x)))
		if(any(zz)){
			warning(paste("taxa not represented in 'tips':\n\t", paste(names(res)[zz], collapse="\n\t"), sep=""))
		}
		res=res[!zz]
		clade_defs=clade_defs[!zz]

		## BUILD clade-level mapping to each row in 'phy'
		gg=matrix("", nrow=Ntip(phy), ncol=length(res))
		for(i in 1:length(res)){
			nm=names(clade_defs)[i]
			cur_tips=res[[i]]
			zz=tips%in%cur_tips
			gg[zz,i]=nm
		}
		colnames(gg)=names(clade_defs)
		tt=as.matrix(cbind(tt, gg))
	}

	if(is.null(tt)){
		if(!is.null(nn<-phy$node.label)){
			names(nn)=Ntip(phy)+1:length(nn)
			nn=nn[nn!=""]
			if(length(nn)){
				tt=matrix("", nrow=Ntip(phy), ncol=length(nn))
				dd=.cache.descendants(phy)$tips
				for(i in 1:length(nn)){
					tt[phy$tip.label%in%phy$tip.label[dd[[as.integer(names(nn[i]))]]],i]=nn[i]
				}
				colnames(tt)=nn

			} else {
				return(NULL)
			}
		} else {
			return(NULL)
		}
	}
	phy_mapping=tt
	rownames(phy_mapping)=orig_phy$tip.label
	return(as.data.frame(as.matrix(phy_mapping), stringsAsFactors=FALSE))
}



## FUNCTIONS ##

.fill.taxonomy=function(taxonomy){
	push.taxon=function(taxonomy, indices, column){
		for(n in which(indices)){
			taxonomy[n,column]=paste(taxonomy[n,min((column:ncol(taxonomy))[!is.na(taxonomy[n,column:ncol(taxonomy)])])],paste("R",column,sep=""),sep=".")
		}
		return(taxonomy)
	}

	for(c in 1:(ncol(taxonomy)-1)){
		ii=is.na(taxonomy[,c])
		if(any(ii)) taxonomy=push.taxon(taxonomy, ii, c)
	}
	return(taxonomy)
}

phylo.lookup=function(taxonomy, ncores=NULL) {
    # GENERAL FUNCTION: convert taxonomic 'lookup' table to phylogeny
    # lookup is data.frame with ranks as columns
    # rowlabels of 'taxonomy' are assumed to be tips of the phylogeny
    # rank corresponds to (numeric) position in rank names (R to L) to use in building taxonomic tree; if null, all ranks are used
    # NOTE: taxonomic groups in 'lookup' MUST be ordered by exclusivity from R to L -- e.g., genera must precede families must precede orders
	labels=rownames(taxonomy)
    mm=nrow(taxonomy)
    tax=as.matrix(taxonomy, ncol=ncol(taxonomy))
	rownames(tax)=labels
	occurrences=table(tax)
	occurrences=occurrences[names(occurrences)!="" & occurrences>1]
    if(any(labels%in%names(occurrences))) stop(paste("The following appear to be both ranks and row.names of 'taxonomy':\n\t", paste(labels[labels%in%names(occurrences)], collapse="\n\t"), sep=""))

    if(!any(occurrences==mm)){
        tax=cbind(as.matrix(taxonomy, ncol=ncol(taxonomy)),root="root")
        rownames(tax)=labels
        occurrences=c(occurrences, root=mm)
    }

    f=.get.parallel(ncores)

    #clds=f(names(occurrences), function(x) {
    #		tmp=which(tax==x, arr.ind=TRUE)
	#	xcol=max(tmp[,"col"])
	#	dat=as.matrix(tax[xrow<-tmp[,"row"],1:xcol], ncol=xcol)
	#	tab=occurrences[names(occurrences)!=x]
	#	unique(sapply(1:nrow(dat), function(idx){
	#		cur=dat[idx,]
	#		   if(any(cur%in%names(tab)->xx)){
	#				ht=tab[names(tab)%in%cur]
	#				return(names(ht)[max(which(ht==max(ht)))])
	#		   } else {
	#				return(rownames(dat)[idx])
	#		   }
	#	}))
	#})


    clds=f(names(occurrences), function(x) {
		tmp=which(tax==x, arr.ind=TRUE)[,"row"]
        taxinx=rownames(tax)[tmp]
        return(taxinx)
	})
    names(clds)=names(occurrences)

    nn=sapply(clds, length)
    clds=clds[order(nn)]

    ## exclude DUPLICATES
    eq=f(1:length(clds), function(idx){
        others=clds[-which(names(clds)==names(clds)[idx])]
        vv=sapply(others, function(x) setequal(clds[[idx]], x))
        if(any(vv)) return(names(which(vv))) else return()
    })
    drop=c()
    tmp=unique(c(unlist(eq)))
    if(length(tmp)){
        mm=match(tmp, names(clds))
        tmp=tmp[order(mm)]
        mm=mm[order(mm)]
        names(mm)=tmp

        for(i in 1:length(mm)){
            ww=which(sapply(eq, function(x) names(mm[i])%in%x))
            tt=names(clds)[ww]
            a=which(names(clds)==names(mm[i]))
            b=match(tt, names(clds))
            if(any(b<a)) drop=c(drop, names(mm[i]))
        }

        clds=clds[-which(names(clds)%in%drop)]
    }

    resolve_definition=function(taxon, clades){
        cldother=clades[-which(names(clades)==taxon)]
        tips=clades[[taxon]]

        resolver=function(tips, cld){

            members=c()
            if(!length(cld)) {
                members=c(members, tips)
            } else {
                tt=sapply(cld, function(y){
                    all(y%in%tips)
                })

                if(any(tt)){
                    nm=names(cld)[max(which(tt))]
                    members=c(members, nm)
                    ftips=tips[-which(tips%in%cld[[nm]])]
                    fcld=cld[which(tt)]
                    fcld=fcld[!names(fcld)%in%nm]
                    members=c(members, resolver(ftips, fcld))
                } else {
                    members=c(members, tips)
                }

            }
            return(members)
        }

        def=resolver(tips, cldother)
        def
    }

    defs=lapply(1:length(clds), function(idx){
        taxon=names(clds)[idx]
        resolve_definition(taxon, clds)
    })
    names(defs)=names(clds)

    op=options()
    op$expressions=max(op$expressions, 500000)
    options(op)

	tmp=phylo.clades(defs, ncores=ncores)
	phy=tmp[[length(defs)]]
	tt=table(phy$tip.label)
	if(any(tt>1)){
		warning(paste("The following tips occur multiply, suggesting non-monophyly of subtending groups:\n\t", paste(names(tt[tt>1]), collapse="\n\t"), sep=""))
		warning("Offending tips removed from labeled phylogeny")
		phy=.drop.tip(phy, names(tt[tt>1]))
	}
	phy$root.edge=0
	phy
}


name.check <- function(phy, data, data.names = NULL) {

	if (is.null(data.names)) {
		if (is.vector(data)) {
			data.names <- names(data);
		} else {
			data.names <- rownames(data);
		}
	}
	t <- phy$tip.label;
	r1 <- t[is.na(match(t, data.names))];
	r2 <- data.names[is.na(match(data.names, t))];

	r <- list(sort(r1), sort(r2));

	names(r) <- cbind("tree_not_data", "data_not_tree")
	if (length(r1) == 0 && length(r2) == 0) {
		return("OK");
	} else {
		return(r);
	}
}

cherries <- function(phy){
    cache=.cache.descendants(phy)
	N=Ntip(phy)
	nds=which(tabulate(phy$edge[phy$edge[,2]<=N,1])==2)
    if(length(nds)){
        mat=matrix(NA, nrow=length(nds), ncol=2)
        for(i in 1:length(nds)) mat[i, ]=phy$tip.label[cache$tips[[nds[i]]]]
        rownames(mat)=nds
    } else {
        mat=NULL
    }
    return(mat)
}


tips <- function(phy, node)
{
	if(node<=Ntip(phy)) return(phy$tip.label[node])
	dd=.get.descendants.of.node(node, phy, tips=TRUE)
	phy$tip.label[dd]
}

span.phylo=function(phy){
    desc=.cache.descendants(phy)
    N=Ntip(phy)
    labs=c(phy$tip.label, phy$node.label)

    ff=function(node){
        if(length(node)>1) stop("Supply 'node' as a single value")
        if(!is.numeric(node)) node=which(labs==node)

        if(node<=N) return(NULL)
        dd=desc$fdesc[[node]]
        labs[sapply(dd, function(x) desc$tips[[x]][1])]
    }
    ff

}

# data.names is optional, and will replace the names or rownames
# of data when matching data to the tree
# ?!? data.names is not involved here - JWB

# if sort is T, data will have rows in the same order
# as the taxon names in phy$tip.label

treedata <- function(phy, data, sort=FALSE, warnings=TRUE) {

	dm=length(dim(data))

	if (is.vector(data)) {
		data<-as.matrix(data)
	}
	if (is.factor(data)) {
		data<-as.matrix(data)
	}
	if (is.array(data) & length(dim(data))==1) {
		data<-as.matrix(data)
	}

#	if (is.null(data.names)) {
	if (is.null(rownames(data))) {
		stop("names for 'data' must be supplied")
#JME				data.names<-phy$tip.label
#JME				if(warnings)
#JME					cat("Warning: no tip labels, order assumed to be the same as in the tree\n")
	} else {
		data.names<-rownames(data)
	}
#	}
	nc<-name.check(phy, data)
	if (is.na(nc[[1]][1]) | nc[[1]][1]!="OK") {
		if (length(nc[[1]]!=0)) {
			phy=.drop.tip(phy, as.character(nc[[1]]))
			if (warnings) {
				warning(paste("The following tips were not found in 'data' and were dropped from 'phy':\n\t",
							  paste(nc[[1]], collapse="\n\t"), sep=""))
#JME			print(nc[[1]])
#JME			cat("\n")
			}
		}

		if(length(nc[[2]]!=0)) {
			m<-match(data.names, nc[[2]])
			data=as.matrix(data[is.na(m),])
			data.names<-data.names[is.na(m)]
			if(warnings) {
				warning(paste("The following tips were not found in 'phy' and were dropped from 'data':\n\t",
							  paste(nc[[2]], collapse="\n\t"), sep=""))
#JME			print(nc[[2]])
#JME			cat("\n")
			}
		}
 	}
	order<-match(data.names, phy$tip.label)

	rownames(data)<-phy$tip.label[order]

	if (sort) {
    	index <- match(phy$tip.label, rownames(data))
   		data <- as.matrix(data[index,])
	}
	if (dm==2){
		data <- as.matrix(data)
	}

	phy$node.label=NULL

	return(list(phy=phy, data=data))
}

## GENERIC
rescale=function(x, ...) UseMethod("rescale")

## GENERIC

## GENERIC
argn.default=function(x, ...){
	attr(x, "argn")
}

## GENERIC
argn.mkn=function(x, ...){
	attr(x, "argn")
}

## HESSIAN COMPUTATION of CONFIDENCE INTERVALS
.bnd.hessian=function(m, p, s, prob=0.05){
	prob=min(c(1-prob, prob))
	dm=unique(dim(m))
	if(length(dm)!=1) {warning("FAILURE in Hessian CI computation: 'm' must be a square matrix"); return(NA)}
	if(dm!=length(p) | dm!=length(s)) {warning("FAILURE in Hessian CI computation: 'p' and 's' must be of equal length to the dimensions of 'm'"); return(NA)}
	if(!all(s%in%c("exp","nat"))) {warning("FAILURE in Hessian CI computation: 's' must indicate the space of 'p': expecting either 'exp' or 'nat' as elements"); return(NA)}

	qq=qnorm(1-(prob/2))
	dd=try(sqrt(diag(solve(m))),silent=TRUE)
	if(inherits(dd, "try-error")){
		warning(paste("ERROR in inversion of Hessian matrix:", gsub(">", "", gsub("<", "", gsub("\n", "", toString(attributes(dd)$condition))))))
		return(NA)
	}
	bb=qq*dd
	res=sapply(1:length(s), function(idx) {
			   if(s[idx]=="exp"){

			   return(c(lb=exp(log(p[idx])-bb[idx]), ub=exp(log(p[idx])+bb[idx])))
			   } else {
			   return(c(lb=p[idx]-bb[idx], ub=p[idx]+bb[idx]))
			   }
			   })
	rownames(res)=c("lb", "ub")
	colnames(res)=names(p)
	res

}


## MKN CONSTRAINT FUNCTIONS
.er.matrix=function(m){
	m[]=1
	diag(m)=0
	m
}

.sym.matrix=function(m){
	m[]=1:length(m)
	ww=which(upper.tri(m),arr.ind=TRUE)
	m[ww]=m[ww[,2:1]]
	diag(m)=0
	m
}

.ard.matrix=function(m){
	m[]=1:length(m)
	diag(m)=0
	m
}

.meristic.matrix=function(m, symmetric=TRUE){
	m[]=0
	k=nrow(m)
	idx=rbind(cbind(1:(k-1), 2:k), cbind(2:k, 1:(k-1)))
	if(symmetric) base=.sym.matrix(m) else base=.ard.matrix(m)
	m[idx]=base[idx]
	m
}

.repars=function(pars, expected){
    if(!length(pars)==length(expected)) stop(paste("The following 'pars' are expected:\n\t", paste(expected, collapse="\n\t", sep=""), sep=""))
	if(all(!is.null(nm<-names(pars)))){
		if(!all(nm%in%expected)) stop(paste("The following 'pars' are unexpected:\n\t", paste(nm[!nm%in%expected], collapse="\n\t", sep=""), sep=""))
		if(length(unique(nm))!=length(expected)) stop(paste("The following 'pars' are expected:\n\t", paste(expected, collapse="\n\t", sep=""), sep=""))
		mm=match(expected, nm)
		return(pars[mm])
	} else {
		return(pars)
	}
}

## EXTENSION of diversitree:::constrain()
## create 'standard' constraint likelihood function given a patterned matrix
.constrain.k=function(f, model=c("ER", "SYM", "ARD", "meristic"), ...){
	e=environment(f)
	k=e$k

	m=matrix(1:k^2,nrow=k,byrow=TRUE)
	model=match.arg(model,c("ER", "SYM", "ARD", "meristic"))
	mm=switch(model,
			  ER=.er.matrix(m),
			  SYM=.sym.matrix(m),
			  ARD=.ard.matrix(m),
			  meristic=.meristic.matrix(m, ...)
			  )
	return(.constrain.m(f, mm))
}

## EXTENSION of diversitree:::constrain()
## create constrained likelihood function by supplying matrix
.constrain.m=function(f, m){
# if entry is 0, assumed to constrain rate to 0
	k=ncol(m)
	if(nrow(m)!=k) stop("supply 'm' as a square matrix")
	idx <- cbind(rep(1:k, each = k - 1), unlist(lapply(1:k, function(i) (1:k)[-i])))
	npar <- k * (k - 1)
	par=cbind(idx, m[idx], duplicated(m[idx]))
	colnames(par)=c("row","col","parm","dup")
	upar=unique(par[,"parm"])
	sdig=nchar(k)
	sfor=paste("%0",sdig,"d",sep="")
	spar=apply(par[,c("row","col")], 1, function(x) paste("q", paste(sapply(x, function(y) sprintf(sfor, y)), collapse=""),sep=""))
	res=unlist(sapply(1:nrow(par), function(idx){
					  curspar=spar[idx]
					  if(par[idx,"parm"]==0){
					  return(paste(curspar, 0, sep="~") )
					  } else if(par[idx,"dup"]==1){
					  mm=min(which(par[,"parm"]==par[idx,"parm"]))
					  return(paste(curspar, spar[mm], sep="~") )
					  } else {
					  return(NULL)
					  }
					  }))
	return(.constrain(f, formulae=res))
}

# (smart) starting pt for optimization
.ou.smartstart=function(dat, bounds){
	vv=var(dat)
    xb=max(bounds)
    nb=min(bounds)
    atry=seq(-8,4,by=2)
    s=sample(1:length(atry),1)
    if(s==1) {
        aa=nb
    } else if(s==length(atry)) {
        aa=xb
    } else {
        aa=vv*2*exp(atry[s])
        if(aa>xb) aa=xb
        if(aa<nb) aa=nb
    }
    if(is.na(aa)) aa=0.1
    aa
}

# (smart) starting pt for optimization
.bm.smartstart=function(phy, dat){
	ss=mean(pic(dat, phy)^2)
	ss
}

# compute path length from root to tip
.paths.cache=function(cache){

    if(is.null(cache$ordering) || cache$ordering!="postorder"){
        stop("'cache' should be postordered")
    }

	n <- cache$n.tip
	pp <- cache$adesc[-c(1:n)]
	e1 <- cache$edge[, 1]
	e2 <- cache$edge[, 2]
	EL <- cache$edge.length
	xx <- numeric(n + cache$n.node)
	for (i in length(e1):1) {
		var.cur.node <- xx[e1[i]]
		xx[e2[i]] <- var.cur.node + EL[i]
		j <- i - 1L
		while (e1[j] == e1[i] && j > 0) {
			left <- if (e2[j] > n)
			pp[[e2[j] - n]]
			else e2[j]
			right <- if (e2[i] > n)
			pp[[e2[i] - n]]
			else e2[i]
			j <- j - 1L
		}
	}
	xx[1:n]
}

# compute path length from root to tip
.paths.phylo=function(phy, ...){

    ## much from ape:::vcv.phylo()
	phy <- reorder(phy, "postorder")

    FUN=function(vcv=FALSE){
        n <- length(phy$tip.label)
        pp <- .cache.descendants(phy)$tips
        e1 <- phy$edge[, 1]
        e2 <- phy$edge[, 2]
        EL <- phy$edge.length
        xx <- numeric(n + phy$Nnode)
        if(vcv) vmat=matrix(0, n, n)
        for (i in length(e1):1) {
            var.cur.node <- xx[e1[i]]
            xx[e2[i]] <- var.cur.node + EL[i]
            if(vcv){
                j <- i - 1L
                while (e1[j] == e1[i] && j > 0) {
                    left=pp[[e2[j]]]
                    right=pp[[e2[i]]]
                    vmat[left, right] <- vmat[right, left] <- var.cur.node
                    j <- j - 1L
                }
            }
        }
        if(vcv) {
            diags <- 1 + 0:(n - 1) * (n + 1)
            vmat[diags] <- xx[1:n]
            colnames(vmat)<-rownames(vmat)<-phy$tip.label
            return(vmat)
        } else {
            return(xx[1:n])
        }
    }

    FUN(...)
}


.heights.cache=function (cache)
{
    if(is.null(cache$ordering) || cache$ordering!="postorder"){
        stop("'cache' should be postordered")
    }

	n <- cache$n.tip
	n.node <- cache$n.node
	xx <- numeric(n + n.node)
	for (i in nrow(cache$edge):1) xx[cache$edge[i, 2]] <- xx[cache$edge[i, 1]] + cache$edge.length[i]
	root = ifelse(is.null(cache$root.edge), 0, cache$root.edge)
	depth = max(xx)
	tt = depth - xx
	idx = 1:length(tt)
	dd = cache$edge.length[idx]
	mm = match(1:length(tt), c(cache$edge[, 2], n + 1))
	dd = c(cache$edge.length, root)[mm]
	ss = tt + dd
	res = cbind(ss, tt)
	rownames(res) = idx
	colnames(res) = c("start", "end")
	res = data.frame(res)
	res
}

## tree transformation
rescale.phylo=function(x, model=c("BM", "OU", "EB", "nrate", "lrate", "trend", "lambda", "kappa", "delta", "white", "depth"), ...){

	phy=x

	model=match.arg(model, c("BM", "OU", "EB", "nrate", "lrate", "trend", "lambda", "kappa", "delta", "white", "depth"))

	if(!"phylo"%in%class(phy)) stop("supply 'phy' as a 'phylo' object")

	FUN=switch(model,
               BM=.bm.phylo(phy),
			   OU=.ou.phylo(phy),
			   EB=.eb.phylo(phy),
               nrate=.nrate.phylo(phy),
               lrate=.lrate.phylo(phy),
			   trend=.trend.phylo(phy),
			   lambda=.lambda.phylo(phy),
			   kappa=.kappa.phylo(phy),
			   delta=.delta.phylo(phy),
			   white=.white.phylo(phy),
			   depth=.depth.phylo(phy)
			   )
	class(FUN)=c("transformer", "function")
	dots=list(...)
	if(!missing(...)) {
        if(!all(names(dots)%in%argn(FUN))) stop(paste("The following parameters are expected:\n\t", paste(argn(FUN), collapse="\n\t", sep=""), sep=""))
        return(FUN(...))
    } else {
    	return(FUN)
    }
}


# tree transformation
.white.cache=function(cache){
	N=cache$n.tip
	cache$len[]=0
	cache$len[1:N]=1

	z=function(){
		cache
	}
	return(z)
}

# tree transformation
.bm.phylo=function(phy){
    el=phy$edge.length
    z=function(sigsq){
        phy$edge.length=el*sigsq
        phy
    }
    attr(z,"argn")="sigsq"
    return(z)
}

# tree transformation
.white.phylo=function(phy){
	N=Ntip(phy)
	phy$edge.length[]=0
	phy$edge.length[phy$edge[,2]<=N]=1
    el=phy$edge.length
	z=function(sigsq=1){
		phy$edge.length=el*sigsq
        phy
	}
    attr(z,"argn")="sigsq"
	return(z)
}



# tree transformation
.null.cache=function(cache) {
	z=function(){
		cache
	}
	return(z)
}

# tree transformation
.null.phylo=function(phy) {
	z=function(){
		phy
	}
	return(z)
}

# tree transformation -- rates by time
.nrate.phylo=function(phy){
	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t=Tmax-ht$end
	ht$e=ht$start-ht$end
	ht$a=ht$t-ht$e
	ht$rS=ht$a/Tmax
	ht$rE=ht$t/Tmax

	dd=phy$edge[,2]

	relscale.brlen=function(start, end, len, dat){
		ss=start<=dat[,"time"]
		strt=min(which(ss))

		ee=dat[,"time"]<end
		etrt=max(which(ee))+1

		bl=numeric()
        fragment=numeric()
		marker=start
		for(i in strt:etrt){
			fragment=c(fragment, (nm<-(min(c(end, dat[i, "time"]))))-marker)
			bl=c(bl,dat[i, "rate"])
			marker=nm
		}
        fragment=fragment/(sum(fragment))
        sclbrlen=numeric()
        for(i in 1:length(bl)) sclbrlen=c(sclbrlen, len*fragment[i]*bl[i])
		sc=structure(as.numeric(sclbrlen), names=strt:etrt)
		return(sc)
	}


	z=function(time, rate, sigsq=1){
		if(any(time>1) | any(time<0)) stop("supply 'time' as a vector of relative time:\n\tvalues should be in the range 0 (root) to 1 (present)")
		if(any(rate<0)) stop("'rate' must consist of positive values")
		if(length(time)!=length(rate)) stop("'time' and 'rate' must be of equal length")
        phy$edge.length=phy$edge.length*sigsq
		ordx=order(time)
		time=time[ordx]
		rate=rate[ordx]
		dat=cbind(time=c(0,time, 1), rate=(c(1, 1, rate)))
		rs=sapply(dd, function(x) as.numeric(sum(relscale.brlen(ht$rS[x], ht$rE[x], ht$e[x], dat))))
		phy$edge.length=rs
		phy
	}
	attr(z, "argn")=c("time", "rate", "sigsq")
	return(z)
}

# tree transformation -- rates by clade

.lrate.phylo=function(phy){
    N=Ntip(phy)
    n=Nnode(phy)
    cache=.cache.tree(phy)
    cache$phy=phy
    vv=rep(1, N+n-1)

	z=function(node, rate, sigsq=1){


        shifts=c(sort(node[node>N]), node[node<=N])
        mm=match(shifts, node)
        rates=rate[mm]

        phy$edge.length=phy$edge.length*sigsq
        for(i in 1:length(shifts)) vv=.assigndescendants(vv, shifts[i], rates[i], exclude=shifts, cache=cache)
        phy$edge.length=vv*phy$edge.length
        phy
    }

	attr(z, "argn")=c("node", "rate", "sigsq")
	return(z)
}


# tree transformation
.lambda.cache=function(cache){
	N=cache$n.tip
	paths=.paths.cache(cache)

	z=function(lambda){
		if(lambda<0) stop("'lambda' must be positive valued")

		bl=cache$len*lambda
		bl[1:N]=bl[1:N]+(paths-(paths*lambda))
		cache$len=bl
		if(any(cache$len<0,na.rm=TRUE)){
			warning("negative branch lengths encountered:\n\tlambda may be too large")
		}
		cache
	}
	attr(z,"argn")="lambda"
	return(z)
}

# tree transformation
.lambda.phylo=function(phy){
	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:N, phy$edge[,2])
	ht$e=ht$start-ht$end
	paths=.paths.phylo(phy)

	z=function(lambda, sigsq=1){
		if(lambda<0) stop("'lambda' must be positive valued")

		bl=phy$edge.length*lambda
		bl[mm]=bl[mm]+(paths-(paths*lambda))
		phy$edge.length=bl
		if(any(phy$edge.length<0)){
			warning("negative branch lengths encountered:\n\tlambda may be too large")
		}
        phy$edge.length=phy$edge.length*sigsq
		phy
	}
	attr(z,"argn")=c("lambda", "sigsq")
	return(z)
}



# tree transformation
.delta.cache=function(cache){
	ht=.heights.cache(cache)
	N=cache$n.tip
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), cache$edge[,2])
	ht$t=Tmax-ht$end
	ht$e=ht$start-ht$end
	ht$a=ht$t-ht$e

	z=function(delta, rescale=TRUE){
		if(delta<0) stop("'delta' must be positive valued")
		bl=(ht$a+ht$e)^delta-ht$a^delta
		cache$len=bl

		if(rescale){
			scl=Tmax^delta
			cache$len=(cache$len/scl)*Tmax
		}
		cache
	}
	attr(z,"argn")="delta"
	return(z)
}

# tree transformation
.delta.phylo=function(phy){
	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t=Tmax-ht$end
	ht$e=ht$start-ht$end
	ht$a=ht$t-ht$e

	z=function(delta, sigsq=1, rescale=TRUE){
		if(delta<0) stop("'delta' must be positive valued")
		bl=(ht$a+ht$e)^delta-ht$a^delta
		phy$edge.length=bl[phy$edge[,2]]

		if(rescale){
			scl=Tmax^delta
			phy$edge.length=(phy$edge.length/scl)*Tmax
		}
        phy$edge.length=phy$edge.length*sigsq
		phy
	}
	attr(z,"argn")=c("delta", "sigsq")
	return(z)
}



# tree transformation
.trend.cache=function(cache){

	ht=.heights.cache(cache)
	N=cache$n.tip
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), cache$edge[,2])
	ht$head=Tmax-ht$end[cache$edge[mm,1]] # age
	ht$tail=ht$head+(ht$start-ht$end)


	z=function(slope){
		ht$br=1+ht$head*slope
		ht$er=1+ht$tail*slope
        scl=sapply(1:nrow(ht), function(idx){
            if(idx==N+1) return(NA)
            if(ht$br[idx]>0 & ht$er[idx]>0) {
                return((ht$br[idx]+ht$er[idx])/2)
            } else if(ht$br[idx]<0 & ht$er[idx]<0) {
                return(0)
            } else {
                si=-1/slope
                return(ht$br[idx]*(si-ht$head[idx])/(2*(ht$tail[idx]-ht$head[idx])))
            }
        })
		cache$len=cache$len*scl
		cache
	}
	attr(z,"argn")="slope"
	return(z)
}

# tree transformation
.trend.phylo=function(phy){

	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$head=Tmax-ht$end[phy$edge[mm,1]] # age
	ht$tail=ht$head+(ht$start-ht$end)


	z=function(slope, sigsq=1){
        # begin (age): head
        # end: tail
		ht$br=1+ht$head*slope
		ht$er=1+ht$tail*slope
        scl=sapply(1:nrow(ht), function(idx){
            if(idx==N+1) return(NA)
            if(ht$br[idx]>0 & ht$er[idx]>0) {
                return((ht$br[idx]+ht$er[idx])/2)
            } else if(ht$br[idx]<0 & ht$er[idx]<0) {
                return(0)
            } else {
                si=-1/slope
                return(ht$br[idx]*(si-ht$head[idx])/(2*(ht$tail[idx]-ht$head[idx])))
            }
        })

		phy$edge.length=phy$edge.length*scl[phy$edge[,2]]
        phy$edge.length=phy$edge.length*sigsq
		phy
	}
	attr(z,"argn")=c("slope", "sigsq")
	return(z)
}



# tree transformation
.ou.cache=function(cache){
	ht=.heights.cache(cache)
	N=cache$n.tip
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), cache$edge[,2])
	ht$t1=Tmax-ht$end[cache$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1
	z=function(alpha){
		if(alpha<0) stop("'alpha' must be positive valued")
		if (alpha == 0){
      			bl = ht$t2-ht$t1
    		} else {
      			bl = (1/(2 * alpha)) * exp(-2 * alpha * (Tmax - ht$t2)) *
        		-(expm1(-2 * alpha * ht$t2)) - (1/(2 * alpha)) *
        		exp(-2 * alpha * (Tmax - ht$t1)) * -(expm1(-2 *
                                                     		alpha * ht$t1))
    			}
		cache$len=bl
		cache
	}
	attr(z,"argn")="alpha"
	return(z)
}

# tree transformation
.ou.phylo=function(phy){
	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t1=Tmax-ht$end[phy$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1
	z=function(alpha, sigsq=1){
		if(alpha<0) stop("'alpha' must be positive valued")
		bl=(1/(2*alpha))*exp(-2*alpha * (Tmax-ht$t2)) * (1 - exp(-2 * alpha * ht$t2)) - (1/(2*alpha))*exp(-2*alpha * (Tmax-ht$t1)) * (1 - exp(-2 * alpha * ht$t1))
		phy$edge.length=bl[phy$edge[,2]]
        phy$edge.length=phy$edge.length*sigsq
		phy
	}
	attr(z,"argn")=c("alpha", "sigsq")
	return(z)
}



# tree transformation
.eb.cache=function(cache){

	ht=.heights.cache(cache)
	N=cache$n.tip
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), cache$edge[,2])
	ht$t1=Tmax-ht$end[cache$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1

	z=function(a){
		if(a==0) return(cache)
		bl = (exp(a*ht$t2)-exp(a*ht$t1))/(a)
		cache$len=bl
		cache
	}
	attr(z,"argn")="a"
	return(z)
}

# tree transformation
.eb.phylo=function(phy){


	ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t1=Tmax-ht$end[phy$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1

	z=function(a, sigsq=1){
		if(a==0) return(phy)
		bl = (exp(a*ht$t2)-exp(a*ht$t1))/(a)
		phy$edge.length=bl[phy$edge[,2]]
        phy$edge.length=phy$edge.length*sigsq
		phy
	}
	attr(z,"argn")=c("a", "sigsq")
	return(z)
}



# tree transformation
.kappa.cache=function(cache){

	z=function(kappa){
		if(kappa<0) stop("'kappa' must be positive valued")

		cache$len=cache$len^kappa
		cache
	}
	attr(z,"argn")="kappa"
	return(z)
}

# tree transformation
.kappa.phylo=function(phy){

	z=function(kappa, sigsq=1){
		if(kappa<0) stop("'kappa' must be positive valued")

		phy$edge.length=phy$edge.length^kappa
        phy$edge.length=phy$edge.length*sigsq
		phy
	}
	attr(z,"argn")=c("kappa", "sigsq")
	return(z)
}

# tree transformation
.depth.phylo=function(phy){
	orig=max(heights.phylo(phy))
	z=function(depth){
		phy$edge.length <- (phy$edge.length/orig) * depth
		if(!is.null(phy$root.edge)) phy$root.edge=(phy$root.edge/orig) * depth
		phy
	}
	attr(z,"argn")="depth"
	z
}

white.mkn <- function(dat) {
	tt <- table(dat);
	n <- sum(tt);
	p <- tt/n;
	ll <- sum(log(p^tt));

	k <- length(tt) - 1;
	opt <- list(lnL = ll, method = "MLE", k = k);
	opt <- .aic(opt, n);
	#opt <- .aic(lnL = opt$lnL, n = n, k = k); # changed .aic args
	return(opt);
}

drop.extinct <- function (phy, tol=NULL) {
    if (!"phylo" %in% class(phy)) {
        stop("\"phy\" is not of class \"phylo\".");
    }
    if (is.null(phy$edge.length)) {
        stop("\"phy\" does not have branch lengths.");
    }
    if (is.null(tol)) {
        tol <- min(phy$edge.length)/100;
    }
    aa <- is.extinct(phy=phy, tol=tol);
    if (length(aa) > 0) {
        #cat("Dropping", length(aa), "taxa:", aa, "\n", sep=" ");
        phy <- .drop.tip(phy, aa);
    }
    return(phy);
}

# return tip.labels, so that tree ordering is not an issue
is.extinct <- function (phy, tol=NULL) {
    if (!"phylo" %in% class(phy)) {
        stop("\"phy\" is not of class \"phylo\".");
    }
    if (is.null(phy$edge.length)) {
        stop("\"phy\" does not have branch lengths.");
    }
    if (is.null(tol)) {
        tol <- min(phy$edge.length)/100;
    }
    phy <- reorder(phy);
    xx <- numeric(Ntip(phy) + phy$Nnode);
    for (i in 1:length(phy$edge[,1])) {
        xx[phy$edge[i,2]] <- xx[phy$edge[i,1]] + phy$edge.length[i];
    }
    aa <- max(xx[1:Ntip(phy)]) - xx[1:Ntip(phy)] > tol;
    if (any(aa)) {
        return(phy$tip.label[which(aa)]);
    } else {
        return(NULL);
    }
}


drop.random<-function (phy, n)
{
    if (class(phy) != "phylo")
	stop("object \"phy\" is not of class \"phylo\"")
    nb.tip <- length(phy$tip.label)
    if (n > nb.tip)
	return(NULL)
    cut <- sample(1:nb.tip, n)
    r <- .drop.tip(phy, cut)
    return(r)
}
