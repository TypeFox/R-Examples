.build_classification=function(species){
	.split_lab=function(label){
		lab=gsub(".", "_", label, fixed=TRUE)
		lab=gsub(" ", "_", lab, fixed=TRUE)
		lab=unlist(strsplit(lab, "_", fixed=TRUE))
		lab=lab[lab!=""]
		lab
	}

	data.frame(genus=sapply(species, function(s) .split_lab(s)[1]), species=species, stringsAsFactors=FALSE)
}


.build_calibrations=function(dat, scion, scion_desc=NULL, tol=0){
#	dat: branching times from stock; rows are 1:(Ntip(stock)+Nnode(stock))
#			time                             hash
#	1001 352.234677 3a4adb7cc0d4a51b9012dfb5615b3d71
#	1002 243.269677 33757769ee61bde8dd5574ae35b47053

#	scion: phylo tree with 'hash' object -- to be scaled from stock

	fetch_spanning=function(phy, nd, desc){
		# desc: a list from 1:(Ntip(phy)+Nnode(phy)) of tips descended from 'nd'
		if(nd<=Ntip(phy)) return(NULL)
		dd=.get.desc.of.node(nd,phy)[1:2]
		tt=sapply(dd, function(x) return(desc[[x]][1]))
		return(phy$tip.label[sort(tt)])
	}

	if(is.null(scion_desc)) scion_desc=.cache.descendants(scion)$tips

	N=Ntip(scion)
	stock_times=dat
	scion_hash=scion$hash
	df=data.frame(MRCA=scion_hash[(N+1):length(scion_hash)], MaxAge=NA, MinAge=NA, taxonA=NA, taxonB=NA, valid=FALSE, stringsAsFactors=FALSE)
	for(i in 1:nrow(df)){
		if(!is.na(hash.cur<-df$MRCA[i])){
			if(hash.cur%in%stock_times$hash){
				node.idx=i+N
				df[i,c("MaxAge","MinAge")]<-age.idx<-stock_times$time[match(hash.cur, stock_times$hash)]
				df[i,c("taxonA","taxonB")]<-taxa.idx<-fetch_spanning(scion, node.idx, scion_desc)
				if(age.idx>tol & all(!is.na(taxa.idx))) df[i,"valid"]=TRUE
			}
		}
	}
	df=df[df$valid,]
	return(df[,-which(names(df)=="valid")])
}

congruify.phylo=function(reference, target, taxonomy=NULL, tol=0, scale=c(NA, "PATHd8"), ncores=NULL){
    ## adding requirement for ncbit
    ## require(ncbit)

    stock=reference
    scion=target
#	stock: a time-calibrated phylogeny with tip-labels that can be treated as an exemplar for clades in 'scion'
#		-- e.g., tip.label in 'stock' might be "Pinus" while in 'scion' we might have "Pinus_cembra"
#		-- tips in 'stock' can be contained in 'scion' (FIXME: is this true?)
#	megaphylogeny: a rooted phylogeny that is to be time-scaled based on 'stock'
#	taxonomy: linkage table between tipsets for 'stock' and 'scion'; if empty, one is attempted to be built by 'scion' labels
#		-- if NULL, 'stock' tips must correspond to tips in 'scion'... e.g., A, B, C in 'stock'; A_sp1, B_sp2, C_sp3 in 'scion'
#		-- rownames of taxonomy must be tips in megaphylogeny

	## functions
	method=match.arg(unname(sapply(scale, toString)), c("NA", "PATHd8"))

	hashes.mapping <- function (phy, taxa, mapping){
	## GENERAL FUNCTION: find hash tag for every edge in tree (using phy$tip.label or predefined set of 'taxa')
	# returns list of hash tags from node 1:(Nnode(phy)+Ntip(phy))
	# taxa: set of species used to create hash keys
	# mapping: named list -- names are tips in 'phy'; elements are tips represented by each tip in 'phy' (and should also be present in 'taxa')
		mapping=mapping[names(mapping)%in%phy$tip.label]
		if(is.null(taxa)) stop("Must supply 'tips'.")
		if(!all(names(mapping)%in%phy$tip.label)) stop("'mapping' must be named list with names in tip labels of 'phy'.")

		mapping=mapping[match(names(mapping), phy$tip.label)]
		descendants <- .cache.descendants(phy)$tips
		hashes <- sapply(descendants, function(desc) .hash.tip(unlist(mapping[desc]), taxa))
		empty=.md5(integer(length(taxa)))
		hashes[hashes==empty]=NA
		phy$hash=hashes
		phy=.uniquify_hashes(phy)

		return(phy)
	}

	times.mapping=function(phy, taxa, mapping){
	#	mapping: named list -- names are tips in 'phy'; elements are tips represented by each tip in 'phy' (and should also be present in 'taxa')
	#	taxa: species that are represented by tips in 'stock'

		# find hash tags for stock 'phy'
		stock=hashes.mapping(phy, taxa, mapping)
		tmp=heights.phylo(stock)
		tmp$hash=stock$hash
		dat=data.frame(time=tmp[,"end"], hash=stock$hash, stringsAsFactors=FALSE)
		dat$hash[1:Ntip(stock)]=NA 	# destroy keys that are associated with tips

		return(list(stock=stock,dat=dat))
	}

	smooth_scion=function(stock, scion, scion_desc, taxa, spp, tol=0.01, method=c("PATHd8", NA)){
		method=match.arg(toString(method), c("NA", "PATHd8"))
		if(!is.ultrametric(stock)) warning("Supplied 'stock' is non-ultrametric.")
		stock_tmp=times.mapping(stock, taxa, spp)
        stock=stock_tmp$stock
        stock_dat=stock_tmp$dat
		calibration=.build_calibrations(stock_dat, scion, scion_desc, tol=tol)
		if(!nrow(calibration)) {
			warning("No concordant branches reconciled between 'stock' and 'scion'; ensure that 'tax' involves rownames found as tip labels in 'scion'")
			return(NA)
		}
		if(method=="PATHd8") {
			phy=PATHd8.phylo(scion, calibration, base=".tmp_PATHd8", rm=FALSE)
			phy$hash=c(rep("", Ntip(phy)), phy$node.label)
			phy$hash[phy$hash==""]=NA
		} else if(method=="NA"){
			phy=NULL
		}

		stock$node.label=stock$hash[(Ntip(stock)+1):max(stock$edge)]
		stock$node.label[is.na(stock$node.label)]=""
		return(list(phy=phy, calibrations=calibration, reference=stock, target=scion))
	}

	## end functions

	## PROCESSING ##
	classification=taxonomy
	unfurl=FALSE

	if(class(stock)=="phylo") {
		stock=list(stock)
		unfurl=TRUE
	}
	if(is.null(classification)) {
		classification=as.data.frame(unique(as.matrix(.build_classification(scion$tip.label)),MARGIN=2))
	}

	tips=unique(unlist(lapply(stock, function(x) x$tip.label)))
	spp=lapply(tips, function(o) {
			   x=rownames(classification)[which(classification==o, arr.ind=TRUE)[,1]]
			   x=x[x%in%scion$tip.label]
	})
	names(spp)=tips
	taxa=unique(unlist(spp))

	scion=hashes.phylo(scion, taxa, ncores)
	scion_desc=.cache.descendants(scion)$tips
    if(is.null(scion$edge.length)) scion$edge.length=numeric(nrow(scion$edge)) ## JME 01302013

	f=lapply
	results=f(1:length(stock), function(i) {
			  phy=stock[[i]]
			  smooth_scion(phy, scion, scion_desc, taxa, spp, tol=tol, method=method)
	})

	if(unfurl) results=results[[1]]
	return(results)
}

## END CONGRUIFICATION FUNCTIONS ##


write.treePL=function(phy, calibrations, nsites, min=1e-4, base="", opts=list(smooth=100, nthreads=8, optad=0, opt=1, cvstart=1000, cviter=3, cvend=0.1, thorough=TRUE)){
#	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB' from .build_calibrations
#	MRCA							MaxAge     MinAge                                  taxonA                                  taxonB
#	c65bacdf65aa29635bec90f3f0447c6e 352.234677 352.234677                          Inga_chartacea             Encephalartos_umbeluziensis
#	d4bc6557dccbd4e8e18b867979f34f8e 243.269677 243.269677                          Inga_chartacea                     Nuphar_sagittifolia
#	5e30c16531694aff7d94da3b35806677 217.632627 217.632627                          Inga_chartacea                  Schisandra_glaucescens

	if(file.exists(inp<-paste(base,"infile",sep="."))) unlink(inp)
	if(file.exists(int<-paste(base,"intree",sep="."))) unlink(inp)

	poss=list(
			  cv="numeric",
			  collapse="boolean",
			  checkconstraints="boolean",
			  cvstart="numeric",
			  cvstop="numeric",
			  cvmultstep="numeric",
			  verbose="boolean",
			  lftemp="numeric",
			  pltemp="numeric",
			  plcool="numeric",
			  lfstoptemp="numeric",
			  plstoptemp="numeric",
			  lfrtstep="numeric",
			  plrtstep="numeric",
			  thorough="boolean",
			  lfiter="integer",
			  pliter="integer",
			  cviter="integer",
			  ldfsimaniter="integer",
			  plsimaniter="integer",
			  cvsimaniter="integer",
			  calcgrad="numeric",
			  paramverbose="boolean",
			  prime="boolean",
			  opt="boolean",
			  optad="boolean",
			  optcvad="boolean",
			  moredetail="boolean",
			  moredetailad="boolean",
			  moredetailcvad="boolean",
			  randomcv="boolean",
			  ftol="numeric",
			  xtol="numeric",
			  mapspace="boolean",
			  nthreads="integer"
	)
	if(length(opts)==0) {
		print(poss)
		stop("No 'opts' specified")
	}

# correct small branch lengths
	z=phy$edge.length[which(phy$edge.length>0)]
	if(any(z<min)){
		scl=min/min(z)
		phy$edge.length=phy$edge.length*scl
	}
	write.tree(phy, file=int)

##	check appropriateness of constraints ##
#	check for 'calibrations' and 'phy' congruence
#	if(!is.null(phy)){
#		check=function(t, phy) all(t%in%phy$tip.label)
#		a=check(calibrations$taxonA, phy)
#		b=check(calibrations$taxonB, phy)

#		if(!all(c(a,b))) stop("Some calibrations not encountered in tree.")
#	}

##	build r8s file
#	calibrations$fixage=ifelse(calibrations$MinAge==calibrations$MaxAge, TRUE, FALSE)
	constraints<-constraintnames<-character(nrow(calibrations))
	for(i in 1:nrow(calibrations)){
		cal=calibrations[i,]
		taxon=cal$MRCA
		desc=c(cal$taxonA, cal$taxonB)

		txt1=ifelse(!is.na(cal$MinAge), paste("min =", taxon, cal$MinAge, sep=" "), "")
		txt2=ifelse(!is.na(cal$MaxAge), paste("max =", taxon, cal$MaxAge, sep=" "), "")
		txt=paste(txt1,txt2,sep="\n")

		constraints[i]=txt
		constraintnames[i]=paste("mrca =", taxon, desc[1], desc[2], sep=" ")
	}
	infile=list(
				tree=paste("treefile = ", int, sep=""),
				ns=paste("numsites = ", nsites, sep=""),
				names=paste(unlist(constraintnames), collapse="\n"),
				mrca=paste(unlist(constraints), collapse="\n"),
				out=paste("outfile = ", paste(base, "dated", "tre", sep="."), sep=""),
				opt=paste(names(opts), opts, sep="=", collapse="\n")
				)

	inp=paste(base,"infile",sep=".")
	writeLines(paste(infile,collapse="\n\n"), con=inp)
	attr(inp, "method")="treePL"
	return(inp)
}

write.r8s=function(phy=NULL, calibrations, base="", blformat=c(lengths="persite", nsites=1, ultrametric="no", round="yes"), divtime=c(method="NPRS", algorithm="POWELL"), describe=c(plot="chrono_description")){
#	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB' from .build_calibrations
#	MRCA							MaxAge     MinAge                                  taxonA                                  taxonB
#	c65bacdf65aa29635bec90f3f0447c6e 352.234677 352.234677                          Inga_chartacea             Encephalartos_umbeluziensis
#	d4bc6557dccbd4e8e18b867979f34f8e 243.269677 243.269677                          Inga_chartacea                     Nuphar_sagittifolia
#	5e30c16531694aff7d94da3b35806677 217.632627 217.632627                          Inga_chartacea                  Schisandra_glaucescens

	if(file.exists(inp<-paste(base,"infile",sep="."))) unlink(inp)

##	check appropriateness of constraints ##
#	check for 'calibrations' and 'phy' congruence
#	if(!is.null(phy)){
#		check=function(t, phy) all(t%in%phy$tip.label)
#		a=check(calibrations$taxonA, phy)
#		b=check(calibrations$taxonB, phy)

#		if(!all(c(a,b))) stop("Some calibrations not encountered in tree.")
#	}

##	build r8s file
#	calibrations$fixage=ifelse(calibrations$MinAge==calibrations$MaxAge, TRUE, FALSE)
	constraints<-constraintnames<-character(nrow(calibrations))
	for(i in 1:nrow(calibrations)){
		cal=calibrations[i,]
		taxon=cal$MRCA
		desc=c(cal$taxonA, cal$taxonB)

		txt=paste(paste("\tfixage taxon=", taxon,sep=""), paste("age=", cal$MinAge, ";\n", sep=""), sep=" ")

		constraints[i]=txt
		constraintnames[i]=paste("\tMRCA", taxon, desc[1], paste(desc[2], ";\n", sep=""), sep=" ")
	}
	#	phy$node.label=NULL
	infile=paste(c(
				"#nexus\n",
				"begin trees;\n",
				paste("tree r8s = ", write.tree(phy), "\n", sep=""),
				"end;\n",
				"begin r8s;\n",
				paste("\tblformat ", paste(names(blformat), blformat, collapse=" ", sep="="), ";\n", sep=""),
				names=paste(unlist(constraintnames), collapse=""),
				mrca=paste(unlist(constraints), collapse=""),
				"\tcollapse;\n",
				paste("\tdivtime ", paste(names(divtime), divtime, collapse=" ", sep="="), ";\n", sep=""),
				paste("\tdescribe ", paste(names(describe), describe, sep="="), ";\n", sep="", collapse=""),
				"end;"
			),collapse="")

	inp=paste(base,"infile",sep=".")
	writeLines(paste(infile,collapse="\n\n"), con=inp)
	attr(inp, "method")="r8s"
	return(inp)
}


write.pathd8=function(phy, calibrations, base=""){
#	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB'

	if(file.exists(inp<-paste(base,"infile",sep="."))) unlink(inp)

	##	check appropriateness of constraints ##
	#	check for 'calibrations' and 'phy' congruence
	check=function(t, phy) all(t%in%phy$tip.label)
	a=check(calibrations$taxonA, phy)
	b=check(calibrations$taxonB, phy)

	if(!all(c(a,b))) stop("Some calibrations not encountered in tree.")

	##	build PATHd8 file
	calibrations$fixage=ifelse(calibrations$MinAge==calibrations$MaxAge, TRUE, FALSE)
	constraints<-constraintnames<-character(nrow(calibrations))
	for(i in 1:nrow(calibrations)){
		cal=calibrations[i,]
		taxon=cal$MRCA
		desc=c(cal$taxonA, cal$taxonB)
		if(cal$fixage) {
			txt=paste("mrca:", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("fixage=", cal$MinAge, ";", sep=""), sep=" ")
		} else {
			txt1=paste("mrca:", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("minage=", cal$MinAge, ";", sep=""), sep=" ")
			txt2=paste("mrca:", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("maxage=", cal$MaxAge, ";", sep=""), sep=" ")
			txt=paste(txt1,txt2,sep="\n")
		}
		constraints[i]=txt
		constraintnames[i]=paste("name of mrca: ", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("name=", cal$MRCA, ";", sep=""), sep=" ")
	}
	phy$node.label=NULL
	infile=list(tree=write.tree(phy),
				mrca=paste(unlist(constraints), collapse="\n"),
				names=paste(unlist(constraintnames), collapse="\n")
				)

	inp=paste(base,"infile",sep=".")
	writeLines(paste(infile,collapse="\n\n"), con=inp)
	attr(inp, "method")="pathd8"
	return(inp)
}

PATHd8.phylo=function(phy, calibrations=NULL, base="", rm=TRUE){
#	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB'
#		-- if NULL, simple ultrametricization of 'phy' is performed

	phy$node.label=NULL
	if(!is.null(calibrations)){
		infile=write.pathd8(phy, calibrations, base)
	} else {
		infile=paste(base, "infile", sep=".")
		write.tree(phy, infile)
	}
	smooth.file=paste(base, "smoothed.tre", sep=".")
	parsed.outfile=paste(base, "pathd8.out", sep=".")
	outfile=paste(base, "pathd8.orig.out", sep=".")
	if(file.exists(outfile)) unlink(outfile)
	if(!system("which PATHd8", ignore.stdout=TRUE)==0) stop("Install 'PATHd8' before proceeding.")
	system(paste("PATHd8 -n", infile, "-pn >", outfile, sep=" "))
	system(paste("grep \"d8 tree\" ", outfile, ">", parsed.outfile, sep=" "))
	smoothed=read.tree(parsed.outfile)
	if(rm & base=="") {
		unlink(parsed.outfile)
		unlink(smooth.file)
		unlink(outfile)
		unlink(infile)
	}
	return(smoothed)
}
