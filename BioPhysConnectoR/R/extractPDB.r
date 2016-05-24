#############################################
#   This code is subject to the license as stated in DESCRIPTION.
#   Using this software implies acceptance of the license terms:
#    - GPL 2
#
#   (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2009.
#
#  keul(AT)bio.tu-darmstadt.de
#
#
#  http://www.kay-hamacher.de
#############################################

extractPDB<-function(file.name, verbose=TRUE){

	# ### Begin of the original bio3d functions "atom2xyz", "atom.select" and "read.pdb" as provided in bio3d 1.0-6 under GPL version2 by Grant, Rodrigues, ElSawy, McCammon, Caves, (2006) {Bioinformatics} 22, 2695--2696.
	atom2xyz<-function(num) {
	num3 <- num*3
	c(t(matrix(c(((num3) - 2),
				((num3) - 1),
				(num3)), ncol=3)))
	}
	atom.select<-function(pdb, string=NULL, verbose=TRUE, rm.insert=FALSE) {

	if (missing(pdb)) {
		stop("atom.select: must supply 'pdb' object, e.g. from 'read.pdb'")
	}

	pdb.bounds <- function(nums) {

		# find the 'bounds' (i.e. the
		# start, end and length) of a
		# concetive range of residue
		# or atom numbers.

		nums   <- as.numeric(nums)
		bounds <- nums[1]
		diff.i <- 1; j <- 1
		nums.start <- nums[1]
		store.inds <- NULL
		# also store ind of 1st atom of new res
		for (i in 2:length(nums)) {
		if (nums[i] != nums[j]) { # for resno
			if ((nums[i] - diff.i)!= nums.start) {
			bounds <- c(bounds,nums[i-1],nums[i])
			nums.start <- nums[i]
			diff.i <- 1
			} else { diff.i <- diff.i + 1 }
			store.inds <- c(store.inds,i)
		}
		j<-j+1
		}

		bounds<-c(bounds, nums[length(nums)])
		bounds<-matrix( bounds, ncol=2, byrow=TRUE,
					dimnames=list( NULL, #c(1:(length(bounds)/2),
						c("start","end")) )

		bounds<-cbind(bounds,length=(bounds[,2]-bounds[,1])+1)

		return(list(bounds=bounds,r.ind=store.inds))

	}

	sel.txt2nums <- function (num.sel.txt) {

		# Splitting functions for numbers

		num1<-unlist(strsplit(num.sel.txt, split=","))
		num2<-suppressWarnings( as.numeric(num1) )
		# comma split still may have "10:100" = NA in num2
		tosplit <- num1[ which(is.na(num2)) ]
		num3 <- unlist(strsplit(tosplit, split=":"))
		# pair-up num3 to make num4
		num4<-NULL; i<-1
		while (i < length(num3) ) {
		num4 <- c(num4, as.numeric(num3[i]):
					as.numeric(num3[i+1]) )
		i<-i+2
		}
		# join and order num2 with num4
		return( sort(unique(c(na.omit(num2),num4))) )
	}

	sel.txt2type <- function (type.sel.txt) {

		# Splitting functions for characters

		type1 <- unlist(strsplit(type.sel.txt, split=","))
		# split on coma and remove white space
		return( gsub(" ","",type1) )
	}


	if (is.null(string)) {
		## Summary if called without a selection string
		sum.segid <- unique(pdb$atom[,"segid"])
		sum.chain <- unique(pdb$atom[,"chain"])
		sum.rnum  <- pdb.bounds(pdb$atom[,"resno"])
		sum.resno <- sum.rnum$bounds
		sum.resid <- table(pdb$atom[sum.rnum$r.ind,"resid"])
		sum.eleno <- pdb.bounds(pdb$atom[,"eleno"])$bounds
		sum.elety <- table(pdb$atom[,"elety"])

		cat(" * Structure Summary *",sep="\n")
		cat("---- segid ----",sep="\n");
		print(sum.segid)
		cat("---- chain ----",sep="\n");
		print(sum.chain)
		cat("---- resno ----",sep="\n")
		print(sum.resno)
		cat("---- resid ----",sep="\n");
		print(sum.resid)
		cat("---- eleno ----",sep="\n");
		print(sum.eleno)
		cat("---- elety ----",sep="\n");
		print(sum.elety)
	} else {

	# string shortcuts
		if (string=="calpha" || string=="CA") {
		string= "//////CA/"
		}
		if (string=="cbeta" || string=="CB") {
		string= "//////N,CA,C,O,CB/"
		}
		if (string=="backbone" || string=="back") {
		string= "//////N,CA,C,O/"
		}
		if (string=="all") {
		string= "///////"
		}
		## - Addation Jan 17 2008
		if (string=="h") {
		h.atom <- which( substr(pdb$atom[,"elety"], 1, 1) %in% "H" )
		match <- list(atom=h.atom, xyz=atom2xyz(h.atom))
		class(match) <- "select"
		if(verbose)
			cat(paste(" *  Selected", length(h.atom), "hydrogen atoms *\n"))
		return(match)
		}
		if (string=="noh") {
		noh.atom <- which( !substr(pdb$atom[,"elety"], 1, 1) %in% "H" )
		match <- list(atom=noh.atom, xyz=atom2xyz(noh.atom))
		class(match) <- "select"
		if(verbose)
			cat(paste(" *  Selected", length(noh.atom), "non-hydrogen atoms *\n"))
		return(match)
		}


		## - end addation

		# main function
		sel <- unlist(strsplit(string, split = "/"))

		if (sel[1] == "") { # full selection string (starts with "/")
		sel <- sel[-1]
		if(length(sel) != 6) {
			print("missing elements, should be:")
			print("/segid/chain/resno/resid/eleno/elety/")
		}

		names(sel) <- c("segid","chain", "resno","resid","eleno","elety")
		#print(sel)

		blank <- rep(TRUE, nrow(pdb$atom) )
		sel.inds <- NULL

		# SEGID
		if(sel["segid"] != "") {
			sel.inds <- cbind(sel.inds,
							segid=is.element( pdb$atom[,"segid"],
								sel.txt2type( sel["segid"] )) )
		} else {  sel.inds <- cbind(sel.inds, segid=blank)  }

		# CHAIN
		if(sel["chain"] != "") {
			sel.inds <- cbind(sel.inds,
							chain=is.element( pdb$atom[,"chain"],
								sel.txt2type( sel["chain"] )) )
		} else { sel.inds <- cbind(sel.inds, chain=blank)  }

		# RESNO
		if(sel["resno"] != "") {

			rn <- sel.txt2nums( sel["resno"] )
			if(is.numeric(rn) & length(rn)==0) {
			# check for R object
			rn <- get(gsub(" ","",sel["resno"]))

			}

			sel.inds <- cbind(sel.inds,
							resno=is.element( as.numeric(pdb$atom[,"resno"]),
								rn))
								#sel.txt2nums( sel["resno"] )) )
		} else {  sel.inds <- cbind(sel.inds, resno=blank)  }

		# RESID
		if(sel["resid"] != "") {
			sel.inds <- cbind(sel.inds,
							resid=is.element(pdb$atom[,"resid"],
								sel.txt2type( sel["resid"] )) )
		} else {  sel.inds <- cbind(sel.inds, resid=blank)  }

		# ELENO
		if(sel["eleno"] != "") {
			sel.inds <- cbind(sel.inds,
							eleno=is.element(as.numeric(pdb$atom[,"eleno"]),
								sel.txt2nums( sel["eleno"] )) )
		} else {  sel.inds <- cbind(sel.inds, eleno=blank)  }

		# ELETY
		if(sel["elety"] != "") {
		##  cat( sel["elety"] ,"\n" ) ### glob2rx
		#if(any(i <- grep("*", sel["elety"]))) {
		#  print("WARN: no wild card '*' matching, yet")
		#}

			sel.inds <- cbind(sel.inds,
							elety=is.element(pdb$atom[,"elety"],
								sel.txt2type( sel["elety"] )) )
		} else {  sel.inds <- cbind(sel.inds, elety=blank)  }

		match.inds <- ( (apply(sel.inds, 1, sum, na.rm=TRUE)==6) )

		if (rm.insert) { # ignore INSERT records
			insert <- which(!is.na(pdb$atom[,"insert"]))
			match.inds[insert] <- FALSE
		}
		# return XYZ indices
		xyz.inds <- matrix(1:length( pdb$atom[,c("x","y","z")] ),nrow=3,byrow=FALSE)
		xyz.inds <- as.vector(xyz.inds[,match.inds])

		if (verbose) {
			sel <- rbind( sel, apply(sel.inds, 2, sum, na.rm=TRUE) )
			rownames(sel)=c("Stest","Natom"); print(sel)
			cat(paste(" *  Selected a total of:",sum(match.inds),
					"intersecting atoms  *"),sep="\n")
		}

		match <- list(atom=which(match.inds), xyz=xyz.inds)
		class(match) <- "select"
		return(match)
		}
	}
	}
	read.pdb<-function (file, maxlines=50000, multi=FALSE,rm.insert=FALSE, rm.alt=TRUE, het2atom=FALSE, verbose=TRUE) {

	if(missing(file)) {
		stop("read.pdb: please specify a PDB 'file' for reading")
	}
	if(!is.numeric(maxlines)) {
		stop("read.pdb: 'maxlines' must be numeric")
	}
	if(!is.logical(multi)) {
		stop("read.pdb: 'multi' must be logical TRUE/FALSE")
	}

	# PDB FORMAT v2.0:    colpos,  datatype,    name,      description
	atom.format <- matrix(c(-6,     NA,          NA,       # (ATOM)
							5,     'numeric',   "eleno",   # atom_no
							-1,     NA,          NA,        # (blank)
							4,     'character', "elety",   # atom_ty
							1,     'character', "alt",     # alt_loc
							4,     'character', "resid",   # res_na
							1,     'character', "chain",   # chain_id
							4,     'numeric',   "resno",   # res_no
							1,     'character', "insert",  # ins_code
							-3,     NA,           NA,       # (blank)
							8,     'numeric',   "x",       # x
							8,     'numeric',   "y",       # y
							8,     'numeric',   "z",       # z
							6,     'numeric',   "o",       # o
							6,     'numeric',   "b",       # b
							-6,     NA,           NA,       # (blank)
							4,     'character', "segid"    # seg_id
							), ncol=3, byrow=TRUE,
						dimnames = list(c(1:17), c("widths","what","name")) )

	split.string <- function(x) {
		# split a string 'x'
		x <- substring(x, first, last)
		x[nchar(x) == 0] <- as.character(NA)
		x
	}
	is.character0 <- function(x){length(x)==0 & is.character(x)}

	trim <- function (s) {
		# Remove leading and traling
		# spaces from character strings
		s <- sub("^ +", "", s)
		s <- sub(" +$", "", s)
		s[(s=="")]<-NA
		s
	}


	# finds first and last (substr positions)
	widths <-  as.numeric(atom.format[,"widths"]) # fixed-width spec
	drop.ind <- (widths < 0) # cols to ignore (i.e. -ve)
	widths <- abs(widths)    # absolute vales for later
	st <- c(1, 1 + cumsum( widths ))
	first <- st[-length(st)][!drop.ind] # substr start
	last <- cumsum( widths )[!drop.ind] # substr end

	# read n lines of PDB file
	raw.lines  <- readLines(file, n = maxlines)
	type <- substring(raw.lines,1,6)

	# check number of END/ENDMDL records
	raw.end <- sort(c(which(type == "END"),
						which(type == "ENDMDL")))

	if (length(raw.end) > 1) {
		print("PDB has multiple END/ENDMDL records")
		if (!multi) {
		print("multi=FALSE: taking first record only")
		raw.lines <- raw.lines[ (1:raw.end[1]) ]
		type <- type[ (1:raw.end[1]) ]
		} else {
		print("multi=TRUE: 'read.dcd' will be quicker!")
		}
	}
	if ( length(raw.end) !=1 ) {
		if (length(raw.lines) == maxlines) {
		# have not yet read all the file
		print("You may need to increase 'maxlines'")
		print("check you have all data in $atom")
		}
	}

	# split by record type
	raw.header <- raw.lines[type == "HEADER"]
	raw.seqres <- raw.lines[type == "SEQRES"]
	raw.helix  <- raw.lines[type == "HELIX "]
	raw.sheet  <- raw.lines[type == "SHEET "]
	raw.atom   <- raw.lines[type == "ATOM  "]
	het.atom   <- raw.lines[type == "HETATM"]
	all.atom   <- raw.lines[type %in% c("ATOM  ","HETATM")]

	# also look for "TER" records
	rm(raw.lines)

	if (verbose) {
		if (!is.character0(raw.header)) { cat(" ", raw.header, "\n") }
	}
	seqres <- unlist(strsplit( trim(substring(raw.seqres,19,80))," "))

	helix  <- list(start = as.numeric(substring(raw.helix,22,25)),
					end   = as.numeric(substring(raw.helix,34,37)),
					chain = trim(substring(raw.helix,20,20)),
					type  = trim(substring(raw.helix,39,40)))

	sheet  <- list(start = as.numeric(substring(raw.sheet,23,26)),
					end   = as.numeric(substring(raw.sheet,34,37)),
					chain = trim(substring(raw.sheet,22,22)),
					sense = trim(substring(raw.sheet,39,40)))

	# format ATOM records as a character matrix
	if (het2atom) {
		atom <- matrix(trim(sapply(all.atom, split.string)), byrow=TRUE,
					ncol=nrow(atom.format[ !drop.ind,]),
					dimnames = list(NULL, atom.format[ !drop.ind,"name"]) )
	} else {
		atom <- matrix(trim(sapply(raw.atom, split.string)), byrow=TRUE,
					ncol=nrow(atom.format[ !drop.ind,]),
					dimnames = list(NULL, atom.format[ !drop.ind,"name"]) )
	}

	# Alt records with m[,"alt"] != NA
	if (rm.alt) {
		if ( sum( !is.na(atom[,"alt"]) ) > 0 ) {
		## Edited: Mon May  4 17:41:11 PDT 2009 to cope with
		## both numeric and character ALT records
		first.alt <- sort( unique(na.omit(atom[,"alt"])) )[1]
		cat(paste("   PDB has ALT records, taking",first.alt,"only, rm.alt=TRUE\n"))
		alt.inds <- which( (atom[,"alt"] != first.alt) ) # take first alt only
		if(length(alt.inds)>0)
			atom <- atom[-alt.inds,]
		}
	}
	# Insert records with m[,"insert"] != NA
	if (rm.insert) {
		if ( sum( !is.na(atom[,"insert"]) ) > 0 ) {
		cat("   PDB has INSERT records, removing, rm.insert=TRUE\n")
		insert.inds <- which(!is.na(atom[,"insert"])) # rm insert positions
		atom <- atom[-insert.inds,]
		}
	}
	het <- matrix(trim(sapply(het.atom, split.string)), byrow=TRUE,
					ncol=nrow(atom.format[ !drop.ind,]),
					dimnames = list(NULL, atom.format[ !drop.ind,"name"]) )

	output<-list(atom=atom,
				het=het,
				helix=helix,
				sheet=sheet,
				seqres=seqres,
				xyz=as.numeric(t(atom[,c("x","y","z")])),
				calpha = as.logical(atom[,"elety"]=="CA"))

	class(output) <- "pdb"
	return(output)

	}
	# ### End of bio3d functions

	p<-read.pdb(file.name,maxlines=5000000,verbose=verbose);
	seq<-p$seqres;
	n1<-length(seq);
	pu<-atom.select(p,string="//////CA/",verbose=verbose);
	n<-length(pu$atom);
	if(n!=n1)print("Length of the sequence extracted from SEQRES and the number of CA atoms in ATOM differ.");

	coords<-matrix(data=p$xyz[pu$xyz],ncol=3,byrow=TRUE);
	b<-as.numeric(p$atom[pu$atom,"b"]);
	seq2<-p$atom[pu$atom,"resid"]
	chains<-summary(as.factor(p$atom[pu$atom,5]));
	ret<-list(pdb=p,seq=seq,lseq=n1,lca=n,caseq=seq2,coords=coords,b=b,chains=chains)
	return(ret)
	}
