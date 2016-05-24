## ================================
## Recode a sequence object
## ================================

recodef <- function(x, recodes, otherwise=NULL, na=NULL){
	prev_al <- levels(x)
	used_al <- c(unlist(recodes), na)
	## Seek for coding instruction not in previous levels
	inexistant_al <- which(is.na(match(used_al, prev_al)))
	if(length(inexistant_al)>0 && !is.numeric(x)) {
		## Check that the listed inexistant state is not NA
		if(length(inexistant_al)>1 || !is.na(used_al[inexistant_al])) {
			warning(" [!] Some state are not defined in the levels of x: ", paste(used_al[inexistant_al], sep=" "))
		}
	}
	alph <- c(names(recodes), otherwise)
	# if(!is.null(na)){
		# alph <- c(alph, NA)
	# }
	## Adding all previous alphabet that may appear if otherwise==NULL
	if(is.null(otherwise)) {
		alph <- c(alph, prev_al[is.na(match(prev_al, used_al))])
	}
	allAlphabet <- alph
	recodeFunc <- function(fac) {
		if(is.null(otherwise)){
			fac2 <- as.character(fac)
		}else {
			fac2 <- rep(as.character(otherwise), length(fac))
		}
		if(!is.null(na)){
			fac2[fac %in% na] <- NA
		}
		for(newcode in names(recodes)){
			fac2[fac %in% recodes[[newcode]]] <- newcode
		}
		return(factor(fac2, levels=allAlphabet))
	
	}
	return(recodeFunc(x))
}

seqrecode <- function(seqdata, recodes, otherwise=NULL, labels=NULL, cpal=NULL){
	if (!inherits(seqdata,"stslist")) {
		stop(" [!] data is not a sequence object, see seqdef function to create one")
	}
	seqdatar <- seqdata
	prev_al <- alphabet(seqdata)
	used_al <- unlist(recodes)
	## Seek for coding instruction not in previous alphabet
	inexistant_al <- is.na(match(used_al, c(prev_al, attr(seqdata, "nr"), attr(seqdata, "void"))))
	if(sum(inexistant_al)>0) {
		warning(" [!] Some state are not defined in the alphabet of seqdata:", paste(used_al[inexistant_al], sep=" "))
	}
	al <- c(names(recodes), otherwise)
	alph <- NULL
	
	for(a in al){
		if(!a %in% c(attr(seqdata, "void"), attr(seqdata, "nr"))){
			alph <- c(alph, a)
		}
	}
	## Adding all previous alphabet that may appear if otherwise==NULL
	otheralph <- prev_al[is.na(match(prev_al, used_al))]
	if(is.null(otherwise)) {
		alph <- c(alph, otheralph)
	}
	allAlphabet <- c(alph, attr(seqdata, "nr"), attr(seqdata, "void"))
	recodeFunc <- function(fac) {
		if(is.null(otherwise)){
			fac2 <- as.character(fac)
		}else {
			fac2 <- rep(as.character(otherwise), length(fac))
		}
		for(newcode in names(recodes)){
			fac2[fac %in% recodes[[newcode]]] <- newcode
		}
		return(factor(fac2, levels=allAlphabet))
		
	}
	
	for(i in 1:ncol(seqdata)){
		seqdatar[, i] <- recodeFunc(seqdata[,i])
	}
	if(is.null(labels)){
		labels <- alph
	}
	if(is.null(cpal)){
		cpal <- rep("", length(alph))
		for(i in 1:length(alph)){
			if(alph[i] %in% names(recodes)) {
				a <- recodes[[i]][1]
			} else if (!is.null(otherwise)) {
				a <- otheralph[1]
			} else {
				a <- alph[i]
			}
			if(a ==attr(seqdata, "nr")){
				cpal[i] <- attr(seqdata, "missing.color")
			} else if(a==attr(seqdata, "void")){
				warning("[!] no previous color defined for the 'void' (", attr(seqdata, "void"), ") state, setting it to 'white'")
				cpal[i] <- "white"
			}else{
				cpal[i] <- attr(seqdata, "cpal")[match(a, alphabet(seqdata))]
			}
		}
	}
	attr(seqdatar,"alphabet") <- alph
	attr(seqdatar,"labels") <- labels
	attr(seqdatar,"cpal") <- cpal
	return(seqdatar)
	
}