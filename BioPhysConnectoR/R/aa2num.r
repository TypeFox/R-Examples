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

aa2num<-function(seq,offset=0,code=0:19,verbose=FALSE){
    # ### Begin of the original bio3d function "aa123" as provided in bio3d 1.0-6 under GPL version2 by Grant, Rodrigues, ElSawy, McCammon, Caves, (2006) {Bioinformatics} 22, 2695--2696.
	aa123<-function(aa){

	# convert one-letter IUPAC amino-acid code into
	# three-letter PDB style, for instance "A" into "ALA".

	aa1 <- c("-","X",
			"A","C","D","E","F","G",
			"H","I","K","L","M","N","P","Q",
			"R","S","T","V","W","Y")
	aa3 <- c("---","UNK",
			"ALA", "CYS", "ASP", "GLU", "PHE", "GLY",
			"HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN",
			"ARG", "SER", "THR", "VAL", "TRP", "TYR")

	convert <- function(x) {
		if(is.na(x)) return(NA)
		if (all(x != aa1)) {
		warning("Unknown one letter code for aminoacid")
	#      return(NA)
		return("UNK")
		}
		else {
		return(aa3[which(x == aa1)])
		}
	}
	return(as.vector(unlist(sapply(aa, convert))))
	}
	# ### End of bio3d function

	if(verbose){
		show.code(code,offset)
		cat("\noffset:",offset,"\n\n")
		}
	if(length(unlist(strsplit(seq[1],split="")))==1){
		seq<-aa123(seq)
		}
	seqcode<-function(seq,i=1,code=0:19){
		switch(seq,
	 	  "CYS"={seq=i+code[1]},
		  "MET"={seq=i+code[2]},
		  "PHE"={seq=i+code[3]},
		  "ILE"={seq=i+code[4]},
		  "LEU"={seq=i+code[5]},
		  "VAL"={seq=i+code[6]},
		  "TRP"={seq=i+code[7]},
		  "TYR"={seq=i+code[8]},
		  "ALA"={seq=i+code[9]},
		  "GLY"={seq=i+code[10]},
		  "THR"={seq=i+code[11]},
		  "SER"={seq=i+code[12]},
		  "ASN"={seq=i+code[13]},
		  "GLN"={seq=i+code[14]},
		  "ASP"={seq=i+code[15]},
		  "GLU"={seq=i+code[16]},
		  "HIS"={seq=i+code[17]},
		  "ARG"={seq=i+code[18]},
		  "LYS"={seq=i+code[19]},
		  "PRO"={seq=i+code[20]});
		return(seq)
		}
	ret<-as.numeric(apply(as.array(seq),1,seqcode,offset,code))
	return(ret)
	}
