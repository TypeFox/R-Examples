#=====================================================================================
#
#       Filename:  convert.snp.affymetrix.R
#
#    Description:  Function for converting affymetrix data to GenABEL raw data format
#
#        Version:  1.0
#        Created:  16-Apr-2008
#       Revision:  none
#				last modification: 16-Apr-2008
#
#         Author:  Maksim V. Struchalin, Yurii S. Aulchenko
#        Company:  ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.
#          Email:  m.struchalin@erasmusmc.nl, i.aoultchenko@erasmusmc.nl
#
#=====================================================================================


"convert.snp.affymetrix" <-
function(dir, map, outfile, skipaffym) {

if (!is(dir,"character"))
	{
	stop("Wrong data class: the first argument should be character")
	}


if (!is(map,"character"))
	{
	stop("Wrong data class: the second argument should be character")
	}


if (!is(outfile,"character"))
	{
	stop("Wrong data class: the third argument should be character")
	}






command <- paste("ls ", dir,sep="")
fileslist <- system(command, intern = T)


files_amount <- length(fileslist)
if(files_amount == 0) stop("\ndirectory ", dir, " is empty or doesn't exsist\n")




alleleID_raw <- alleleID.char2raw()
alleleID_names_char <- names(alleleID_raw)
names(alleleID_raw) <- NULL
alleleID_amount <- length(alleleID_raw)


d <- .C("convert_snp_affymetrix_C", as.character(dir),
				as.character(fileslist),
			 	as.integer(files_amount),
			 	as.character(map),
			 	as.character(outfile),
			 	as.integer(skipaffym),
			 	as.character(alleleID_names_char),
			 	as.raw(alleleID_raw),
			 	as.integer(alleleID_amount))


}
