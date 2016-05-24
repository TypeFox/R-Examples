input.txt <- function(tab.txt, mod)
	{ 
	
  # object=character (TRUE/FALSE)
  # exports mod of data [,2] numeric == distance data, factor == pop.comparisons

	# tab.txt=normal input format:
	#	individual population locus1.allele1 locus1.allele2 locus2.allele1 locus2.allele2 locus3.allele1 locus3.allele2
	#	p1.1	population1	123	123	233	245	176	145
	#	p1.2	population1	123	123	235	245	123	123
	#	p1.3	population1	145	123	233	245	176	176
	#	.	.		.	.	.	.	.	.

			{tab <- read.table(tab.txt, header=TRUE)}
		

    if (mod=="dist") {message("Information on distances are loaded. Pairwise relatedness is analysed by linear regression.","\n","\n")}
    if (mod=="pop") {message("Column two contains information on populations. Following populations are used for calculations:","\n",levels(tab[,2]),"\n","\n")}
    if (mod=="ref.pop") {message("Custom reference populations are loaded from parameter reference.pop.","\n","\n")}
    return(tab)
    
	}
