#=====================================================================================
#
#       Filename:  drop.out.points.R 
#
#    Description:  R function allow you to drop points which is out of n*sigma. Where sigma is standart diviation.
#
#        Version:  1.0
#        Created:  22-Apr-2009
#       Revision:  none
#       
#
#         Author:  Maksim V. Struchalin
#        Company:  ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.
#          Email:  m.struchalin@erasmusmc.nl
#				 license:  GPL (>=2)
#
#=====================================================================================


"drop.out.points" <-
function(trait=c(), n=3) {

if(length(trait) == 0)
	{
	print("Function drop.out.points(...) drop out points which are beyond of n*sigma. Where sigma is standart diviation.")
	print("Usage: drop.out.points(trait), where trait is simle numerical vector, n - number of sigma")
	print("All values which are beyond n*sigma are replaced by NA")
	print("author: M.Struchalin, m.struchalin@erasmusmc.nl")
	print("company: ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.")
	print("license: GPL (>=2)")
	stop("no input parameters")
	}




if(!is.numeric(trait))
	{
	stop("trait must be numeric")
	}




#start exclusion

trait_passed <- c(T,F)

while(dim(table(trait_passed)) != 1)
	{
	trait_mean <- mean(trait, na.rm=T)
	sd <- sqrt(var(trait, na.rm=T))
	trait_passed  <- (trait < trait_mean + sd*n & trait > trait_mean - sd*n)
	trait[!trait_passed] <- NA
	}

trait
}
