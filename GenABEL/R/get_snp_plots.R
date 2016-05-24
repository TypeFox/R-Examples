#=====================================================================================
#
#       Filename:  get_snp_plots.R 
#
#    Description:  Function get.snp.distr.plot(...) allow you to make box plots and histograms for traits in each genogroup for given SNP
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


"get.snp.distr.plot" <-
function(df=data.frame(), plotname="", snpname="snpname", chromosome="chr", traitname="trait", filename="") {
#df contains trait and genotypes (0,1,2 only)
#like here
#
#trait snp
#1.2 0
#1.1 1
#0.7 3
#.....

if(dim(df)[1] == 0 || dim(df)[2] == 0)
	{
	print("Function get.snp.distr.plot(...) allow you to make box plots and histograms for traits in each genogroup for given SNP")
	print("Usage: get.snp.distr.plot(df)")
	print("Input parameteres are data.frame object with 2 columns (trait, snp). Other parameters are default.")
	print("If there are no filename then output is on screen.")
	print("author: M.Struchalin, m.struchalin@erasmusmc.nl")
	print("company: ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.")
	print("license: GPL (>=2)")
	stop("no input parameters")
	}



df <- na.omit(df)


if(class(df$trait) != "numeric" || class(df$snp) != "numeric")
	{
	stop("trait and snp variables have to have numeric type")
	}


#Checking of input parameters
snp_types <- sort(as.numeric(na.omit(unique(df$snp))))
snp_types_num <- length(snp_types)



if(snp_types_num<=1)
	{
	stop(paste("SNP contains wrong amount of different types\ntypes amount is ", snp_types_num, sep=""))
	}


id_num <- dim(df)[1]

if(id_num <=0)
	{
	stop(paste("Wrong id amount (=", id_num, ")", sep=""))
	}



if(filename != "")
	{
	bitmap(filename, type="jpeg", height=12, width=12, res=200)
	}



#start drawing



par(mfrow=c(snp_types_num-1,2))

bartlett_list <- list()
list_num <- 1
group <- df$trait[df$snp == 0]
if(length(group) > 1)
	{
	bartlett_list[[list_num]] <- group
	list_num <- list_num + 1
	}




bartlet <- bartlett.test(df$trait ~ df$snp)
pval <- format(bartlet$p.value, digits=4)

main <- paste(plotname, "\n", "bartlett.test pval=", pval, sep="")
xlab <- paste(snpname, ", chr=", chromosome, sep="")
boxplot(df$trait ~ df$snp, main=main, 
				xlab=xlab, ylab=traitname)

for(gen in 1:snp_types_num)
	{
	gen_group <- snp_types[gen]
	var <- var(df$trait[df$snp==gen_group])
	var <- format(var, digits=4)
	main <- paste(plotname, "\n", "var=",var, "\n", table(df$snp==gen_group)["TRUE"], " ids", sep="")
	plot(density(df$trait[df$snp==gen_group]), main=main, xlab=gen_group)
	}





if(filename != "")
	{
	dev.off()
	}






}
