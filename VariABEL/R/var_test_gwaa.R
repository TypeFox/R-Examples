#=====================================================================================
#
#           Filename:  var.test.gwaa.R
#
#        Description:  VarABEL functions: test each SNP in gwaa data for homogeneity of trait's variance among genotypes.
#
#            Version:  0.1
#            Created:  ---
#           Revision:  none
#  last modification:  26-Apr-2010
#
#             Author:  Maksim V. Struchalin
#        Modified by:  Maksim V. Struchalin, 11-Jan-2009
#            Company:  ErasmusMC, Epidemiology Department, The Netherlands.
#              Email:  m.struchalin@erasmusmc.nl
#            license:  GPL (>=2)
#
#=====================================================================================





#_____________________________________________________________________________________________________________________________
#Available for user function for genome-wide testing variance homogeneity of trait's distribution.
#
"var_test_gwaa" <- function(formula, genodata, phenodata, genodata_info=NULL, testname="svlm", analysis_type="AAvsABvsBB")
{


if(testname == "bartlett") {testname <- 0}
else if(testname == "levene") {testname <- 1}
else if(testname == "likelihood") {testname <- 2}
else if(testname == "kolmogorov_smirnov") {testname <- 3}
else if(testname == "svlm") {testname <- 4}
else {stop(paste(testname, "is unsupportive type of test. Only levene, and svlm are supported."))}


if(analysis_type == "AAvsABvsBB") {analysis_type <- 0}
else if(analysis_type == "AAvsABandBB") {analysis_type <- 1}
else if(analysis_type == "ABvsAAandBB") {analysis_type <- 2}
else if(analysis_type == "BBvsAAandAB") {analysis_type <- 3}
else {stop(paste(analysis_type, "is unsupportive type of analysis. Only AAvsABvsBB, AAvsABandBB, ABvsAAandBB, and BBvsAAandAB are supported."))} 

MAR <- 2
OUT <- "R" 
FUN <- "variance_homogeneity_test_C_wrapper"


idnum <- 0

if(is(formula, "formula"))
	{
	design_matrix <- model.matrix(formula, data = phenodata)
#	design_matrix <- design_matrix_
	trat_name <- as.character(formula[[2]])
	trait <- phenodata[,trat_name]
#	idnum <- length(trait)
	design_matrix_df <- data.frame(design_matrix)

#Put NAs back to the matrix
	if(dim(design_matrix_df)[1] != dim(phenodata)[1])
		{
		design_matrix_df_with_nas <- matrix(rep(NA, dim(phenodata)[1]*dim(design_matrix_df)[2]), ncol=dim(design_matrix_df)[2], nrow=dim(phenodata)[1])
		design_matrix_df_with_nas <- data.frame(design_matrix_df_with_nas)
		formula_vars <- all.vars(formula)
		na_bool <- rep(T, dim(phenodata)[1])
		for(term in formula_vars)
			{
			na_bool <- na_bool & !is.na(phenodata[,term])
			}
		colnames(design_matrix_df_with_nas) <- colnames(design_matrix_df)
		design_matrix_df_with_nas[na_bool,] <- design_matrix_df	
		design_matrix_df <- design_matrix_df_with_nas
		}
			
	design_matrix_df$snp <- 0 #will be filled into iterator
	idnum <- dim(design_matrix_df)[1]
	
	
	}
#	else if(is(formula, "numeric") || is(formula, "integer") || is(formula, "double"))
	else if(is.character(formula)) 
	{
	trait <- phenodata[,formula]
	idnum <- length(trait)
	design_matrix_df <- data.frame(X.Intercept=rep(1,idnum), snp=rep(0,idnum))
	}

#design_matrix_geno_means_df <- data.frame(X.Intercept=rep(1,idnum), snp=rep(0,idnum))


p <- dim(design_matrix_df)[2] #covariates number + 1

gtNrow <- NA
gtNcol <- NA

genodata_info_df <- data.frame()


if(class(genodata) == "snp.data")
	{
	print("reading genotype info...")
	genodata_info_df <- data.frame(name=genodata@snpnames,
				 												 coding=as.character(genodata@coding),
																 strand=as.character(genodata@strand),
																 chromosome=as.character(genodata@chromosome),
				 												 map=genodata@map
																 )
#	gtNrow <- genodata@nsnps
#	gtNcol <- genodata@nids
	gtNrow <- genodata@nids
	gtNcol <- genodata@nsnps
	genodata <- as.raw(genodata@gtps)
  }	
	else if (class(genodata)=="databel") 
	{
	gtNrow <- dim(genodata)[1]
	gtNcol <- dim(genodata)[2]
	genodata <- genodata@data
	
	if(!is.null(genodata_info))
		{	
		if(file.exists(genodata_info))
			{
			print("reading genotype's info...")
			genodata_info_df <- read.table(genodata_info, stringsAsFactors=F, header=T)
#			print(c(dim(genodata_info_df)[1], gtNcol))
			if(dim(genodata_info_df)[1] != gtNcol) 
				{
				print(paste("warning: File ", genodata_info, " contains information about ", dim(genodata_info_df)[1],
																												 " SNPs but genotype data contains ", gtNcol, ". Ignoring of the info file.", sep=""))
				genodata_info_df <- data.frame(name=sub("^", "snp", as.character(1:gtNcol)))
				}		
			}
		else
			{
			print(paste("warning: File ", genodata_info, "does not exist. Trying to run analysis withouth SNP's information..."))
			genodata_info_df <- data.frame(name=sub("^", "snp", as.character(1:gtNcol)))
			}
		}
	else
		{
		genodata_info_df <- data.frame(name=sub("^", "snp", as.character(1:gtNcol)))
		}

	}
	else
 	{
	stop(paste("genodata class not recognized ('",class(genodata),"')",sep=""))
	}


if(testname == 1) #Levene's
	{
	if(is(formula, "formula"))
		{
		trait_regression <- lm(formula, data=phenodata)
		trait <- trait_regression$residuals
		}
	}	


#print(design_matrix_df)
#print(dim(design_matrix_df))


print("Start variance analysis...")
results_C <- .Call("iterator", genodata,
														 as.integer(gtNrow), as.integer(gtNcol),
													 	 as.character(FUN),
														 as.character(OUT), 
														 as.integer(MAR),
														 as.integer(1),
														 as.integer(18), #iterator additional inputa parameters number
														 as.double(trait),
														 as.double(data.matrix(design_matrix_df)),
														 as.double(data.matrix(design_matrix_df)),
														 as.integer(p),
														 as.integer(analysis_type),
														 as.integer(testname),
														 double(p+1),#betas 
														 double(p+1),#se 
														 double(1),#chi2
														 integer(1),#df
														 double(idnum),#residuals
														 double(idnum),#qty 
														 integer(p),#jpvt
														 double(p),#qraux
														 double(2*p),#work
														 double(p*p),#v
														 double(p*p),#x_for_ch2inv
														 integer(idnum)
														 )

print("Variance analysis done.")

#print(genodata_info_df)
#print(results_C)

results_df <- data.frame(genodata_info_df, results_C)
#print(results_df)

cov_names <- colnames(design_matrix_df)

#print(cov_names)

output_column_names <- c(colnames(genodata_info_df), "chisq", "df", "Intercept_effect", "Intercept_sd")


#print(output_column_names)

for(i in 2:p)
	{
	output_column_names <- c(output_column_names, paste(cov_names[i], "_eff", sep=""))
	output_column_names <- c(output_column_names, paste(cov_names[i], "_se", sep=""))
	}

output_column_names <- c(output_column_names, "snp_eff_dispertion")
output_column_names <- c(output_column_names, "snp_se_dispertion")

colnames(results_df) <- output_column_names

if(testname != 4)
	{
	results_df <- results_df[,1:7]
	}

return(results_df)
}




