library("MatrixEQTL");

# Number of columns (samples)
n = 25;

# Number of covariates
nc = 10;

# Generate the standard deviation of the noise
noise.std = 0.1 + rnorm(n)^2;

# Generate the covariates
cvrt.mat = 2 + matrix(rnorm(n*nc), ncol = nc);

# Generate the vectors with single genotype and expression variables
snps.mat = cvrt.mat %*% rnorm(nc) + rnorm(n);
gene.mat = cvrt.mat %*% rnorm(nc) + rnorm(n) * noise.std + 
						1 + 0.5 * snps.mat + snps.mat * cvrt.mat[,nc];
# Create 3 SlicedData objects for the analysis
snps1 = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
gene1 = SlicedData$new( matrix( gene.mat, nrow = 1 ) );
cvrt1 = SlicedData$new( t(cvrt.mat) );

# Produce no output files
filename = NULL; # tempfile()

# Call the main analysis function
me = Matrix_eQTL_main(
	snps = snps1, 
	gene = gene1, 
	cvrt = cvrt1, 
	output_file_name = filename, 
	pvOutputThreshold = 1, 
	useModel = modelLINEAR_CROSS, 
	errorCovariance = diag(noise.std^2), 
	verbose = TRUE,
	pvalue.hist = FALSE );

# Pull Matrix eQTL results - t-statistic and p-value
beta = me$all$eqtls$beta;
tstat = me$all$eqtls$statistic;
pvalue = me$all$eqtls$pvalue;
rez = c(beta = beta, tstat = tstat, pvalue = pvalue);
# And compare to those from the linear regression in R
{
	cat("\n\n Matrix eQTL: \n"); 
	print(rez);
	cat("\n R summary(lm()) output: \n");
	lmdl = lm( gene.mat ~ snps.mat + cvrt.mat + snps.mat*cvrt.mat[,nc],
					   weights = 1/noise.std^2 );
	lmout = tail(summary(lmdl)$coefficients,1)[,c(1,3,4)];
	print( tail(lmout) );
}

# Results from Matrix eQTL and "lm" must agree
stopifnot(all.equal(lmout, rez, check.attributes=FALSE));
