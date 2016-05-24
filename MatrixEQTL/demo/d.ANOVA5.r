library("MatrixEQTL");
			
anova.groups = 5;
options(MatrixEQTL.ANOVA.categories = anova.groups);

# Number of columns (samples)
n = 100;
# Number of covariates
nc = 10;
# Generate the standard deviation of the noise
noise.std = 0.1 + rnorm(n)^2;
# Generate the covariates
cvrt.mat = 2 + matrix(rnorm(n*nc), ncol = nc);

# Generate the vectors with single genotype and expression variables
snps.mat = floor(runif(n, min = 0, max = anova.groups));
gene.mat = 1 + (snps.mat==1) + cvrt.mat %*% rnorm(nc) + rnorm(n) * noise.std;

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
  useModel = modelANOVA, 
  errorCovariance = diag(noise.std^2), 
  verbose = TRUE,
  pvalue.hist = FALSE );

# Pull Matrix eQTL results - t-statistic and p-value

fstat = me$all$eqtls$statistic;
pvalue = me$all$eqtls$pvalue;
rez = c( Fstat = fstat, pvalue = pvalue);
# And compare to those from ANOVA in R
{
  cat("\n\n Matrix eQTL: \n"); 
  print(rez);
  cat("\n R anova(lm()) output: \n");
	lmdl = lm( gene.mat ~ cvrt.mat + factor(snps.mat), 
						 weights = 1/noise.std^2 );
	lmout = anova(lmdl)[2, c("F value","Pr(>F)")];
  print( lmout );
}

# Results from Matrix eQTL and "lm" must agree
stopifnot(all.equal(lmout, rez, check.attributes=FALSE));
