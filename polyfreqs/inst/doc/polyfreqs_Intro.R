## ------------------------------------------------------------------------
# Simulation setup. Feel free to mess around with these values.
library(polyfreqs)
n_ind <- 30
n_loci <- 5
lambda <- 15
error <- 0.01
allele_freqs <- runif(n_loci, 0.1, 0.9)
ploidy <- 4

dat <- sim_reads(allele_freqs, n_ind, lambda, ploidy, error)

## ------------------------------------------------------------------------
names(dat)
dat$genos[1:2,1:5]
dat$tot_read_mat[1:2,1:5]
dat$ref_read_mat[1:2,1:5]

## ---- eval=FALSE, results='hide'-----------------------------------------
#  # Default run. At a minimum the function uses the total reads matrix,
#  # the reference reads matrix and the ploidy level. This should only take a couple
#  # of minutes.
#  
#  out1 <- polyfreqs(dat$tot_read_mat, dat$ref_read_mat, ploidy=4)

## ------------------------------------------------------------------------
itr<-100000
thn<-100
ofile<-"polyfreqs_100k-mcmc.out"

## ---- results='hide'-----------------------------------------------------
out2 <- polyfreqs(dat$tot_read_mat, 
                  dat$ref_read_mat, 
                  ploidy=4,iter=itr, 
                  thin=thn, outfile=ofile)

## ---- warning=FALSE, message=FALSE---------------------------------------
library(coda)

# Read in the MCMC output as a table.
p_table<-read.table("polyfreqs_100k-mcmc.out",header=T,row.names=1)

# Convert the table to an `mcmc` object.
p_post<-mcmc(p_table[251:1000,])

# check that effective sample sizes are greater than 200
sum(effectiveSize(p_post) < 200)

# If any are less than 200, we can see which ones they are
colnames(p_post[,effectiveSize(p_post) < 200])

## ---- fig.width=7--------------------------------------------------------
plot(p_post[,1])

## ---- eval=FALSE---------------------------------------------------------
#  par(ask=T)
#  traceplot(p_post)
#  
#  # Remember to reset the `ask` plotting parameter to FALSE.
#  par(ask=F)

## ---- fig.width=6, fig.height=5, fig.align='center'----------------------
# Take the mean across loci using the apply function.
# Take 25% burn-in as well
multi_het_obs <- apply(out2$het_obs[251:1000,],1,mean)
multi_het_exp <- apply(out2$het_exp[251:1000,],1,mean)

# Check for convergence of the multi-locus estimate
effectiveSize(multi_het_obs)
effectiveSize(multi_het_exp)

# Plot histograms to visually compare the estimates
hist(multi_het_exp, col="blue", main="Heterozygosity", xlab="", ylim=c(0,200))
hist(multi_het_obs, col="red", add=T)
legend(x="topright", 
       c("expected","observed"), 
       col=c("blue","red"), 
       fill=c("blue","red"), bty="n")

# Gets means and 95% highest posterior density (HPD) intervals
list("mean_exp" = mean(multi_het_exp), 
     "95HPD_exp" = quantile(multi_het_exp, c(0.025, 0.975)), 
     "mean_obs" = mean(multi_het_obs), 
     "95HPD_obs" = quantile(multi_het_obs, c(0.025, 0.975)))

## ---- fig.width=7--------------------------------------------------------
# Read in the table using the code below if you haven't already done so.
# p_table <- read.table("polyfreqs_100k-mcmc.out", header=T, row.names=1)

pps <- polyfreqs_pps(as.matrix(p_table[251:1000,]), 
                     dat$tot_read_mat, 
                     dat$ref_read_mat, 
                     ploidy, error)

names(pps)

plot(density(pps$ratio_diff[,1]), main="PP ratio distribution")
abline(v=0)

## ---- eval=FALSE---------------------------------------------------------
#  # Example of running `polyfreqs` with a 10,000 locus matrix that is split in half
#  # and uses the `col_header` parameter to distinguish between runs.
#  
#  polyfreqs(dat$tot_read_mat[,1:5000],
#            dat$ref_read_mat[,1:5000],
#            col_header="run1")
#  
#  polyfreqs(dat$tot_read_mat[,5001:10000],
#            dat$ref_read_mat[,5001:10000],
#            col_header="run2")

