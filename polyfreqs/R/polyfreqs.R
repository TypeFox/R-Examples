#' Bayesian population genomics in autopolyploids
#'
#' \code{polyfreqs} implements a Gibbs sampling algorithm to perform Bayesian inference on the allele frequencies (and other quantities) in a population of autopolyploids.
#' It is the main function for conducting inference with the \code{polyfreqs} package.
#'
#' Data sets run through \code{polyfreqs} must be of class "matrix" with row names representing the names of the individuals sampled.
#' The simplest way to get data into R for running an analysis is to format the total read matrix and reference read matrix as tab delimited text files with the first column containing the individual names and one column after that with the read counts for each locus. These data can then be read in using the \code{read.table} function with the \code{row.names} argument set equal to 1.
#' An optional tab delimited list of locus names can be included as the first row and are treated as column headers for each locus (set \code{header=T} in the \code{read.table} function).
#' When running the \code{polyfreqs}, there are a number of options that control what the function returns.
#' To estimate genotypes and print posterior genotype samples to file, set the \code{genotypes} argument to \code{TRUE} and select a name for the output directory \code{geno_dir} (defaults to "\code{genotypes}").
#' \code{polyfreqs} also prints the current MCMC generation (with a frequency set by the \code{print_freqs} argument) to the R console so that users can track run times.
#' This print can be turned off by setting \code{quiet=TRUE}. More details on using \code{polyfreqs} can be found in the introductory vignette.
#'
#' @author Paul Blischak
#' @param tM Total reads matrix: matrix containing the total number of reads mapping to each locus for each individual.
#' @param rM Reference reads marix: matrix containing the number of reference reads mapping to each locus for each individual.
#' @param ploidy The ploidy level of individuals in the population (must be >= 2).
#' @param iter The number of MCMC generations to run (default=100,000).
#' @param thin Thins the MCMC output by sampling everything \code{thin} generations (default=100).
#' @param print Frequency of printing the current MCMC generation to stdout (default=1000).
#' @param burnin Percent of the posterior samples to discard as burn-in (default=20).
#' @param error The level of sequencing error. A fixed constant (default=0.01).
#' @param genotypes Logical variable indicating whether or not to print the values of the genotypes sampled during the MCMC (default=FALSE).
#' @param geno_dir File path to directory containing the posterior samples of genotypes output by \code{\link{polyfreqs}} (default = "genotypes").
#' @param col_header Optional column header tag for use in running loci in parallel (default="").
#' @param outfile The name of the ouput file that samples from the posterior distribution of allele frequencies are written to (default="polyfreqs-mcmc.out").
#' @param quiet Suppress the printing of the current MCMC generation to stdout (default=FALSE).
#' @return Returns a list of 3 (4 if \code{genotypes=TRUE}) items:
#' \describe{
#'  \item{\code{posterior_freqs}}{A matrix of the posterior samples of allele frequencies. These are also printed to the file with the name given by the \code{outfile} argument.}
#'  \item{\code{map_genotypes}}{If \code{genotypes=TRUE}, then a fourth item will be returned as a matrix containing the maximum \emph{a posteriori} genotype estimates accounting for burn-in.}
#'  \item{\code{het_obs}}{Matrix of posterior samples of observed heterozygosity.}
#'  \item{\code{het_exp}}{Matrix of posterior samples of expected heterozygosity.}
#'  }
#'
#' @references Blischak PD, LS Kubatko and AD Wolfe. Accounting for genotype uncertainty in the estimation of allele frequencies in autopolyploids. \emph{In revision}.
#'
#' @examples
#' data(total_reads)
#' data(ref_reads)
#' polyfreqs(total_reads,ref_reads,4,iter=100,thin=10)
#'
#' @useDynLib polyfreqs
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit quantile rbinom rpois runif
#' @importFrom utils read.table

#' @export
polyfreqs <- function(tM, rM, ploidy, iter=100000, thin=100, burnin=20, print=1000, error=0.01, genotypes=FALSE, geno_dir="genotypes", col_header="", outfile="polyfreqs-mcmc.out", quiet=FALSE){

  # Check that input matrices are valid.
  stopifnot(is.matrix(tM))
  stopifnot(is.matrix(rM))

  if(genotypes){
    d.success<-dir.create(paste("./",geno_dir,"/",sep=""))
    if(!d.success){
      stop("Attempt to make genotypes directory failed.")
    }

    # Print column headers for the output files.
    cat("iter", paste("p_",col_header,"_loc",1:ncol(tM),sep=""), sep="\t", file=outfile)
    cat("\n",file=outfile, append=TRUE)

    ## Initialize genotype matrix from discrete uniform(0,ploidy)
    ## Initialize allele frequency vector from uniform(0,1)
    ## Replace entries in genotype matrix with missing data in
    ## total read and reference read files (tM = 0).
    missing_data<-(tM==0)
    gM_init<-matrix(sample(0:ploidy,nrow(tM)*ncol(tM),replace=TRUE), nrow(tM), ncol(tM))
    gM_init[missing_data]<-NA
    pV_init<-runif(ncol(tM))
    pV_mat <- matrix(NA, nrow=iter/thin, ncol=ncol(tM))
    het_obs <- matrix(NA, nrow=iter/thin, ncol=ncol(tM))
    het_exp <- matrix(NA, nrow=iter/thin, ncol=ncol(tM))



    rnames<-row.names(tM)
    for(i in 1:nrow(tM)){
      cat("iter", paste("g_",col_header,"_loc",1:ncol(tM),sep=""),sep="\t", file=paste("./",geno_dir,"/",rnames[i],"_g-mcmc.out",sep=""))
      cat("\n",file=paste("./",geno_dir, "/", rnames[i],"_g-mcmc.out",sep=""), append=TRUE)
    }

    # Start MCMC

    if(!quiet){
    cat("Starting MCMC...\n\n")
    }

    for(k in 1:iter){

      # Sample from full conditional on genotypes first.
      # Then reassign the values.
      gM<-sample_g(tM, rM, gM_init, pV_init, ploidy, error)
      gM_init<-gM
      gM_init[missing_data]<-NA

      # Sample from the full conditional on allele frequencies
      # given draw from genotypes.
      pV<-sample_p(tM, gM_init, ploidy)
      pV_init<-pV

      # Print every 'thin' generation of the MCMC.
      if(k %% thin == 0){
        cat(k, pV_init, sep="\t",file=outfile, append=TRUE)
        index <- k/thin
        pV_mat[index,] <- pV_init
        cat("\n",file=outfile, append=TRUE)
        print_g(k,gM_init,tM,geno_dir)
        het_obs[index,] <- point_Hobs(gM_init, ploidy)
        het_exp[index,] <- point_Hexp(pV_init, gM_init, ploidy)
      }

      if(k %% print == 0 && !quiet){
        cat("MCMC generation ", k, "\n")
      }
    }

    map_genotypes <- get_map_genotypes(tM, burnin, geno_dir)

    return(list(map_genotypes=map_genotypes,
                posterior_freqs=pV_mat,
                het_obs=het_obs,
                het_exp=het_exp))

  } else {
    # Print column headers for the output files.
    cat("iter", paste("p_",col_header,"_loc",1:ncol(tM),sep=""), sep="\t", file=outfile)
    cat("\n",file=outfile, append=TRUE)

    ## Initialize genotype matrix from discrete uniform(0,ploidy)
    ## Initialize allele frequency vector from uniform(0,1)
    ## Replace entries in genotype matrix with missing data in
    ## total read and reference read files (tM = 0).
    missing_data<-(tM==0)
    gM_init<-matrix(sample(0:ploidy,nrow(tM)*ncol(tM),replace=TRUE), nrow(tM), ncol(tM))
    gM_init[missing_data]=0
    pV_init<-runif(ncol(tM))
    pV_mat <- matrix(NA, nrow=iter/thin, ncol=ncol(tM))
    het_obs <- matrix(NA, nrow=iter/thin, ncol=ncol(tM))
    het_exp <- matrix(NA, nrow=iter/thin, ncol=ncol(tM))


    # Start MCMC
    cat("Starting MCMC...\n\n")
    for(k in 1:iter){

      # Sample from full conditional on genotypes first.
      # Then reassign the values.
      gM<-sample_g(tM, rM, gM_init, pV_init, ploidy, error)
      gM_init<-gM
      gM_init[missing_data]<-NA

      # Sample from the full conditional on allele frequencies
      # given draw from genotypes.
      pV<-sample_p(tM, gM_init, ploidy)
      pV_init<-pV

      # Print every 'thin' generation of the MCMC.
      if(k %% thin == 0){
        cat(k, pV_init, sep="\t",file=outfile, append=TRUE)
        cat("\n",file=outfile, append=TRUE)
        index <- k/thin
        pV_mat[index,] <- pV_init
        het_obs[index,] <- point_Hobs(gM_init, ploidy)
        het_exp[index,] <- point_Hexp(pV_init, gM_init, ploidy)
      }

      if(k %% print == 0 && !quiet){
        cat("MCMC generation ", k, "\n")
      }
    }

    return(list(posterior_freqs=pV_mat,
                het_obs=het_obs,
                het_exp=het_exp))

  }
}
