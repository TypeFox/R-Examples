##' Main function of the RVPedigree package
##'
##' The RVPedigree function is the main function of the RVPedigree
##' used package.
##'
##' Under the hood this function calls \code{\link{ASKAT.region}},
##'     \code{\link{NormalizedASKAT.region}},
##'     \code{\link{VCC1.region}}, \code{\link{VCC2.region}} or
##'     \code{\link{VCC3.region}}, depending on the \code{method}
##'     parameter specified by the user.
##' @title RVPedigree main function
##' @inheritParams read.haplo
##' @param method character, selects the method to use for the
##'     association testing. Can be one of the following:
##' \itemize{
##' \item \code{"ASKAT"} (default)
##' \item \code{"NASKAT"}, normalized ASKAT
##' \item \code{"VCC1"}, VC-C1
##' \item \code{"VCC2"}, VC-C2
##' \item \code{"VCC3"}, VC-C3
##' }
##' @param y vector of phenotype data (one entry per individual), of
##'     length \eqn{n}.
##' @param X matrix of covariates including intercept (dimension:
##'     \eqn{n \times p}, with \eqn{p} the number of covariates)
##' @param Phi Relationship matrix (i.e. twice the kinship matrix); an
##'     \eqn{n \times n} square symmetric positive-definite matrix.
##' @param regions a data frame with details of the genomic regions in
##'     which the association test specified by the \code{method}
##'     parameter should be run. The data frame should have one row
##'     per region and (at least) four columns with the following
##'     names:
##'     \itemize{
##'     \item \code{Name}: Name of the region (e.g. \code{Gene 01})
##'     \item \code{Chr}: Chromosome on which the region is located.
##'     \item \code{StartPos}: The base pair coordinate at which the
##'     region starts
##'     \item \code{EndPos}: The base pair coordinate at which the
##'     region ends.
##'     }
##'     Any other columns will be ignored.
##' @param weights optional numeric vector of genotype weights. If
##'     this option is not specified, the beta distribution is used
##'     for weighting the variants, with each weight given by
##'     \eqn{w_i = dbeta(f_i, 1, 25)^2}, with \eqn{f_i} the minor
##'     allele frequency (MAF) of variant \eqn{i}. This default is the
##'     same as used by the
##'     \href{https://cran.r-project.org/package=SKAT}{\code{SKAT}
##'     package}. This vector is used as the diagonal of the
##'     \eqn{m \times m} matrix \eqn{W}, with \eqn{m} the number of
##'     variants.
##' @param Nperm (integer) The number of permutations to be done to
##'     calculate the empirical p-value if the VCC2 or VCC3 method is
##'     used. For other methods this parameter is ignored (default:
##'     100).
##' @param pvalThreshold (numeric) Threshold for the association
##'     p-value. Regions with a p-value below this threshold will not
##'     be present in the output data frame (default: 0.1).
##' @param VCC3afterVCC1 (logical) Boolean value that indicates
##'     whether the VC-C3 method should automatically be run on the
##'     variants passing the p-value threshold set using the
##'     \code{pvalThreshold} parameter (default: FALSE).
##' @param Ncores (integer) Number of processor (CPU) cores to be used
##'     in parallel when doing running the association analysis. If
##'     the number of regions is larger than the number of cores, then
##'     each region gets to use maximum one core. If the number of
##'     cores is larger than the number of regions and the VCC2 or
##'     VCC3 methods are selected, the remaining cores are distributed
##'     among the regions to parallelize the permutations used to
##'     determine the p-value (default: 1).
##' @return A data frame containing results of the association test
##'     specified by the \code{method} parameter for each region in
##'     the data frame specified by the \code{regions} parameter. The
##'     output data frame contains the following columns:
##'     \itemize{
##'     \item \code{Score.Test}: the score of the given association test
##'     \item \code{P.value}: the p-value of the association test
##'     \item \code{N.Markers}: the number of markers in the region
##'     \item \code{regionname}: Name of the regions/genes on which you
##'     are running the association tests
##'     }
##'     Note that regions that do not contain any genetic variants
##'     will be removed from the output.
##' @author Lennart C. Karssen
##' @export
RVPedigree <- function(method="ASKAT",
                       y=NULL,
                       X=NULL,
                       Phi=NULL,
                       filename=NULL,
                       type="bed",
                       regions=NULL,
                       weights=NULL,
                       Nperm=100,
                       pvalThreshold=0.1,
                       VCC3afterVCC1=FALSE,
                       Ncores=1)
{
    ## Parameter checks
    check_method(method)
    y <- check_pheno(y)
    check_covariates(X, y)
    check_relmatrix(Phi)
    check_files(filename, type)
    check_regions(regions)
    check_weights(weights)
    check_Ncores(Ncores)
    check_VCC3afterVCC1(VCC3afterVCC1, method)

    ## Compute the eigen vectors and eigen values of the relationship
    ## matrix.
    Phi.eig <- eigen(Phi)
    U       <- Phi.eig$vectors
    S       <- Phi.eig$values

    ## Determine the name of the map file based on the file type.
    switch(type,
           bed={
               mapfile <- sub(".bed$", ".bim", filename)
           },
           ped={
               mapfile <- sub(".ped$", ".map", filename)
           }
           )
    map <- readMapFile(mapfile)

    ## Compute the null model (excl. genotypes) first and only once
    message("Estimating null model...")

    switch(method,
           ASKAT={
               H0 <- Estim.H0.ASKAT(y=y, X=X, S=S, U=U)
           },
           NASKAT={
               H0 <- Estim.H0.NormalizedASKAT(y=y, X=X, S=S, U=U)
           },
           VCC1={
               H0 <- Estim.H0.VCC(y=y, X=X, S=S, U=U)
           },
           VCC2={
               H0 <- Estim.H0.VCC(y=y, X=X, S=S, U=U)
           },
           VCC3={
               H0 <- Estim.H0.VCC(y=y, X=X, S=S, U=U)
           }
           )
    message("Null model estimated.")


    ## If we have more regions than cores, don't parallelise the
    ## permutations in the pvalue.VCC2/3 functions. If there are more
    ## CPU cores than regions the left over CPU cores are distributed
    ## amongst the regions. See notes in Bitbucket Issue #14.
    Nregions <- nrow(regions)
    if(Nregions >= Ncores)
    {
        Ncores.pvalue  <- 1
        Ncores.regions <- Ncores
    } else {
        Ncores.regions  <-  Nregions
        Ncores.leftover <- Ncores - Ncores.regions
        Ncores.pvalue   <- 1 + Ncores.leftover %/% Ncores.regions
    }

    if (Ncores > 1)
    {
        message("Using ", Ncores.regions,
                " cores to analyse the regions in parallel")
        if (method == "VCC2" | method == "VCC3")
            message("Using ", Ncores.pvalue,
                    " cores per region for p-value permutations")
    }

    ## Set up parallel environment for distributing the analysis over
    ## multiple regions. By default Ncores=1, so no problems are
    ## expected.
    registerDoParallel(Ncores.regions)

    message("Starting association analysis of the ", Nregions,
            " regions...")
    ## ## Allocate an empty result data frame to which we can rbind the
    ## ## outcomes of each region
    ## results <- data.frame()
    region <- 0
    results <- foreach(region=1:Nregions, .combine='rbind', .verbose=FALSE) %dopar% {
        regionname <- regions[region, "Name"]
        chr        <- regions[region, "Chr"]
        startpos   <- regions[region, "StartPos"]
        endpos     <- regions[region, "EndPos"]

        switch(method,
               ASKAT={
                   region.results <- ASKAT.region(y=y,
                                                  X=X,
                                                  Phi=Phi,
                                                  type=type,
                                                  filename=filename,
                                                  map=map,
                                                  chr=chr,
                                                  startpos=startpos,
                                                  endpos=endpos,
                                                  regionname=regionname,
                                                  U=U,
                                                  S=S,
                                                  RH.Null=H0,
                                                  weights=weights
                                                  )
                   if (!is.na(region.results$P.value) &
                       region.results$P.value <= pvalThreshold) {
                       region.results
                   }
               },
               NASKAT={
                   region.results <- NormalizedASKAT.region(y=y,
                                                            X=X,
                                                            Phi=Phi,
                                                            type=type,
                                                            filename=filename,
                                                            map=map,
                                                            chr=chr,
                                                            startpos=startpos,
                                                            endpos=endpos,
                                                            regionname=regionname,
                                                            U=U,
                                                            S=S,
                                                            RH.Null=H0,
                                                            weights=weights
                                                            )
                   if (!is.na(region.results$P.value) &
                       region.results$P.value <= pvalThreshold) {
                       region.results
                   }
               },
               VCC1={
                   region.results <- VCC1.region(y=y,
                                                 X=X,
                                                 Phi=Phi,
                                                 type=type,
                                                 filename=filename,
                                                 map=map,
                                                 chr=chr,
                                                 startpos=startpos,
                                                 endpos=endpos,
                                                 regionname=regionname,
                                                 U=U,
                                                 S=S,
                                                 RH.Null=H0,
                                                 weights=weights
                                                 )
                   if (!is.na(region.results$P.value) &
                       region.results$P.value <= pvalThreshold) {
                       region.results
                   }
               },
               VCC2={
                   region.results <- VCC2.region(y=y,
                                                 X=X,
                                                 Phi=Phi,
                                                 type=type,
                                                 filename=filename,
                                                 map=map,
                                                 chr=chr,
                                                 startpos=startpos,
                                                 endpos=endpos,
                                                 regionname=regionname,
                                                 U=U,
                                                 S=S,
                                                 RH.Null=H0,
                                                 weights=weights,
                                                 Nperm=Nperm,
                                                 Ncores=Ncores.pvalue
                                                 )
                   if (!is.na(region.results$P.value) &
                       region.results$P.value <= pvalThreshold) {
                       region.results
                   }
               },
               VCC3={
                   region.results <- VCC3.region(y=y,
                                                 X=X,
                                                 Phi=Phi,
                                                 type=type,
                                                 filename=filename,
                                                 map=map,
                                                 chr=chr,
                                                 startpos=startpos,
                                                 endpos=endpos,
                                                 regionname=regionname,
                                                 U=U,
                                                 S=S,
                                                 RH.Null=H0,
                                                 weights=weights,
                                                 Nperm=Nperm,
                                                 Ncores=Ncores.pvalue
                                                 )
                   if (!is.na(region.results$P.value) &
                       region.results$P.value <= pvalThreshold) {
                       region.results
                   }
               },
               {
                   # This part is executed if the method parameter
                   # doesn't match any of the above.
                   stop("Unkown method specified; the method parameter ",
                        "should be one of: 'ASKAT', 'NASKAT', 'VCC1', ",
                        "'VCC2', or 'VCC3'" )
               }
               )                        # End of switch over methods
    }                                   # End of foreach over regions


    ## If the VCC3afterVCC1 parameter is TRUE, automatically run
    ## VCC3 on the output regions of VCC1.
    if (method == "VCC1" & VCC3afterVCC1) {
        message("'VCC3afterVCC1' option set, starting VCC3 run...")
        rowsBelowThres <- which(rownames(results) %in% regions$Name)
        newRegions <- regions[rowsBelowThres, ]
        results <- RVPedigree(method="VCC3",
                              y=y,
                              X=X,
                              Phi=Phi,
                              type=type,
                              filename=filename,
                              regions=newRegions,
                              weights=weights,
                              Nperm=Nperm,
                              pvalThreshold=pvalThreshold,
                              VCC3afterVCC1=FALSE,
                              Ncores=Ncores
                              )
    }

    return(results)
}
