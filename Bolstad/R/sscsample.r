#' Simple, Stratified and Cluster Sampling
#' 
#' Samples from a fixed population using either simple random sampling,
#' stratitified sampling or cluster sampling.
#' 
#' 
#' @param size the desired size of the sample
#' @param n.samples the number of repeat samples to take
#' @param sample.type the sampling method. Can be one of "simple",
#' "stratified", "cluser" or 1, 2, 3 where 1 corresponds to "simple", 2 to
#' "stratified" and 3 to "cluster"
#' @param x a vector of measurements for each unit in the population. By
#' default x is not used, and the builtin data set sscsample.data is used
#' @param strata a corresponding vector for each unit in the population
#' indicating membership to a stratum
#' @param cluster a corresponding vector for each unit in the population
#' indicating membership to a cluster
#' @return A list will be returned with the following components:
#' \item{samples}{a matrix with the number of rows equal to size and the number
#' of columns equal to n.samples. Each column corresponds to a sample drawn
#' from the population} \item{s.strata}{a matrix showing how many units from
#' each stratum were included in the sample} \item{means}{a vector containing
#' the mean of each sample drawn}
#' @author James M. Curran, Dept. of Statistics, University of Auckland. Janko
#' Dietzsch, Proteomics Algorithm and Simulation,Zentrum f. Bioinformatik
#' Tuebingen Fakultaet f. Informations- und Kognitionswissenschaften,
#' Universitaet Tuebingen
#' @keywords misc
#' @examples
#' 
#' ## Draw 200 samples of size 20 using simple random sampling
#' sscsample(20,200)
#' 
#' ## Draw 200 samples of size 20 using simple random sampling and store the
#' ## results. Extract the means of all 200 samples, and the 50th sample
#' res = sscsample(20,200)
#' res$means
#' res$samples[,50]
#' 
#' @export
sscsample  =  function (size, n.samples, sample.type = "simple",
                        x = NULL, strata = NULL,
                        cluster = NULL){
    ## Written initially by:
    ## James M. Curran,
    ## Dept. of Statistics, University of Auckland
    ## Auckland, New Zealand
    ##
    ## Modified, corrected and improved by:
    ## Janko Dietzsch
    ## Proteomics Algorithm and Simulation
    ## Zentrum f. Bioinformatik Tuebingen
    ## Fakultaet f. Informations- und Kognitionswissenschaften
    ## Universitaet Tuebingen
    ## R. Mark Sharp
    ## Southwest National Primate Center
    ## Southwest Foundation for Biomedical Research

    group.idx.by.name  =  function(name,names.vec,idx){
        return(idx[names.vec == name])
    }

    draw.stratum  =  function(thresholds) {
        r  =  runif(1)
        for (i in 1:length(thresholds))
            if (r < thresholds[i]){
                stratum  =  i; break;
            }
        return(stratum)
    }

    # data(sscsample.data)
    #sscsample.data = sscsample.data

    if (is.null(x))
        x  =  sscsample.data$income

    nx  =  length(x)

    if (size > nx)
        stop("Sample size must be less than population size")

    if (is.null(strata))
        strata  =  sscsample.data$ethnicity

    strata.names  =  unique(strata)
    n.strata  =  length(strata.names)

    if (nx != length(strata))
        stop("The length of the strata and data vectors must be equal")

    if (is.null(cluster))
        cluster  =  sscsample.data$neighborhood

    n.clusters  =  length(unique(cluster))

    if (nx != length(cluster))
        stop("The length of the cluster and data vectors must be equal")

    samples  =  matrix(0, nrow = size, ncol = n.samples)

    if(sample.type == "stratified" | sample.type == 2){
        idx.vec  =  1:nx
        stratified.data  =  lapply(strata.names,group.idx.by.name,
                                  names.vec=strata,idx=idx.vec)
        names(stratified.data)  =  strata.names

        sample.strata.size  =  size * sapply(stratified.data, length) / nx
        sample.strata.units  =  floor(sample.strata.size) ## integer part of units
        sample.strata.fractions  =  sample.strata.size - sample.strata.units ## determine the fractional unit parts for every stratum
        sample.unit.residuals  =  sum(sample.strata.fractions) ## how many remaining units are determined by fractions

        ## prepare the random draw of the residual units
        if (sample.unit.residuals > 0){
            for (i in 2:length(sample.strata.fractions))
                sample.strata.fractions[i]  =  sample.strata.fractions[i] + sample.strata.fractions[i-1]
            sample.strata.thresholds  =  sample.strata.fractions / sample.unit.residuals
        }
    }else if (sample.type == "cluster" | sample.type == 3) {
        cluster.names  =  unique(cluster)
        cluster.names  =  sort(cluster.names) ## clustered.data should be ordered to be useful inside the 'sampling-loop'
        idx.vec  =  1:nx
        clustered.data  =  lapply(cluster.names,group.idx.by.name,names.vec=cluster,idx=idx.vec)
        names(clustered.data)  =  cluster.names
    }
    for (r in 1:n.samples) {
        if (sample.type == "simple" | sample.type == 1){
            sample.idx  =  sample(1:nx, size)
        }else if (sample.type == "stratified" | sample.type == 2){
            for (stratum in 1:n.strata) { ## Sample the whole units from all strata
                if (stratum == 1)
                    sample.idx  =  sample(stratified.data[[stratum]],
                                         sample.strata.units[stratum])
                else
                    sample.idx  =  c(sample.idx, sample(stratified.data[[stratum]],
                                                        sample.strata.units[stratum]))
            }
            if (sample.unit.residuals > 0) { ## Are there fractional parts?
                for (i in 1:sample.unit.residuals) { ## sample the residual units randomly but according the fractions
                    selected.stratum  =  draw.stratum(sample.strata.thresholds)
                    repeat { ## draw a unit that was not already selected
                        draw.idx  =  sample(stratified.data[[selected.stratum]],1)
                        if (! draw.idx %in% sample.idx) break ## Have we already sampled this unit?
                    }
                    sample.idx  =  c(sample.idx, draw.idx)
                }
            }
        }
        else if (sample.type == "cluster" | sample.type == 3) {
            ## This part samples as many clusters as necessary to reach the specified
            ## sampling size.
            sample.idx  =  vector(mode="numeric")
            temp.cluster.names  =  cluster.names
            while (size > length(sample.idx)) {
                sampled.cluster.name  =  sample(temp.cluster.names,1)
                rest  =  size - length(sample.idx)
                if (length(clustered.data[[sampled.cluster.name]]) <= rest) {
                    sample.idx  =  c(sample.idx, clustered.data[[sampled.cluster.name]])
                } else {
                    sample.idx  =  c(sample.idx, sample(clustered.data[[sampled.cluster.name]], rest))
                }
                temp.cluster.names  =  temp.cluster.names[temp.cluster.names != sampled.cluster.name]
            }
        }
        else stop(paste("Unknown sampling sample.type :", sample.type))
        samples[, r]  =  sample.idx
    }
    means  =  rep(0, n.samples)
    s.strata  =  matrix(0, nrow = n.samples, ncol = n.strata)
    sample.out  =  matrix(0, nrow = size, ncol = n.samples)


    for (r in 1:n.samples) {
        idx  =  samples[, r]
        means[r]  =  mean(x[idx])
        for (j in 1:n.strata)
            s.strata[r, j]  =  sum(strata[idx] == strata.names[j])
        sample.out[, r]  =  x[idx]
     }
        
    results = list(samples = samples, s.strata = s.strata, means = means)
    class(results) = "sscsamp"
    
 
    return(results)
}
