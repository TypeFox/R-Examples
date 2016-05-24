#' Species synchrony
#' @description Calculates the degree synchrony in species abundances within a community over time.
#' Includes the option for two different synchrony metrics. The first, developed by Loreau and de Mazancourt (2008), compares the variance of the aggregated community with the variance of individual components.
#' The second, developed by Gross et al. (2014), compares the average correlation of each individual species with the rest of the aggregated community.
#'
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column
#' @param species.var The name of the species column
#' @param abundance.var The name of the abundance column
#' @param replicate.var The name of the optional replicate column
#' @param metric The synchrony metric to return:
#' \itemize{
#'  \item{"Loreau": }{The default metric, calculates synchrony following Loreau and de Mazancourt (2008).}
#'  \item{"Gross": }{Calculates synchrony following Gross et al. (2014).}
#' }
#' @return The \code{synchrony} function returns a numeric synchrony value unless a replication column is specified in the input data frame.
#' If replication is specified, the function returns a data frame with the following attributes:
#' \itemize{
#'  \item{synchrony: }{A numeric column with the synchrony values.}
#'  \item{replicate.var: }{A column that shares the same name and type as the replicate.var column in the input data frame.}
#' }
#' @details
#' The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' If multiple replicates are included in the data frame, that column should be specified with replicate.var. Each replicate should reflect a single experimental unit - there must be a single abundance value per species within each time point and replicate.
#' @references
#' Gross, Kevin, Bradley J. Cardinale, Jeremy W. Fox, Andrew Gonzalez, Michel Loreau, H. Wayne Polley, Peter B. Reich, and Jasper van Ruijven. (2014) "Species richness and the temporal stability of biomass production: A new analysis of recent biodiversity experiments." The American Naturalist 183, no. 1: 1-12. doi:10.1086/673915.
#'
#' Loreau, Michel, and Claire de Mazancourt. (2008) "Species synchrony and its drivers: Neutral and nonneutral community dynamics in fluctuating environments." The American Naturalist 172, no. 2: E48-66. doi:10.1086/589746.
#' @examples
#' data(knz_001d)
#' synchrony(knz_001d[knz_001d$subplot=="A_1",],
#'                      time.var = "year", 
#'                      species.var = "species",
#'                      abundance.var = "abundance") # for one subplot
#' \dontrun{
#' synchrony(knz_001d, 
#'              time.var = "year", 
#'              species.var = "species",
#'              abundance.var = "abundance",
#'              replicate.var = "subplot") # across all subplots
#'              
#' synchrony(knz_001d,  
#'              time.var = "year", 
#'              species.var = "species",
#'              abundance.var = "abundance",
#'              replicate.var = "subplot", 
#'              metric="Gross") # With Gross et al. (2014) metric.
#' }
#' @importFrom stats aggregate as.formula sd
#' @export
synchrony <- function(df, time.var, 
                      species.var, 
                      abundance.var, 
                      metric = "Loreau", replicate.var = NA) {

  check_numeric(df, time.var, abundance.var)
  
  # remove unused levels if replicate.var is a factor
  if(is.na(replicate.var) == TRUE){
    check_single_onerep(df, time.var, species.var)
        output <- synch_onerep(df, time.var, species.var, abundance.var, metric)
    } else {
      df[replicate.var] <- if(is.factor(df[[replicate.var]])) {
        factor(df[[replicate.var]])
      } else {
        df[replicate.var]
      }

      # check there  is more than one species
      check_multispp(df, species.var, replicate.var)
      
      # check there unique species x time combinations
      check_single(df, time.var, species.var, replicate.var)
      
      # sort and apply synchrony to all replicates
      df <- df[order(df[[replicate.var]]),]
      X <- split(df, df[replicate.var])
      out <- lapply(X, FUN=synch_onerep, time.var, species.var, abundance.var, metric)
      reps <- unique(df[replicate.var])
      output <- cbind(reps, do.call("rbind", out))
      names(output) = c(replicate.var, "synchrony")
      row.names(output) <- NULL
    
      }
    
  # result
  return(output)
  }

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

#' A function to calculate species synchrony over time within one replicate
#'
#' @param df A dataframe containing rep, time, species and abundance columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @param metric The synchrony metric to return. The default, "Loreau", returns synchrony as calculated by Loreau and de Mazancourt 2008.
#'        The alternative, "Gross", returns synchrony as calculated by Gross et al. 2014
#' @return output The degree of species synchrony. If "Loreau", 1 is perfect synchrony and 0 is perfect asynchrony.
#'        If "Gross", 1 is perfect synchrony and -1 is perfect asynchrony.
synch_onerep <- function(df, time.var, species.var, abundance.var,
                         metric = "Loreau") {
    metric = match.arg(metric, choices = c("Loreau", "Gross")) # for partial argument matching

    #remove any species that were never present.
    df <- subset(df, abundance.var > 0)
    
    #fill in 0s
    spplist <- unique(df[species.var])
    yearlist <- unique(df[time.var])
    fulllist <- expand.grid(species.var = spplist[[species.var]],
                            time.var = yearlist[[time.var]])
    # recapture original names
    names(fulllist) = c(species.var, time.var)
    df2 <- merge(df[c(species.var, time.var, abundance.var)], fulllist, all.y = T)
    df2[is.na(df2)] <- 0

    if (metric == "Loreau") {
      #calculate community variance
      XTformula <- as.formula(paste(abundance.var, "~", time.var, sep = ""))
      XT <- aggregate(XTformula, data = df2, sum)
      #do this within rep
      varXT <- var(XT[abundance.var])

      #calculate species variance
      sdSppformula <- as.formula(paste(abundance.var, "~", species.var,
                                       sep = ""))
      sdSpp <- aggregate(sdSppformula, data = df2, sd)
      varSpp <- sum((sdSpp[abundance.var])) * sum(sdSpp[abundance.var])

      #calculate synchrony
      synchrony <- as.numeric(varXT/varSpp)

    } else {
      if (metric == "Gross") {
        corout <- as.data.frame(cbind(species.var = as.character(),
                                      "sppcor" =  as.numeric()))

        # check to see if there are species which do not vary within a subplot
        nonvary <- apply(transpose_community(df, time.var, species.var, abundance.var), 2, sd)
        if(any(nonvary == 0)) {
          warning("One or more species has non-varying abundance within a subplot and has been omitted")

          # remove non-varying species from the spplist for this subplot
          spplist <- data.frame(species =
                      spplist[[1]][is.na(match(as.character(unlist(spplist)), names(nonvary)[which(nonvary == 0)]))]
                        )
          }

        for (i in 1:nrow(spplist)){
          myspp <- as.character(spplist[[1]][i])
          focalspp <- df2[which(df2[species.var] == myspp),]
          com.focalspp <- transpose_community(focalspp, time.var, species.var, abundance.var)
          otherspp <- df2[which(df2[species.var] != myspp),]
          com.otherspp <- transpose_community(otherspp, time.var, species.var, abundance.var)
          agg.otherspp <- rowSums(com.otherspp)
          sppcor <- stats::cor(agg.otherspp, com.focalspp)
          subout <- as.data.frame(cbind(myspp, sppcor))
          names(subout) = c(species.var, "sppcor")
          subout$sppcor <- as.numeric(as.character(subout$sppcor))
          corout <- rbind(corout, subout)
        }
        
        #average correlation for the community
        synchrony <- mean(corout$sppcor)
      }
    }

    return(synchrony)
}

