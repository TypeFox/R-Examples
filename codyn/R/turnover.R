#' @title Species Turnover 
#'@description Computes species turnover between time periods as the proportion of species either gained or lost relative to the total number of species observed across both time periods.
#'Includes an option to compute turnover as just the proportion of species gained (i.e., "appearances") or lost (i.e., "disappearances").
#'
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 
#' @param metric The turnover metric to return:
#' \itemize{
#'  \item{total: }{The default metric, calculates summed appearances and disappearances relative to total species richness across both time periods.}
#'  \item{appearance: }{Calculates the number of species that appeared in the second time period relative to total species richness across both time periods.}
#'  \item{disappearance: }{Calculates the number of species that disappeared in the second time period relative to total species richness across both time periods.}
#' }
#' @return The turnover function returns a data frame with the following attributes:
#' \itemize{
#'  \item{turnover: }{A numeric column with the turnover values. The name of this column is the same as the specified metric (default is "total").}
#'  \item{time.var: }{A column containing the second time point; the name and type of this column is the same as the time.var column in the input dataframe.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if replication is specified.}
#' }
#' @details
#' The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' If multiple replicates are included in the data frame, that column should be specified with replicate.var. Each replicate should reflect a single experimental unit - there must be a single abundance value per species within each time point and replicate.
#' @references
#' Cleland, Elsa E., Scott L. Collins, Timothy L. Dickson, Emily C. Farrer, Katherine L. Gross, Laureano A. Gherardi, Lauren M. Hallett, et al. (2013) "Sensitivity of grassland plant community composition to spatial vs. temporal variation in precipitation." Ecology 94, no. 8: 1687-96.
#' @examples 
#'  data(knz_001d)
#'
#'  # Calculate relative total turnover within replicates
#'  total.res <- turnover(df=knz_001d,  
#'                time.var = "year", 
#'                species.var = "species",
#'                abundance.var = "abundance", 
#'                replicate.var="subplot")
#'  
#'  # Calculate relative species appearances within replicates
#'  appear.res <- turnover(df=knz_001d, 
#'                time.var = "year", 
#'                species.var = "species",
#'                abundance.var = "abundance",  
#'                replicate.var="subplot", 
#'                metric="appearance")
#'  
#'  # Calculate relative species disappearances within replicates
#'  disappear.res <- turnover(df=knz_001d, 
#'                time.var = "year", 
#'                species.var = "species",
#'                abundance.var = "abundance", 
#'                replicate.var="subplot", 
#'                metric="disappearance")
#'  
#' @export
turnover <- function(df, time.var, 
                     species.var, 
                     abundance.var, 
                     replicate.var=NA, 
                     metric="total") {
  
  if(is.na(replicate.var)){
    
    # check there unique species x time combinations
    check_single_onerep(df, time.var, species.var)
    
    # calculate turnover
    output <- turnover_allyears(df, time.var, species.var, abundance.var, metric)
    
    
    } else {
      
      # remove unused levels if replicate.var is a factor
      df[replicate.var] <- if(is.factor(df[[replicate.var]])) {
        factor(df[[replicate.var]])
      } else {
        df[replicate.var]
      }
      
      # check unique species x time x replicate combinations
      check_single(df, time.var, species.var, replicate.var)
      
      # sort and apply turnover to all replicates
      df <- df[order(df[[replicate.var]]),]
      X <- split(df, df[replicate.var])
      out <- lapply(X, FUN=turnover_allyears, time.var, species.var, abundance.var, metric)
      ID <- unique(names(out))
      out <- mapply(function(x, y) "[<-"(x, replicate.var, value = y) ,
                  out, ID, SIMPLIFY = FALSE)
      output <- do.call("rbind", out)
    }
  
  # results
  row.names(output) <- NULL
  return(as.data.frame(output))
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.  
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

#' A function to calculate species turnover between years
#'
#' @param df A dataframe containing time, species and abundance columns
#' @param species.var The name of the species column from df
#' @param time.var The name of the time column from df
#' @param abundance.var The name of the abundance column from df
#' @param metric The turnover metric to return; the default, total, returns summed appearances and disappearances relative to total species richness across both years
#' \itemize{
#'  \item{appearance: }{ returns the number of appearances in the second year relative to total species richness across both years }
#'   \item{disappearance: }{ returns the number of disappearances in the second year relative to the total species richness across both years }
#'   }
#' @return output A dataframe containing the specificed turnover metric and year
turnover_allyears <- function(df, 
                              time.var, 
                              species.var, 
                              abundance.var, 
                              metric=c("total", "disappearance","appearance")) {
    
    # allows partial argument matching
    metric = match.arg(metric) 
  
    # check time and abundance are numeric
    check_numeric(df, time.var, abundance.var)
    
    # sort and remove 0s
    df <- df[order(df[[time.var]]),]
    df <- df[which(df[[abundance.var]]>0),]
    
    ## split data by year
    templist <- split(df, df[[time.var]])
    
    ## create consecutive pairs of time points
    t1 <- templist[-length(templist)]
    t2 <- templist[-1]
    
    ## calculate turnover for across all time points
    out <- Map(turnover_twoyears, t1, t2, species.var, metric)
    output <- as.data.frame(unlist(out))
    names(output)[1] = metric
    
    ## add time variable column
    alltemp <- unique(df[[time.var]])
    output[time.var] =  alltemp[2:length(alltemp)]
    
    # results
    return(output)
}

#' A function to calculate species turnover between two years 
#'
#' @param d1 A dataframe containing a species column from one year
#' @param d2 A dataframe containing a species column from the following year
#' @param species.var The name of the species column in d1 and d2
#' @param metric The turnover metric to return; the default, total, returns summed appearances and disappearances relative to total species richness across both years
#' \itemize{
#'  \item{appearance: }{ returns the number of appearances in the second year relative to total species richness across both years }
#'  \item{disappearance: }{ returns the number of disappearances in the second year relative to the total species richness across both years }
#'  }
#' @return output The specificed turnover metric

turnover_twoyears <- function(d1, d2, 
                              species.var, 
                              metric=c("total", "disappearance","appearance")){
  
  # allows partial argument matching
  metric = match.arg(metric)
  
  # create character vectors of unique species from each df
  d1spp <- as.character(unique(d1[[species.var]]))
  d2spp <- as.character(unique(d2[[species.var]]))
  
  # ID shared species
  commspp <- intersect(d1spp, d2spp)
  
  # count number not present in d2
  disappear <- length(d1spp)-length(commspp)
  
  # count number that appear in d2
  appear <- length(d2spp)-length(commspp)
  
  # calculate total richness
  totrich <- sum(disappear, appear, length(commspp))
 
  # output based on metric 
   if(metric == "total"){
        output <- ((appear+disappear)/totrich)
    } else {
  if(metric == "appearance"){
        output <- appear/totrich
    } else {
  if(metric == "disappearance"){
        output <- disappear/totrich
      }
    }
    }
  
  # results
  return(output)
}



