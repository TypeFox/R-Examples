#' @title Mean Rank Shifts
#' @description A measure of the relative change in species rank abundances, which indicates shifts in relative abundances over time (Collins et al. 2008).
#' Mean rank shifts are calculated as the average difference in species' ranks between consecutive time periods, among species that are present across the entire time series.
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column
#' @param species.var The name of the species column
#' @param abundance.var The name of the abundance column
#' @param replicate.var The name of the optional replicate column
#' @return rank_shift returns a data frame with the following columns:
#' \itemize{
#'  \item{time.var_pair: }{A factor column that returns the two time periods being compared, separated by a dash. The name of this column is the same as the time.var column in the input dataframe followed by "_pair".}
#'  \item{MRS: }{A numeric column with the mean rank shift values.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if replication is specified.}
#' }
#' @details
#' The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' If multiple replicates are included in the data frame, that column should be specified with replicate.var. Each replicate should reflect a single experimental unit - there must be a single abundance value per species within each time point and replicate.
#' @references  Collins, Scott L., Katharine N. Suding, Elsa E. Cleland, Michael Batty, Steven C. Pennings, Katherine L. Gross, James B. Grace, Laura Gough, Joe E. Fargione, and Christopher M. Clark.  (2008) "Rank clocks and plant community dynamics." Ecology 89, no. 12: 3534-41.
#' @examples
#'  # Calculate mean rank shifts within replicates
#'  data(knz_001d)
#'
#'  myoutput <- rank_shift(knz_001d,
#'                      time.var = "year",
#'                      species.var = "species",
#'                      abundance.var = "abundance",
#'                      replicate.var = "subplot")
#'
#'  # Calculate mean rank shifts for a data frame with no replication
#'
#'  myoutput_singlerep <- rank_shift(subset(knz_001d, subplot=="A_1"),
#'                            time.var = "year",
#'                            species.var = "species",
#'                            abundance.var = "abundance")
#' @export
rank_shift <- function(df, time.var,
                            species.var,
                            abundance.var,
                            replicate.var = as.character(NA)) {

    # check time and abundance are numeric
    check_numeric(df, time.var, abundance.var)

  if(is.na(replicate.var)) {

      # check there unique species x time combinations
      check_single_onerep(df, time.var, species.var)

      # calculate mean rank shift
      output <- rank_onerep(df, time.var, species.var, abundance.var)

      } else {

        # check there unique species x time x replicate combinations
        check_single(df, time.var, species.var, replicate.var)

        # remove unused levels if replicate.var is a factor, sort
        df[replicate.var] <- if(is.factor(df[[replicate.var]]) == TRUE) {
                factor(df[[replicate.var]])
            } else {
                df[replicate.var]
            }
        df<-df[order(df[[replicate.var]]),]

        # apply mean rank shift to all replicates
        X <- split(df, df[replicate.var])
        out <- (lapply(X, FUN = rank_onerep, time.var, species.var, abundance.var))
        ID <- unique(names(out))
        out <- mapply(function(x, y) "[<-"(x, replicate.var, value = y) ,
                    out, ID, SIMPLIFY = FALSE)
        output<-do.call("rbind", out)
        row.names(output)<-NULL
      }

    # results
    return(output)
}



############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

#' @title Mean Rank Shifts
#' @description A measure of the relative change in species rank abundances, which indicates shifts in relative abundances over time (Collins et al. 2008).
#' Mean rank shifts are calculated as the average difference in species' ranks between consecutive time periods, among species that are present across the entire time series.
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column
#' @param species.var The name of the species column
#' @param abundance.var The name of the abundance column
#' @param replicate.var The name of the optional replicate column
#' @return mean_rank_shift returns a data frame with the following columns:
#' \itemize{
#'  \item{time.var_pair: }{A factor column that returns the two time periods being compared, separated by a dash. The name of this column is the same as the time.var column in the input dataframe followed by "_pair".}
#'  \item{MRS: }{A numeric column with the mean rank shift values.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if replication is specified.}
#' }
#' @details
#' The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' If multiple replicates are included in the data frame, that column should be specified with replicate.var. Each replicate should reflect a single experimental unit - there must be a single abundance value per species within each time point and replicate.
#' @export
mean_rank_shift <- function(df, time.var,
                            species.var,
                            abundance.var,
                            replicate.var = as.character(NA)){
  .Deprecated("rank_shift")

  rank_shift(df = df, time.var = time.var,
             species.var = species.var, abundance.var = abundance.var,
             replicate.var = replicate.var)
}


#' Function for calculating mean rank shifts
#'
#'
#' This is a function that calculates mean rank shifts
#' @param df dataframe of Community dataset. Must be in 'long' format.
#' @param time.var The time variable
#' @param species.var The species variable
#' @param abundance.var The abundance variable
#' @return a dataframe, showing years compared
rank_onerep <- function(df, time.var,
                     species.var,
                     abundance.var) {
    ## split data by year

  for (v in c(time.var, species.var, abundance.var)){
    assertthat::assert_that(assertthat::has_name(df, v))
  }

  yearlist <- split(df, df[[time.var]])

  ## Compare consecutive pairs of years
  y1 <- yearlist[-length(yearlist)]
  y2 <- yearlist[-1]

  commonspp <- Map(df_intersect, y1, y2, dataname = species.var)

  names(commonspp) <- Map(function(x, y) paste0(x, "-", y), names(y1), names(y2))

  abdname1 <- paste0(abundance.var,"1")
  abdname2 <- paste0(abundance.var,"2")
  rank1 <- ""   # Note: initialized rank1 and rank2 simply to eliminate R CMD check NOTE
  rank2 <- ""
  ranknames <- lapply(commonspp, function(x) cbind(x,
                                                   rank1 = rank(x[[abdname1]]),
                                                   rank2 = rank(x[[abdname2]])
  ))

  rankdiff <- lapply(ranknames,
                     function(x) transform(x, abs_ch_rank = abs(rank2 - rank1)))

  MRS <- sapply(rankdiff, function(x) mean(x$abs_ch_rank))

  data.frame(year_pair = names(MRS), MRS, row.names = NULL)
}


#' Create intersected data frames
#'
#' Create intersections.
#' @param df1 A dataframe
#' @param df2 A dataframe
#' @param dataname The name of the column on which the two datasets will be joined and intersected
df_intersect <- function(df1, df2, dataname = "species") {
  commspp <- intersect(df1[[dataname]], df2[[dataname]])

  ## select out the dataname columsn from df1 and df2
  df1dataname<-data.frame(df1[[dataname]])
  names(df1dataname)=dataname
  df2dataname<-data.frame(df2[[dataname]])
  names(df2dataname)=dataname

  ## rename df1 and df2 columns
  df1[[dataname]]<-NULL
  df2[[dataname]]<-NULL
  names(df1) <- paste0(names(df1),"1")
  names(df2) <- paste0(names(df2),"2")

  df1<-cbind(df1, df1dataname)
  df2<-cbind(df2, df2dataname)
  merge(x = df1[df1[[dataname]] %in% commspp, ],
        y = df2[df2[[dataname]] %in% commspp, ])
}



