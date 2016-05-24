
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2015-04-08 10:47:41 +0200 (Wed, 08 Apr 2015) $
# $Rev: 342 $
# ----------------------------------------------------------------

##########################################################################
##-------------------------ExtractWuxDataFrame--------------------------##
##########################################################################

ExtractWuxDataFrame <- function(datain.df,
                                gcm.subset = NULL,
                                rcm.subset = NULL,
                                em.scn.subset = NULL,
                                subreg.subset = NULL,
                                period.subset = NULL,
                                ref.per.subset = NULL,
                                season.subset = NULL,
                                resolution.subset = NULL,
                                corrected.subset = NULL) {

  ## Extracts a subset dataframe according to the WUX dataframe
  ##
  ## Args:
  ##   datain.df: WUX dataframe
  ##   gcm.subset, rcm.subset, em.scn.subset... : String vectors indicating the subset
  ##
  ## Returns:
  ##   dataout.df: A subset of the WUX dataframe
  ##
  ## Example:
  ##   ## read wux test data
  ##   wux.df <- data(ensembles)
  ##
  ##   ## extract dataframe
  ##   subset.df <- ExtractDataFrame(wux.df, subreg.subset = c("Austria.", "GAR"),
  ##     gcm.subset = c("ECHAM5-r3", "IPSL"))
  ##
  ## History:
  ##   2011-07-14 | original code (geh)
  ##   2014-11-20 | all subset commands changed to '[' operation (thm)

  ## gcm subset
  if (is.null(gcm.subset)) {
    gcm.subset = levels(datain.df$gcm)
  }

  ## rcm subset
  if (is.null(rcm.subset)) {
    rcm.subset = levels(datain.df$rcm)
  }

  ## emission scenario subset
  if (is.null(em.scn.subset)) {
    em.scn.subset = levels(datain.df$em.scn)
  }

  ## subregion subset
  if (is.null(subreg.subset)) {
    subreg.subset = levels(datain.df$subreg)
  }

  ## period subset
  if (is.null(period.subset)) {
    period.subset = levels(datain.df$period)
  }

  ## reference period subset
  if (is.null(ref.per.subset)) {
    ref.per.subset = levels(datain.df$ref.per)
  }

  ## period subset
  if (is.null(period.subset)) {
    period.subset = levels(datain.df$period)
  }

  ## season subset
  if (is.null(season.subset)) {
    season.subset = levels(datain.df$season)
  }

  ## resolution subset
  if (is.null(resolution.subset)) {
    resolution.subset = levels(datain.df$resolution)
  }

  ## corrected subset
  if (is.null(corrected.subset)) {
    corrected.subset = levels(datain.df$corrected)
  }


  dataout.df <- datain.df[datain.df[["gcm"]] %in% gcm.subset &
                          datain.df[["rcm"]] %in% rcm.subset &
                          datain.df[["em.scn"]] %in% em.scn.subset &
                          datain.df[["subreg"]] %in% subreg.subset &
                          datain.df[["period"]] %in% period.subset &
                          datain.df[["ref.per"]] %in% ref.per.subset &
                          datain.df[["season"]] %in% season.subset &
                          datain.df[["resolution"]] %in% resolution.subset &
                          datain.df[["corrected"]] %in% corrected.subset       , , drop = TRUE]

  return(dataout.df)
}
##########################################################################


##########################################################################
##-------------------------ExpandWuxDataFrame---------------------------##
##########################################################################

ExpandWuxDataFrame <- function(datain.df,
                               factor1.name,
                               factor2.name) {

  ## Expands WUX dataframe according to two factors (e.g. RCM and GCM) to get
  ## a balanced matrix-like dataframe. Missing data are indicated as NA.
  ##
  ## Args:
  ##   datain.df: WUX dataframe
  ##   factor1.name: name of the 1st factor
  ##   factor2.name: name of the 2nd factor
  ##
  ## Returns:
  ##   dataout.df: WUX dataframe with all possible combinations according to the
  ##     two indicated factors
  ##
  ## Example:
  ##   ## read wux test data
  ##   wux.df <- data(ensembles)
  ##
  ##   ## expand dataframe
  ##   expanded.df <- ExpandWuxDataFrame(wux.df, "gcm", "rcm")
  ##
  ## History:
  ##   2010-10-27 | original code (thm)
  ##   2011-07-14 | generalisation (geh)

  factor1.levels <- levels(datain.df[ ,factor1.name])
  factor2.levels <- levels(datain.df[ ,factor2.name])

  expanded.df <- expand.grid(factor1.levels, factor2.levels)
  names(expanded.df) <- c(factor1.name, factor2.name)

  dataout.df <- merge(datain.df, expanded.df, all.y = TRUE)

  return(dataout.df)
}
##########################################################################


##########################################################################
##--------------------------reconstruct methos -------------------------##
##########################################################################

reconstruct <- function(x, factor1.name, factor2.name, data.name,
                        method = "LES", iterations.num = 100){
  ## method: LES, Iterative, IterativeCC
  if (method == "LES") {
    rwux.df <- AnovaReconstructLES(x, factor1.name, factor2.name, data.name)
  } else if (method == "Iterative") {
    rwux.df <- AnovaReconstructIterative(x, factor1.name, factor2.name,
                                         data.name, iterations.num = 100)
  } else if (method == "IterativeCC") {
    rwux.df <- AnovaReconstructIterativeCC(x, factor1.name, factor2.name,
                                           data.name, iterations.num = 100)
  }  else {
    stop("UNKNOWN METHOD for \"reconstruct\"")
  }

  class(rwux.df) <- append("wux.df", class(rwux.df))
  class(rwux.df) <- append("rwux.df", class(rwux.df))
  return(rwux.df)
}

##########################################################################
##--------------------------AnovaReconstructLES-------------------------##
##########################################################################

AnovaReconstructLES <- function(datain.df,
                                factor1.name,
                                factor2.name,
                                data.name) {

  ## Performs a simple missing value reconstruction with two factors. The
  ## algorithm follows Deque et al. (2007) but the reconstruction is based
  ## on solving the linear equation system (LES) of the ANOVA instead of
  ## reconstructing iteratively. The main advantages of this method are that
  ## it is much faster and can be more easily extended to more factors than
  ## the original one. However, keep in mind that the results slightly differ
  ## from the iterative procedure proposed by Deque et al. (2007).
  ##
  ## Args:
  ##   datain.df: WUX dataframe
  ##   factor1.name: name of the 1st factor
  ##   factor2.name: name of the 2nd factor
  ##   data.name: name of the variable to be reconstructed
  ##
  ## Returns:
  ##   dataout.df: WUX dataframe containing the reconstructed values
  ##
  ## Example:
  ##   ## read wux test data
  ##   wux.df <- data(ensembles)
  ##
  ##   ## reconstruction of the  missing data
  ##   reconstruct.df <- AnovaReconstructLES(wux.df, factor1.name = "rcm", factor2.name = "gcm",
  ##     data.name = "perc.delta.precipitation_amount")
  ##
  ## History:
  ##   2010-10-27 | original code (geh)
  ##   2014-11-20 | all subset commands changed to '[' operation (thm)
  ##   2014-11-20 | changed hardcoded "eval" statement to explicit model.matrix
  ##                call with a formula object (thm)

  factor1.name <- factor1.name
  factor2.name <- factor2.name
  data.name <- data.name
  datain.df <- droplevels(datain.df)
  season.levels <- levels(datain.df$season)
  subreg.levels <- levels(datain.df$subreg)

  dataout.df <- data.frame()
  for (ii in subreg.levels){
    for (jj in season.levels){

      ## extract data frame
      subset.df <- datain.df[datain.df[["subreg"]] == ii &
                             datain.df[["season"]] == jj, ]
      subset.df <- droplevels(subset.df)

      ## expand data frame
      expanded.df <- ExpandWuxDataFrame(subset.df, factor1.name, factor2.name)
      expanded.df$season <- factor(levels(expanded.df$season))
      expanded.df$subreg <- factor(levels(expanded.df$subreg))

      ## create the design matrix for two-way ANOVA without interaction
      contrasts.arg <- "contr.sum"
      designmatrix.formula <- formula(paste(
                                        "~ ",
                                        factor1.name,
                                        " + ",
                                        factor2.name,
                                        sep = ""))
      contrast.list <- list(contrasts.arg, contrasts.arg)
      names(contrast.list) <- c(factor1.name, factor2.name)
      design.matrix <- model.matrix(designmatrix.formula,
                                    expanded.df,
                                    contrasts = contrast.list
                                    )
      ## previous implicit design matrix call (changed 2014-11-20 by thm)
      ## design.matrix.expr <- paste("design.matrix <- model.matrix(~ ",
      ##                             factor1.name, " + ", factor2.name,
      ##                             ", expanded.df, contrasts = list(",
      ##                             factor1.name, " = ", contrasts.arg,
      ##                             ", ", factor2.name," = ", contrasts.arg,
      ##                             "))",  sep = "")
      ## eval(parse(text = design.matrix.expr))

      ## solve the linear equation system via the pseudoinverse
      data.vector <- expanded.df[ ,data.name]
      model.vector <- as.vector(design.matrix %*%
                                pseudoinverse(design.matrix[-which(is.na(data.vector)),]) %*%
                                data.vector[-which(is.na(data.vector))])

      model.vector[which(is.na(data.vector) == F)] <-
        data.vector[which(is.na(data.vector) == F)]

      expanded.df[, data.name] <- model.vector

      ## concatenate data frames
      dataout.df <- rbind(expanded.df, dataout.df)

    }
  }

  return(dataout.df)
}
##########################################################################


##########################################################################
##----------------------AnovaReconstructIterative-----------------------##
##########################################################################

AnovaReconstructIterative <- function(datain.df,
                                      factor1.name,
                                      factor2.name,
                                      data.name,
                                      iterations.num) {

  ## Performs a simple missing value reconstruction with two factors. The
  ## algorithm follows the iterative procedure of Deque et al. (2007).
  ##
  ## Args:
  ##   datain.df: WUX dataframe
  ##   factor1.name: name of the 1st factor
  ##   factor2.name: name of the 2nd factor
  ##   data.name: name of the variable to be reconstructed
  ##   iterations.num: number of iterations to be performed
  ##
  ## Returns:
  ##   dataout.df: WUX dataframe containing the reconstructed values
  ##
  ## Example:
  ##   ## read wux test data
  ##   wux.df <- data(ensembles)
  ##
  ##   ## reconstruction of the  missing data
  ##   reconstruct.df <- AnovaReconstructIterative(wux.df, data.name = "perc.delta.precipitation_amount",
  ##     factor1.name = "rcm", factor2.name = "gcm", iterations.num = 100)
  ##
  ## History:
  ##   2010-10-27 | original code (geh)
  ##   2014-11-20 | all subset commands changed to '[' operation (thm)

  factor1.name <- factor1.name
  factor2.name <- factor2.name
  data.name <- data.name
  datain.df <- gdata::drop.levels(datain.df)
  season.levels <- levels(datain.df$season)
  subreg.levels <- levels(datain.df$subreg)

  dataout.df <- data.frame()

  ## loop over subregions and seasons
  for (ii in subreg.levels){
    for (jj in season.levels){

      ## extract data frame
      extracted.df <- datain.df[datain.df[["subreg"]] == ii &
                                datain.df[["season"]] == jj, ]
      extracted.df <- droplevels(extracted.df)

      ## expand data frame
      expanded.df <- ExpandWuxDataFrame(extracted.df, factor1.name, factor2.name)
      expanded.df$season <- factor(levels(expanded.df$season))
      expanded.df$subreg <- factor(levels(expanded.df$subreg))

      ## create data matrix
      expanded.df.expr <- paste("expanded.df <- expanded.df[order(expanded.df$",
                                factor1.name, "), ]", sep = "")
      eval(parse(text = expanded.df.expr))

      data.matrix.expr <- paste("data.matrix <- matrix(expanded.df$",
                                data.name, ", nrow = length(levels(expanded.df$",
                                factor2.name, ")), ", "ncol = length(levels(expanded.df$",
                                factor1.name, ")))", sep = "")
      eval(parse(text = data.matrix.expr))

      ## reconstruction of the  missing values
      reconstruct.matrix <- AnovaFillValues(data.matrix, iterations.num)
      expanded.df[, data.name] <- as.vector(reconstruct.matrix)

      ## concatenate data frames
      dataout.df <- rbind(expanded.df, dataout.df)

    }
  }

  return(dataout.df)

}

#########################################################################


##########################################################################
##---------------------AnovaReconstructIterativeCC----------------------##
##########################################################################

AnovaReconstructIterativeCC <- function(datain.df,
                                        factor1.name,
                                        factor2.name,
                                        data.name,
                                        iterations.num) {

  ## Performs a leave one out cross calculation (CC) of the reconstructed data
  ## for two factors following the iterative procedure of Deque et al. (2007).
  ##
  ## Args:
  ##   datain.df: WUX dataframe
  ##   factor1.name: name of the 1st factor
  ##   factor2.name: name of the 2nd factor
  ##   data.name: name of the variable to be reconstructed
  ##   iterations.num: number of iterations to be performed
  ##
  ## Returns:
  ##   dataout.df: WUX dataframe with the reconstructed values based on
  ##     the leave one out cross calculation. A new variable in the dataframe
  ##     is created with the variable name "NameOfVariable_crosscal".
  ##
  ## Example:
  ##   ## read wux test data
  ##   wux.df <- data(ensembles)
  ##
  ##   ## cross calculation of the missing data
  ##   reconstruct.df <- AnovaReconstructIterativeCC(wux.df,
  ##     data.name = "perc.delta.precipitation_amount", factor1.name = "rcm",
  ##     factor2.name = "gcm", iterations.num = 100)
  ##
  ## History:
  ##   2010-10-27 | original code (geh)
  ##   2014-11-20 | all subset commands changed to '[' operation (thm)

  data.name <- data.name
  factor1.name <- factor1.name
  factor2.name <- factor2.name
  datain.df <- gdata::drop.levels(datain.df)
  season.levels <- levels(datain.df$season)
  subreg.levels <- levels(datain.df$subreg)

  dataout.df <- data.frame()

  ## loop over subregions and seasons
  for (ii in subreg.levels) {
    for (jj in season.levels) {

      ## extract data frame
      subset.df <- datain.df[datain.df[["subreg"]] == ii &
                             datain.df[["season"]] == jj , ]
      subset.df <- droplevels(subset.df)

      ## expand data frame
      expanded.df <- ExpandWuxDataFrame(subset.df, factor1.name, factor2.name)
      expanded.df$season <- factor(levels(expanded.df$season))
      expanded.df$subreg <- factor(levels(expanded.df$subreg))

      ## create data matrix
      expanded.df.expr <- paste("expanded.df <- expanded.df[order(expanded.df$",
                                factor1.name, "), ]", sep = "")
      eval(parse(text = expanded.df.expr))

      data.matrix.expr <- paste("data.matrix <- matrix(expanded.df$",
                                data.name, ", nrow = length(levels(expanded.df$",
                                factor2.name, ")), ", "ncol = length(levels(expanded.df$",
                                factor1.name, ")))", sep = "")
      eval(parse(text = data.matrix.expr))

      ## leave one out cross calculation
      crossval.matrix <- data.matrix
      values.matrix <- which(is.na(data.matrix) == FALSE, arr.ind = TRUE)

      for (kk in seq(along = values.matrix[ ,1])) {
        print(kk)
        tmp.matrix <- data.matrix
        tmp.matrix[values.matrix[kk, 1], values.matrix[kk, 2]] <- NA

        crossval.matrix[values.matrix[kk, 1], values.matrix[kk, 2]] <-
          AnovaFillValues(tmp.matrix, iterations.num)[values.matrix[kk, 1], values.matrix[kk, 2]]
      }

      crossval.name <- paste(data.name, "_", "crosscal", sep = "")
      expanded.df[crossval.name] <- as.vector(crossval.matrix)

      ## concatenate data frames
      dataout.df <- rbind(expanded.df, dataout.df)

    }
  }

  return(dataout.df)

}
##########################################################################


##########################################################################
##---------------------------AnovaFillValues----------------------------##
##########################################################################

AnovaFillValues <- function(data.matrix,
                            iterations.num) {

  ## Calculates the reconstructed data for two factors following
  ## the iterative procedure of Deque et al. (2007).
  ##
  ## Args:
  ##   data.matrix: matrix of the two factors with missing values
  ##   iterations.num: number of iterations to be performed
  ##
  ## Returns:
  ##   data.matrix: matrix of the two factors with the reconstructed
  ##     missing values
  ##
  ## Example:
  ##   reconstruct.matrix <- AnovaFillValues(data.matrix, 100)
  ##
  ## History:
  ##   2010-10-27 | original code (geh)

                                        # missing value positions
  NA.pos.matrix <- which(is.na(data.matrix), arr.ind = TRUE)

                                        # loop over indicated iterations and missing positions
  for (ii in seq(along = c(1:iterations.num))) {
    for (jj in seq(along = NA.pos.matrix[ ,1])) {
      data.matrix[NA.pos.matrix[jj,1],NA.pos.matrix[jj,2]] <-
        mean(data.matrix[NA.pos.matrix[jj,1], ], na.rm = TRUE) +
          mean(data.matrix[ ,NA.pos.matrix[jj,2]], na.rm = TRUE) -
            mean(data.matrix, na.rm = TRUE)
    }
  }

  return(data.matrix)

}
##########################################################################


##########################################################################
##------------------------------WuxAnova--------------------------------##
##########################################################################
aovWux <- function(model.formula = formula(model.formula),
                     datain.df) {

  ## Calculates an analysis of variance (ANOVA) based on the
  ## specified model
  ##
  ## Args:
  ##   model.formula: model formula used for aov()
  ##   datain.df: WUX data frame
  ##
  ## Returns:
  ##   dataout.list: list containing the ANOVA results for each subregion
  ##     and season. names of list entries are "subreg = xx;season = yy".
  ##
  ## Example:
  ##   ## read WUX test data
  ##   wux.df <- data(ensembles)
  ##
  ##   ## data reconstruction to obtain a balanced design
  ##   reconstruct.df <- AnovaReconstructLES(wux.df,
  ##     factor1.name = "rcm", factor2.name = "gcm", data.name =
  ##     "perc.delta.precipitation_amount")
  ##
  ##   ## calculate ANOVA
  ##   anova.list <- WuxANOVA(perc.delta.precipitation_amount ~ rcm +
  ##     gcm, reconstruct.df)
  ##
  ## History:
  ##   2010-10-27 | original code (geh)
  ##   2014-11-20 | all subset commands changed to '[' operation (thm)

  datain.df <- droplevels(datain.df)
  season.levels <- levels(datain.df$season)
  subreg.levels <- levels(datain.df$subreg)

  dataout.list <- list()

                                        # loop over subregions and seasons
  for (ii in subreg.levels) {
    for (jj in season.levels) {

                                        # extract data frame
      subset.df <- datain.df[datain.df[["subreg"]] == ii &
                             datain.df[["season"]] == jj , ]
      subset.df <- droplevels(subset.df)

                                        # calculate anova and write result to dataout.list
      tmp.list <- aov(model.formula, subset.df)
      tmp.list.name <- paste("subreg = ", ii, ";", "season = ", jj, sep = "")

      dataout.list[[tmp.list.name]] <- tmp.list

    }
  }

  class(dataout.list) <- append("wux.aov", class(dataout.list))
  return(dataout.list)

}

##########################################################################
