#
# Contains functions for calculating measures
# for a multilabel dataset
#
measures <- function(mld) {
  if(nrow(mld$dataset) > 0) {
    labelsets <- do.call(paste, c(mld$dataset[, mld$labels$index], sep = ""))
    labelsets <- table(as.factor(labelsets))

    list(
      num.attributes = length(names(mld$dataset)) - 2, # 2 columns are measures!
      num.instances = nrow(mld$dataset),
      num.inputs = length(names(mld$dataset)) - 2 - nrow(mld$labels),
      num.labels = nrow(mld$labels),

      # Number of different labelsets
      num.labelsets = length(labelsets),

      # Number of labelsets appearing only once
      num.single.labelsets = sum(labelsets == 1),


      # Maximum frequency on a labelset
      max.frequency = max(labelsets),

      cardinality = mean(mld$dataset$.labelcount),
      density = mean(mld$dataset$.labelcount) / nrow(mld$labels),
      meanIR = mean(mld$labels$IRLbl, na.rm = TRUE),  # Avoid NA IRLbls
      scumble = mean(mld$dataset$.SCUMBLE),
      scumble.cv = sd(mld$dataset$.SCUMBLE)/mean(mld$dataset$.SCUMBLE),
      tcs = log((length(names(mld$dataset)) - 2 - nrow(mld$labels)) * nrow(mld$labels) * length(labelsets))
    )
  } else
    list(
      num.attributes = length(names(mld$dataset)) - 2,
      num.instances = nrow(mld$dataset),
      num.inputs = length(names(mld$dataset)) - 2 - nrow(mld$labels),
      num.labels = nrow(mld$labels),
      num.labelsets = NA,
      num.single.labelsets = NA,
      max.frequency = NA,
      cardinality = NA,
      density = NA,
      meanIR = NA,
      scumble = NA,
      scumble.cv = NA,
      tcs = 0
    )
}

label_measures <- function(data, indexes) {
  if(nrow(data) > 0) {
    msr <- data.frame(
      index = indexes,
      count = colSums(data[indexes])
    )

    msr$freq <- msr$count / nrow(data)
    msr$IRLbl <- max(msr$count) / msr$count
    msr$IRLbl[msr$IRLbl == Inf] <- 0 # Avoid Inf values for labels with 0 samples

    msr
  } else
    data.frame(index = indexes, count = 0, freq = NA, IRLbl = NA)
}

dataset_measures <- function(mld) {
  if(nrow(mld$dataset) > 0) {
    mld$dataset$.labelcount <- rowSums(mld$dataset[, mld$labels$index])

    # Atkinson index (aka SCUMBLE) for IR in each instance
    IRs <- data.frame(t(t(mld$dataset[, mld$labels$index]) * mld$labels$IRLbl))
    IRmeans <- rowSums(IRs) / mld$dataset$.labelcount

    IRs[IRs == 0] <- 1            # Identity element for (R, *)
    IRprod <- Reduce("*", IRs)    # Row products
    mld$dataset$.SCUMBLE <- ifelse(mld$dataset$.labelcount > 0,
                                   1 - (IRprod)^(1/mld$dataset$.labelcount) / IRmeans,
                                   0)

    # lblSCUMBLE: SCUMBLE mean by label - Avoid dividing by 0 if there are not appearances of a label
    mld$labels$SCUMBLE <- ifelse(colSums(mld$dataset[mld$labels$index]) == 0,
                                 0,
                                 sapply(mld$labels$index, function(i) mean(mld$dataset[mld$dataset[, i] != 0, ".SCUMBLE"])))

    # lblSCUMBLE.CV: Coefficient of variation of the corresponding SCUMBLE mean
    ####################################
    # If the mean SCUMBLE of a label is 0, then its standard deviation is 0 (since SCUMBLE is always nonnegative),
    # thus we define its CV as 0.
    ####################################
    mld$labels$SCUMBLE.CV <- ifelse(colSums(mld$dataset[mld$labels$index]) <= 1,
      NA,
      ifelse(mld$labels$SCUMBLE == 0,
        0,
        sapply(mld$labels$index, function(i) sd(mld$dataset[mld$dataset[, i] != 0, ".SCUMBLE"])) / mld$labels$SCUMBLE
      )
    )
  }
  else {
    mld$dataset$.labelcount <- numeric()
    mld$dataset$.SCUMBLE <- numeric()
    mld$labels$SCUMBLE <- numeric()
    mld$labels$SCUMBLE.CV <- numeric()
  }

  mld
}
