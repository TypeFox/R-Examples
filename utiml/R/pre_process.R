#' Fill sparce dataset with 0 or '' values
#'
#' Transform a sparce dataset filling NA values to 0 or '' based on the column
#' type. Text columns with numeric values will be modified to numerical.
#'
#' @family pre process
#' @param mdata The mldr dataset to be filled.
#' @return a new mldr object.
#' @export
#'
#' @examples
#' sparce.toy <- toyml
#' sparce.toy$dataset$ratt10[sample(100, 30)] <- NA
#' complete.toy <- fill_sparce_mldata(sparce.toy)
fill_sparce_mldata <- function(mdata) {
  is.letter <- function(x) {
    grepl("[[:alpha:]]", x)
  }

  attrs <- seq(mdata$measures$num.attributes)
  new.cols <- lapply(mdata$dataset[, attrs], function(col) {
    if (anyNA(col)) {
      # Has NA value
      if (is.numeric(col)) {
        # Numeric value - fill with 0
        col[is.na(col)] <- 0
      }
      else if (any(is.letter(col))) {
        # Text value - fill with ''
        col <- as.character(col)
        col[is.na(col)] <- ""
      }
      else {
        # Text but with numeric values - convert to numeric and fill with 0
        col <- as.numeric(as.character(col))
        col[is.na(col)] <- 0
      }
    }
    col
  })

  dataset <- data.frame(row.names = rownames(mdata$dataset))
  dataset <- cbind(dataset, new.cols)
  mldr::mldr_from_dataframe(dataset, mdata$labels$index, name = mdata$name)
}

#' Normalize numerical attributes
#'
#' Normalize all numerical attributes to values between 0 and 1. The highest
#' value is changed to 1 and the lowest value to 0.
#'
#' @family pre process
#' @param mdata The mldr dataset to be normalized.
#' @return a new mldr object.
#' @export
#'
#' @examples
#' norm.toy <- normalize_mldata(toyml)
normalize_mldata <- function(mdata) {
  data <- mdata$dataset[seq(mdata$measures$num.attributes)]
  for (col in mdata$attributesIndexes) {
    if (is.numeric(data[, col])) {
      data[col] <- utiml_normalize(data[col])
    }
  }
  mldr::mldr_from_dataframe(data, mdata$labels$index, name = mdata$name)
}

#' Remove attributes from the dataset
#'
#' Remove spectified attributes generating a new multi-label dataset.
#'
#' @family pre process
#' @param mdata The mldr dataset to remove labels.
#' @param attributes Attributes indexes or attributes names to be removed.
#' @return a new mldr object.
#' @note If invalid attributes names or indexes were informed, they will be
#'  ignored.
#' @export
#'
#' @examples
#' toyml1 <- remove_attributes(toyml, c("iatt8","iatt9", "ratt10"))
#' toyml2 <- remove_attributes(toyml, 10)
remove_attributes <- function (mdata, attributes) {
  if (mode(attributes) == "character") {
    attributes <- which(colnames(mdata$dataset) %in% attributes)
  }

  use.attributes <- setdiff(seq(mdata$measures$num.attributes), attributes)
  create_subset(mdata, seq(mdata$measures$num.instances), use.attributes)
}

#' Remove labels from the dataset
#'
#' Remove spectified labels generating a new multi-label dataset.
#'
#' @family pre process
#' @param mdata The mldr dataset to remove labels.
#' @param labels Label indexes or label names to be removed.
#' @return a new mldr object.
#' @note If invalid labels names or indexes were informed, they will be ignored.
#' @export
#'
#' @examples
#' toyml1 <- remove_labels(toyml, c("y1","y5"))
#' toyml2 <- remove_labels(toyml, c(11, 15))
remove_labels <- function (mdata, labels) {
  if (mode(labels) == "character") {
    labels <- mdata$labels[labels, "index"]
    labels <- labels[!is.na(labels)]
  }
  else {
    # Only labels index, not attributes index
    labels <- mdata$labels$index[which(mdata$labels$index %in% labels)]
  }

  new.attrs <- setdiff(seq(mdata$measures$num.attributes), labels)
  dataset <- mdata$dataset[new.attrs]
  labels <- which(colnames(dataset) %in% rownames(mdata$labels))

  if (length(labels) <= 1) {
    stop("The pre process procedure result in a single label")
  }

  mldr::mldr_from_dataframe(dataset, labels, name = mdata$name)
}

#' Remove unique attributes
#'
#' Remove the attributes that have a single value for all instances. Empty and
#' NA values are considered different values.
#'
#' @family pre process
#' @param mdata The mldr dataset to remove.
#' @return a new mldr object.
#' @export
#'
#' @examples
#' alt.toy <- toyml
#' alt.toy$dataset$ratt10 <- mean(alt.toy$dataset$ratt10)
#' new.toy <- remove_unique_attributes(alt.toy)
remove_unique_attributes <- function(mdata) {
  labelsIndexes <- c()
  attributesIndexes <- c()

  for (col in seq(mdata$measures$num.attributes)) {
    if (col %in% mdata$labels$index) {
      attributesIndexes <- c(attributesIndexes, col)
      labelsIndexes <- c(labelsIndexes, length(attributesIndexes))
    }
    else {
      if (length(unique(mdata$dataset[, col])) > 1) {
        attributesIndexes <- c(attributesIndexes, col)
      }
    }
  }

  mldr::mldr_from_dataframe(mdata$dataset[attributesIndexes],
                            labelsIndexes,
                            name = mdata$name)
}

#' Remove examples without labels
#'
#' Remove the examples that do not have labels.
#'
#' @family pre process
#' @param mdata The mldr dataset to remove the instances.
#' @return a new mldr object.
#' @export
#'
#' @examples
#' new.toy <- remove_labels(toyml, c(12,14))
#' remove_unlabeled_instances(new.toy)
remove_unlabeled_instances <- function(mdata) {
  labelset <- rep(0, mdata$measures$num.labels)
  rows <- !apply(mdata$dataset[mdata$labels$index] == labelset, 1, all)
  cols <- seq(mdata$measures$num.attributes)

  mldr::mldr_from_dataframe(mdata$dataset[rows, cols],
                            mdata$labels$index,
                            name = mdata$name)
}

#' Remove unusual or very common labels
#'
#' Remove the labels that have smaller number of positive or negative examples
#' based on a specific threshold value.
#'
#' @family pre process
#' @param mdata The mldr dataset to remove the skewness labels.
#' @param t Threshold value. Number of minimum examples positive and negative.
#' @return a new mldr object.
#' @export
#'
#' @examples
#' remove_skewness_labels(toyml, 20)
remove_skewness_labels <- function(mdata, t = 1) {
  labelsIndexes <- c()

  for (col in mdata$labels$index) {
    tbl <- table(mdata$dataset[col])
    if (length(tbl) > 1 && all(tbl > t)) {
      labelsIndexes <- c(labelsIndexes, col)
    }
  }

  if (length(labelsIndexes) <= 1) {
    stop("The pre process procedure result in a single label")
  }

  dataset <- mdata$dataset[sort(c(mdata$attributesIndexes, labelsIndexes))]
  labels <- which(colnames(dataset) %in% rownames(mdata$labels))

  mldr::mldr_from_dataframe(dataset, labels, name = mdata$name)
}

#' Replace nominal attributes
#' Replace the nominal attributes by binary attributes.
#'
#' @family pre process
#' @param mdata The mldr dataset to remove.
#' @param ordinal.attributes Not yet, but it will be used to specify which
#'  attributes need to be replaced.
#' @return a new mldr object.
#' @export
#'
#' @examples
#' new.toy <- toyml
#' new.column <- as.factor(sample(c("a","b","c"), 100, replace = TRUE))
#' new.toy$dataset$ratt10 <- new.column
#' head(replace_nominal_attributes(new.toy))
replace_nominal_attributes <- function(mdata, ordinal.attributes = list()) {
  # TODO ordinal.attributes
  replace_nominal_column <- function(column, column.name = "", type = 1) {
    column <- as.factor(column)
    symbols <- levels(column)
    result <- {}

    if (length(symbols) == 2 && type == 1 && 0 %in% symbols && 1 %in% symbols) {
      result <- cbind(result, as.double(column == 1))
      names <- column.name
    }
    else {
      for (i in seq(length(symbols) - type)) {
        result <- cbind(result, as.double(column == symbols[i]))
      }
      names <- paste(column.name, symbols[seq(length(symbols) - type)], sep="_")
    }

    if (column.name != "") {
      colnames(result) <- names
    }

    result
  }

  dataset <- data.frame(row.names = rownames(mdata$dataset))
  labelIndexes <- c()
  for (col in seq(mdata$measures$num.attributes)) {
    if (is.numeric(mdata$dataset[, col])) {
      dataset <- cbind(dataset, mdata$dataset[col])
      if (col %in% mdata$labels$index) {
        labelIndexes <- c(labelIndexes, ncol(dataset))
      }
    }
    else {
      column <- replace_nominal_column(mdata$dataset[, col],
                                       colnames(mdata$dataset[col]))
      dataset <- cbind(dataset, column)
    }
  }

  mldr::mldr_from_dataframe(dataset, labelIndexes, name = mdata$name)
}
