#' Create a holdout partition based on the specified algorithm
#'
#' This method creates multi-label dataset for train, test, validation or other
#' proposes the partition method defined in \code{method}. The number of
#' partitions is defined in \code{partitions} parameter. Each instance is used
#' in only one partition of divistion.
#'
#' @family sampling
#' @param mdata A mldr dataset.
#' @param partitions A list of percentages or a single value. The sum of all
#'  values does not be greater than 1. If a single value is informed then the
#'  complement of them is applied to generated the second partition. If two or
#'  more values are informed and the sum of them is lower than 1 the partitions
#'  will be generated with the informed proportion. If partitions have names,
#'  they are used to name the return. (Default: \code{c(train=0.7, test=0.3)}).
#' @param method The method to split the data. The default methods are:
#'  \describe{
#'    \item{random}{Split randomly the folds.}
#'    \item{iterative}{Split the folds considering the labels proportions
#'                      individually. Some specific label can not occurs in all
#'                      folds.}
#'    \item{stratified}{Split the folds considering the labelset proportions.}
#'  }
#'  You can also create your own partition method. See the note and example
#'  sections to more details. (Default: "random")
#' @return A list with at least two datasets sampled as specified in partitions
#'  parameter.
#' @references Sechidis, K., Tsoumakas, G., & Vlahavas, I. (2011). On the
#'  stratification of multi-label data. In Proceedings of the Machine
#'  Learning and Knowledge Discovery in Databases - European Conference,
#'  ECML PKDD (pp. 145-158).
#' @note To create your own split method, you need to build a function that
#'  receive a mldr object and a list with the proportions of examples in each
#'  fold and return an other list with the index of the elements for each fold.
#' @export
#'
#' @examples
#' dataset <- create_holdout_partition(toyml)
#' names(dataset)
#' ## [1] "train" "test"
#' #dataset$train
#' #dataset$test
#'
#' dataset <- create_holdout_partition(toyml, c(a=0.1, b=0.2, c=0.3, d=0.4))
#' #' names(dataset)
#' #' ## [1] "a" "b" "c" "d"
#'
#' sequencial_split <- function (mdata, r) {
#'  S <- list()
#'
#'  amount <- trunc(r * mdata$measures$num.instances)
#'  indexes <- c(0, cumsum(amount))
#'  indexes[length(r)+1] <- mdata$measures$num.instances
#'
#'  S <- lapply(seq(length(r)), function (i) {
#'    seq(indexes[i]+1, indexes[i+1])
#'  })
#'
#'  S
#' }
#' dataset <- create_holdout_partition(toyml, method="sequencial_split")
create_holdout_partition <- function (mdata,
                                      partitions = c(train=0.7, test=0.3),
                                      method = c("random", "iterative",
                                                 "stratified")) {
  # Validations
  if (class(mdata) != "mldr") {
    stop("First argument must be an mldr object")
  }

  if (sum(partitions) > 1) {
    stop("The sum of partitions can not be greater than 1")
  }

  holdout.method <- utiml_validate_splitmethod(method[1])

  partitions <- utiml_ifelse(length(partitions) == 1,
                             c(partitions, 1 - partitions),
                             partitions)

  # Split data
  folds <- do.call(holdout.method, list(mdata = mdata, r = partitions))
  names(folds) <- names(partitions)
  ldata <- lapply(folds, function (fold) {
    create_subset(mdata, fold, mdata$attributesIndexes)
  })

  ldata
}

#' Create the k-folds partition based on the specified algorithm
#'
#' This method create the kFoldPartition object, from it is possible create
#' the dataset partitions to train, test and optionally to validation.
#'
#' @family sampling
#' @param mdata A mldr dataset.
#' @param k The number of desirable folds. (Default: 10)
#' @param method The method to split the data. The default methods are:
#'  \describe{
#'    \item{random}{Split randomly the folds.}
#'    \item{iterative}{Split the folds considering the labels proportions
#'                      individually. Some specific label can not occurs in all
#'                      folds.}
#'    \item{stratified}{Split the folds considering the labelset
#'                          proportions.}
#'  }
#'  You can also create your own partition method. See the note and example
#'  sections to more details. (Default: "random")
#' @return An object of type kFoldPartition.
#' @references Sechidis, K., Tsoumakas, G., & Vlahavas, I. (2011). On the
#'  stratification of multi-label data. In Proceedings of the Machine
#'  Learning and Knowledge Discovery in Databases - European Conference,
#'  ECML PKDD (pp. 145-158).
#' @note To create your own split method, you need to build a function that
#'  receive a mldr object and a list with the proportions of examples in each
#'  fold and return an other list with the index of the elements for each fold.
#' @seealso \link[=partition_fold]{How to create the datasets from folds}
#' @export
#'
#' @examples
#' k10 <- create_kfold_partition(toyml, 10)
#' k5 <- create_kfold_partition(toyml, 5, "stratified")
#'
#' sequencial_split <- function (mdata, r) {
#'  S <- list()
#'
#'  amount <- trunc(r * mdata$measures$num.instances)
#'  indexes <- c(0, cumsum(amount))
#'  indexes[length(r)+1] <- mdata$measures$num.instances
#'
#'  S <- lapply(seq(length(r)), function (i) {
#'    seq(indexes[i]+1, indexes[i+1])
#'  })
#'
#'  S
#' }
#' k3 <- create_kfold_partition(toyml, 3, "sequencial_split")
create_kfold_partition <- function (mdata,
                                    k = 10,
                                    method = c("random", "iterative",
                                               "stratified")) {
  if (class(mdata) != "mldr") {
    stop("First argument must be an mldr object")
  }

  if (k < 1 || k > mdata$measures$num.instances) {
    stop("The k value is not valid")
  }

  kfold.method <- utiml_validate_splitmethod(method[1])

  kf <- list(dataset = mdata, k = k)
  kf$fold <- do.call(kfold.method, list(mdata = mdata, r = rep(1/k, k)))
  class(kf) <- "kFoldPartition"

  kf
}

#' Create a random subset of a dataset
#'
#' @family sampling
#' @param mdata A mldr dataset
#' @param instances The number of expected instances
#' @param attributes The number of expected attributes.
#'  (Default: all attributes)
#' @return A new mldr subset
#' @export
#'
#' @examples
#' small.toy <- create_random_subset(toyml, 10, 3)
#' medium.toy <- create_random_subset(toyml, 50, 5)
create_random_subset <- function(mdata, instances,
                                 attributes = mdata$measures$num.inputs) {
  if (instances > mdata$measures$num.instances) {
    stop(paste("The expected number of instances is greater than ",
               mdata$measures$num.instances))
  }
  if (attributes > mdata$measures$num.inputs) {
    stop(paste("The expected number of attributes is greater than ",
               mdata$measures$num.inputs))
  }
  rows <- sample(mdata$measures$num.instances, instances)
  cols <- sample(mdata$attributesIndexes, attributes)
  create_subset(mdata, rows, cols)
}

#' Create a subset of a dataset
#'
#' @family sampling
#' @param mdata A mldr dataset
#' @param rows A vector with the instances indexes (names or indexes).
#' @param cols A vector with the attributes indexes (names or indexes).
#' @return A new mldr subset
#' @note It is not necessary specify the labels attributes because they are
#'  included by default.
#' @export
#'
#' @examples
#' ## Create a dataset with the 20 first examples and the 7 first attributes
#' small.toy <- create_subset(toyml, seq(20), seq(7))
#'
#' ## Create a random dataset with 50 examples and 5 attributes
#' random.toy <- create_subset(toyml, sample(100, 50), sample(10, 5))
create_subset <- function(mdata, rows, cols = NULL) {
  if (mode(cols) == "character") {
    cols <- which(colnames(mdata$dataset[mdata$attributesIndexes]) %in% cols)
  }
  else if (is.null(cols)) {
    cols <- mdata$attributesIndexes
  }
  else {
    cols <- intersect(cols, seq(mdata$measures$num.attributes))
  }

  if (mode(rows) == "character") {
    rows <- intersect(rows, rownames(mdata$dataset))
  }
  else {
    rows <- intersect(rows, seq(mdata$measures$num.instances))
  }

  dataset <- mdata$dataset[rows, sort(unique(c(cols, mdata$labels$index)))]
  labelIndexes <- which(colnames(dataset) %in% rownames(mdata$labels))

  mldr::mldr_from_dataframe(dataset, labelIndices=labelIndexes, name=mdata$name)
}

#' Create the multi-label dataset from folds
#'
#' This is a simple way to use k-fold cross validation.
#'
#' @param kfold A \code{kFoldPartition} object obtained from use of the method
#'  \link{create_kfold_partition}.
#' @param n The number of the fold to separated train and test subsets.
#' @param has.validation Logical value that indicate if a validation
#'  dataset will be used. (Defaul: \code{FALSE})
#' @return A list contained train and test mldr dataset:
#'  \describe{
#'    \code{train}{The mldr dataset with train examples, that inclue all
#'      examples except those that are in test and validation samples}
#'    \code{test}{The mldr dataset with test examples, defined by the
#'      number of the fold}
#'    \code{validation}{Optionally, only if \code{has.validation = TRUE}.
#'      The mldr dataset with validation examples}
#'  }
#' @export
#'
#' @examples
#' folds <- create_kfold_partition(toyml, 10)
#'
#' # Using the first partition
#' dataset <- partition_fold(folds, 1)
#' names(dataset)
#' ## [1] "train" "test"
#'
#' # All iterations
#' for (i in 1:10) {
#'    dataset <- partition_fold(folds, i)
#'    #dataset$train
#'    #dataset$test
#' }
#'
#' # Using 3 folds validation
#' dataset <- partition_fold(folds, 3, TRUE)
#' # dataset$train, dataset$test, #dataset$validation
partition_fold <- function(kfold, n, has.validation = FALSE) {
    if (class(kfold) != "kFoldPartition") {
        stop("Second argument must be an 'kFoldPartition' object")
    }

    if (n < 1 || n > kfold$k) {
        stop(cat("The 'n' value must be between 1 and", kfold$k))
    }

    folds <- kfold$fold[-n]
    if (has.validation) {
        i <- n == length(folds)
        v <- c(1, n)[c(i, !i)]
        folds <- folds[-v]
    }

    mdata <- kfold$dataset
    ldata <- list(
      train = create_subset(mdata, unlist(folds), mdata$attributesIndexes),
      test  = create_subset(mdata, kfold$fold[[n]], mdata$attributesIndexes)
    )

    if (has.validation) {
      ldata$validation <- create_subset(mdata, kfold$fold[[v]],
                                        mdata$attributesIndexes)
    }

    ldata
}

# Internal methods -------------------------------------------------------------

#' Return the name of split method and validate if it is valid
#'
#' @param method The method name
#' @return The correct name of split method
utiml_validate_splitmethod <- function (method) {
  DEFAULT.METHODS <- c("random", "iterative", "stratified")
  method.name <- ifelse(method %in% DEFAULT.METHODS,
                       paste("utiml", method, "split", sep = "_"),
                       method)

  if (!exists(method.name, mode = "function")) {
    stop(paste("The partition method '", method.name,
               "' is not a valid function", sep=''))
  }

  method.name
}

#' Internal Iterative Stratification
#'
#' Create the indexes using the Iterative Stratification algorithm.
#'
#' @param mdata A mldr dataset.
#' @param r Desired proportion of examples in each subset r1, . . . rk.
#' @return A list with k disjoint indexes subsets S1, . . .Sk.
#' @references Sechidis, K., Tsoumakas, G., & Vlahavas, I. (2011). On the
#'  stratification of multi-label data. In Proceedings of the Machine
#'  Learningand Knowledge Discovery in Databases - European Conference,
#'  ECML PKDD (pp. 145-158).
#'
#' @examples
#' \dontrun{
#' # Create 3 partitions for train, validation and test
#' indexes <- utiml_iterative_split(emotions, c(0.6,0.1,0.3))
#'
#' # Create a stratified 10-fold
#' indexes <- utiml_iterative_split(emotions, rep(0.1,10))
#' }
utiml_iterative_split <- function(mdata, r) {
  D <- rownames(mdata$dataset)
  S <- lapply(seq(length(r)), function(i) character())

  # Calculate the desired number of examples at each subset
  cj <- round(mdata$measures$num.instances * r)
  dif <- mdata$measures$num.instances - sum(cj)
  if (dif != 0) {
    cj[seq(abs(dif))] <- cj[seq(abs(dif))] + utiml_ifelse(dif > 0, 1, -1)
  }

  # Calculate the desired number of examples of each label at each subset
  cji <- trunc(sapply(mdata$labels$count, function(di) di * r))
  colnames(cji) <- rownames(mdata$labels)

  # Empty examples (without any labels)
  empty.inst <- apply(mdata$dataset[, mdata$labels$index], 1, sum)
  empty.inst <- as.character(which(empty.inst == 0))
  if (length(empty.inst) > 0) {
    D <- setdiff(D, empty.inst)
    prop.empty <- ceiling(length(empty.inst) / length(r))
    indexes <- rep(seq(length(r)), prop.empty)[seq(length(empty.inst))]
    Dist <- split(empty.inst, indexes)
    for (i in 1:length(Dist)) {
      S[[i]] <- Dist[[i]]
      cj[i] <- cj[i] - length(S[[i]])
    }
  }

  while (length(D) > 0) {
    # Find the label with the fewest (but at least one) remaining examples,
    # Do not use apply because sometimes its returns is a matrix
    Dl <- lapply(mdata$labels$index, function(col) {
      D[which(mdata$dataset[D, col] == 1)]
    })
    names(Dl) <- rownames(mdata$labels)
    Di <- unlist(lapply(Dl, length))
    l <- names(which.min(Di[Di > 0]))

    for (ex in Dl[[l]]) {
      # Find the subset(s) with the largest number of desired examples for ,
      # this label, breaking ties by considering the largest number of desired
      # examples
      m <- which(cji[which.max(cji[, l]), l] == cji[, l])
      if (length(m) > 1) {
        m <- intersect(m, which(cj[m[which.max(cj[m])]] == cj))
        if (length(m) > 1) m <- sample(m)[1]
      }

      S[[m]] <- c(S[[m]], ex)
      D <- D[D != ex]

      # Update desired number of examples
      i <- which(mdata$dataset[ex, mdata$labels$index] == 1)
      cji[m, i] <- cji[m, i] - 1
      cj[m] <- cj[m] - 1
    }
  }

  S <- lapply(S, function(fold) {
    new.fold <- which(rownames(mdata$dataset) %in% fold)
    names(new.fold) <- rownames(mdata$dataset[new.fold, ])
    new.fold
  })

  S
}

#' Random split of a dataset
#'
#' @param mdata A mldr dataset.
#' @param r Desired proportion of examples in each subset r1, . . . rk.
#' @return A list with k disjoint indexes subsets S1, . . .Sk.
#'
#' @examples
#' \dontrun{
#' utiml_random_split(emotions, c(0.6, 0.2, 0.2))
#' }
utiml_random_split <- function(mdata, r) {
  index <- c()
  amount <- round(mdata$measures$num.instances * r)

  dif <- mdata$measures$num.instances - sum(amount)
  for (i in seq(abs(dif))) {
    amount[i] <- amount[i] + sign(dif)
  }

  for (i in seq(length(amount))) {
    index <- c(index, rep(i, amount[i]))
  }

  S <- split(sample(seq(mdata$measures$num.instances)), index)
  for (i in 1:length(S)) {
    names(S[[i]]) <- rownames(mdata$dataset[S[[i]], ])
  }

  S
}

#' Labelsets Stratification
#' Create the indexes using the Labelsets Stratification approach.
#'
#' @param mdata A mldr dataset
#' @param r Desired proportion of examples in each subset, r1, . . . rk
#' @return A list with k disjoint indexes subsets S1, . . .Sk
#' @references Sechidis, K., Tsoumakas, G., & Vlahavas, I. (2011). On the
#'  stratification of multi-label data. In Proceedings of the Machine
#'  Learningand Knowledge Discovery in Databases - European Conference,
#'  ECML PKDD (pp. 145-158).
#'
#' @examples
#' \dontrun{
#' # Create 3 partitions for train, validation and test
#' indexes <- utiml_stratified_split(emotions, c(0.6,0.1,0.3))
#'
#' # Create a stratified 10-fold
#' indexes <- utiml_stratified_split(emotions, rep(0.1,10))
#' }
utiml_stratified_split <- function(mdata, r) {
  D <- sample(mdata$measures$num.instances)
  S <- lapply(1:length(r), function(i) integer())
  labelsets <- apply(mdata$dataset[, mdata$labels$index], 1, paste, collapse="")

  # Calculate the desired number of examples of each labelset at each subset
  cji.aux <- sapply(mdata$labelsets, function(di) di * r)
  cji <- trunc(cji.aux)
  dif <- cji.aux - cji
  rest <- round(apply(dif, 1, sum))
  for (ls in rev(names(mdata$labelsets))) {
    s <- sum(dif[, ls])
    if (s > 0) {
      for (i in seq(s)) {
        fold <- which.max(rest)
        rest[fold] <- rest[fold] - 1
        cji[fold, ls] <- cji[fold, ls] + 1
      }
    }
  }

  for (ex in D) {
    ls <- labelsets[ex]
    fold <- which.max(cji[, ls])
    if (cji[fold, ls] > 0) {
      S[[fold]] <- c(S[[fold]], ex)
      cji[fold, ls] <- cji[fold, ls] - 1
    }
  }

  for (i in seq(length(S))) {
    names(S[[i]]) <- rownames(mdata$dataset[S[[i]], ])
  }

  S
}

#' Print a kFoldPartition object
#'
#' @param x The kFoldPartition object
#' @param ... ignored
#' @export
print.kFoldPartition <- function (x, ...) {
  cat("K Fold Partition", paste("(k = ",x$k,")", sep=''), "\n\n")

  folds <- rbind(lapply(x$fold, length))
  rownames(folds) <- c("Examples:")
  colnames(folds) <- paste("Fold", seq(x$k), sep='_')
  print(folds)
  cat("\n")
  for (i in seq(x$k)) {
    cat(paste("Fold ", i, ":", sep=''), x$fold[[i]], "\n")
  }
}
