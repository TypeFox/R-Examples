context("Pre-process tests")

test_that("Sparce data", {
  df <- data.frame(
    X1 = factor(c("1", "2", rep(NA, 98))),
    X2 = c(1, 2, rep(NA, 98)),
    X3 = factor(c("a", "b", rep(NA, 98))),
    X4 = c("1", "2", rep(NA, 98)),
    X5 = c("a", "b", rep(NA, 98)),
    X6 = c("alfa", "beta", rep(NA, 98))
  )
  df$Label1 <- c(sample(c(0,1), 100, replace = TRUE))
  df$Label2 <- c(sample(c(0,1), 100, replace = TRUE))
  mdata <- mldr::mldr_from_dataframe(df, labelIndices = c(7, 8), name = "testMLDR")

  new.data <- fill_sparce_mldata(mdata)
  expect_equal(as.numeric(new.data$dataset[, 1]),  c(1, 2, rep(0, 98)))
  expect_equal(as.numeric(new.data$dataset[, 2]),  c(1, 2, rep(0, 98)))
  expect_equal(as.character(new.data$dataset[, 3]),  c("a", "b", rep("", 98)))
  expect_equal(as.numeric(new.data$dataset[, 4]),  c(1, 2, rep(0, 98)))
  expect_equal(as.character(new.data$dataset[, 5]),  c("a", "b", rep("", 98)))
  expect_equal(as.character(new.data$dataset[, 6]),
               c("alfa", "beta", rep("", 98)))
  expect_equal(new.data$name, mdata$name)
})

test_that("Normalize data", {
  df <- data.frame(
    X1 = seq(1, 100, by=2),
    X2 = rnorm(100),
    X3 = rnorm(100, 1000, 30),
    X4 = sample(c(stats::runif(90, 0, 1000), rep(NA, 10))),
    X5 = stats::runif(100, -50, 700),
    X6 = c("alfa", "beta", rep("gama", 98))
  )
  df$Label1 <- c(sample(c(0,1), 100, replace = TRUE))
  df$Label2 <- c(sample(c(0,1), 100, replace = TRUE))
  mdata <- mldr::mldr_from_dataframe(df, labelIndices = c(7, 8), name = "testMLDR")

  new.data <- normalize_mldata(mdata)
  for (i in seq(5)) {
    new.col <- as.numeric(new.data$dataset[, i])
    expect_equal(max(new.col, na.rm = TRUE),  1)
    expect_equal(min(new.col, na.rm = TRUE),  0)
    expect_equal(which.max(new.col), which.max(df[, i]))
    expect_equal(which.min(new.col), which.min(df[, i]))
  }
  expect_equal(new.data$dataset[, 6], mdata$dataset[, 6])
  expect_equal(new.data$name, mdata$name)
})

test_that("Remove examples and attributes", {
  df <- data.frame(
    X1 = rep(1, 100),
    X2 = rep(c(1,2), 50),
    X3 = stats::runif(100, 1, 3),
    X4 = rep("XYZ", 100),
    X5 = sample(c("abc", "bcd"), 100, replace = TRUE),
    X6 = c("alfa", "beta", rep("gama", 98)),
    X7 = sample(c(rep(1, 90), rep(NA, 10))),
    X8 = sample(c(rnorm(90), rep(NA, 10)))
  )
  df$Label1 <- rep(0, 100)
  df$Label2 <- sample(c(rep(1, 30), rep(0, 30),
                        sample(c(0,1), 40, replace = TRUE)))
  mdata <- mldr::mldr_from_dataframe(df, labelIndices = c(9, 10), name = "testMLDR")

  new.data <- remove_attributes(mdata, 2)
  expect_equal(new.data$measures$num.attributes, 9)
  expect_named(new.data$dataset[new.data$attributesIndexes],
               c("X1", "X3", "X4", "X5", "X6", "X7", "X8"))
  new.data <- remove_attributes(new.data, c("X3","X6","X8"))
  expect_equal(new.data$measures$num.attributes, 6)
  expect_named(new.data$dataset[new.data$attributesIndexes],
               c("X1", "X4", "X5", "X7"))
  expect_equal(new.data$labels[, c("count","freq")],
               mdata$labels[, c("count","freq")])
  same.data <- remove_attributes(new.data, c("Label1", "ABC"))
  expect_equal(same.data$dataset, new.data$dataset)
  same.data <- remove_attributes(new.data, c(5,7,10))
  expect_equal(same.data$dataset, new.data$dataset)
  expect_equal(new.data$name, mdata$name)
  expect_equal(same.data$name, mdata$name)

  new.data <- remove_unique_attributes(mdata)
  expect_equal(new.data$measures$num.attributes, 8)
  expect_named(new.data$dataset[new.data$attributesIndexes],
               c("X2", "X3", "X5", "X6", "X7", "X8"))
  expect_equal(new.data$name, mdata$name)

  new.data <- remove_unlabeled_instances(mdata)
  has.label <- mdata$dataset$Label2 == 1
  expect_equal(new.data$measures$num.instances, sum(has.label))
  expect_equal(new.data$dataset[mdata$attributesIndexes],
               mdata$dataset[has.label, mdata$attributesIndexes])
  expect_equal(new.data$name, mdata$name)

  df$Label3 <- c(c(1, 1), rep(0, 98))
  df$Label4 <- c(c(0, 0), rep(1, 98))
  df$Label5 <- rep(1, 100)
  df$Label6 <- c(rep(1, 11), rep(0, 89))
  mdata <- mldr::mldr_from_dataframe(df, labelIndices = 9:14, name = "testMLDR")

  new.data <- remove_labels(mdata, 9)
  expect_equal(new.data$measures$num.labels, 5)
  expect_equal(rownames(new.data$labels),
               c("Label2", "Label3", "Label4", "Label5", "Label6"))
  new.data <- remove_labels(new.data, c("Label3","Label5","Label6"))
  expect_equal(new.data$measures$num.labels, 2)
  expect_equal(rownames(new.data$labels), c("Label2", "Label4"))
  expect_equal(new.data$dataset[new.data$attributesIndexes],
               mdata$dataset[mdata$attributesIndexes])
  same.data <- remove_labels(new.data, c("X1", "ABC"))
  expect_equal(same.data$dataset, new.data$dataset)
  same.data <- remove_labels(new.data, c(2,12))
  expect_equal(same.data$dataset, new.data$dataset)
  expect_equal(new.data$name, mdata$name)
  expect_equal(same.data$name, mdata$name)

  new.data <- remove_skewness_labels(mdata)
  expect_equal(new.data$measures$num.labels, 4)
  expect_equal(rownames(new.data$labels),
               c("Label2", "Label3", "Label4", "Label6"))
  expect_equal(new.data$name, mdata$name)

  new.data <- remove_skewness_labels(mdata, 2)
  expect_equal(new.data$measures$num.labels, 2)
  expect_equal(rownames(new.data$labels), c("Label2", "Label6"))
  expect_equal(new.data$name, mdata$name)

  new.data <- remove_skewness_labels(mdata, 10)
  expect_equal(new.data$measures$num.labels, 2)
  expect_equal(rownames(new.data$labels), c("Label2", "Label6"))
  expect_equal(new.data$name, mdata$name)

  expect_error(remove_skewness_labels(mdata, 11))
})

test_that("Replace nominal attributes", {
  df <- data.frame(
    X1 = sample(c("abc", "bcd", "cde"), 100, replace = TRUE),
    X2 = c(1, 2, rep(NA, 98)),
    X3 = factor(c(rep("a", 10), rep(NA, 90))),
    X4 = c("1", "2", rep(NA, 98)),
    X5 = c("alfa", "beta", rep(NA, 98))
  )
  zero.um <- c(rep(1, 30), rep(0, 30))
  df$Label1 <- sample(c(zero.um, sample(c(0,1), 40, replace = TRUE)))
  df$Label2 <- sample(c(zero.um, sample(c(0,1), 40, replace = TRUE)))
  mdata <- mldr::mldr_from_dataframe(df, labelIndices = c(6, 7), name = "testMLDR")

  new.data <- replace_nominal_attributes(mdata)
  expect_equal(new.data$measures$num.attributes, 8)
  expect_equal(colnames(new.data$dataset[,new.data$attributesIndexes]),
               c("X1_abc", "X1_bcd", "X2", "X3_a", "X4_1", "X5_alfa"))
  expect_equal(new.data$dataset[,"X1_abc"], as.numeric(df$X1 == "abc"))
  expect_equal(new.data$dataset[,"X1_bcd"], as.numeric(df$X1 == "bcd"))
  expect_equal(new.data$name, mdata$name)
})

test_that("Alternatives datasets", {
  df <- data.frame(
    Label1 = c(sample(c(0,1), 100, replace = TRUE)),
    Label2 = c(sample(c(0,1), 100, replace = TRUE)),
    Label3 = c(sample(c(0,1), 100, replace = TRUE)),
    X1 = factor(c("1", "2", rep(NA, 98))),
    X2 = c(1, 2, rep(NA, 98)),
    X3 = factor(c("a", "b", rep(NA, 98))),
    X4 = c("1", "2", rep(NA, 98)),
    X5 = c("a", "b", rep(NA, 98)),
    X6 = c("alfa", "beta", rep(NA, 98))
  )
  df[df$Label1 == 0 & df$Label2 == 0,"Label3"] <- 1
  mdata <- mldr::mldr_from_dataframe(df, labelIndices = c(1, 2, 3), name = "testMLDR")

  ndata <- fill_sparce_mldata(mdata)
  expect_equal(ndata$measures, mdata$measures)
  expect_equal(ndata$labels, mdata$labels)
  expect_equal(ndata$name, mdata$name)

  new.data <- remove_labels(mdata, "Label2")
  expect_equal(new.data$measures$num.labels, 2)
  expect_equal(new.data$labels$index, c(1,2))
  expect_equal(rownames(new.data$labels), c("Label1","Label3"))
  expect_equal(new.data$name, mdata$name)

  new.data <- remove_attributes(mdata, c("X3","X6","X8"))
  expect_equal(new.data$measures$num.attributes, 7)
  expect_equal(new.data$labels[c("index","count","freq")],
               mdata$labels[c("index","count","freq")])
  expect_equal(new.data$name, mdata$name)

  ndata <- remove_unique_attributes(ndata)
  expect_equal(ndata$measures, mdata$measures)
  expect_equal(ndata$labels, mdata$labels)
  expect_equal(ndata$name, mdata$name)

  ndata <- remove_skewness_labels(ndata)
  expect_equal(ndata$measures, mdata$measures)
  expect_equal(ndata$labels, mdata$labels)
  expect_equal(ndata$name, mdata$name)

  ndata <- remove_unlabeled_instances(ndata)
  expect_equal(ndata$measures, mdata$measures)
  expect_equal(ndata$labels, mdata$labels)
  expect_equal(ndata$name, mdata$name)

  ndata <- normalize_mldata(ndata)
  expect_equal(ndata$measures, mdata$measures)
  expect_equal(ndata$labels, mdata$labels)
  expect_equal(ndata$name, mdata$name)

  ndata <- replace_nominal_attributes(ndata)
  attrs <- c("num.instances", "num.labels",
             "num.labelsets", "num.single.labelsets",
             "max.frequency", "cardinality", "density")
  expect_equal(ndata$measures[attrs], mdata$measures[attrs])
  expect_equal(ndata$labels, mdata$labels)
  expect_equal(ndata$name, mdata$name)
})
