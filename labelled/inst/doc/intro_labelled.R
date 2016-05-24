## ------------------------------------------------------------------------
library(labelled)

var_label(iris$Sepal.Length) <- "Length of sepal"

## ------------------------------------------------------------------------
var_label(iris) <- list(Petal.Length = "Length of petal", Petal.Width = "Width of Petal")

## ------------------------------------------------------------------------
var_label(iris$Petal.Width)
var_label(iris)

## ------------------------------------------------------------------------
var_label(iris$Sepal.Length) <- NULL

## ---- eval=FALSE---------------------------------------------------------
#  View(iris)

## ------------------------------------------------------------------------
v <- labelled(c(1,2,2,2,3,9,1,3,2,NA), c(yes = 1, no = 3, "don't know" = 8, refused = 9))
v

## ------------------------------------------------------------------------
val_labels(v)
val_label(v, 8)

## ------------------------------------------------------------------------
val_labels(v) <- c(yes = 1, nno = 3, bug = 5)
v
val_label(v, 3) <- "no"
v

## ------------------------------------------------------------------------
val_label(v, 2) <- "maybe"
val_label(v, 5) <- NULL
v

## ------------------------------------------------------------------------
val_labels(v) <- NULL
v

## ------------------------------------------------------------------------
val_label(v, 1) <- "yes"
v

## ------------------------------------------------------------------------
f <- factor(1:3)
f
val_labels(f) <- c(yes = 1, no = 3)
f

## ------------------------------------------------------------------------
df <- data.frame(v1 = 1:3, v2 = c(2, 3, 1), v3 = 3:1)

val_label(df, 1) <- "yes"
val_label(df[, c("v1", "v3")], 2) <- "maybe"
val_label(df[, c("v2", "v3")], 3) <- "no"
val_labels(df)

val_labels(df[, c("v1", "v3")]) <- c(YES = 1, MAYBE = 2, NO = 3)
val_labels(df)
val_labels(df) <- NULL
val_labels(df)
val_labels(df) <- list(v1 = c(yes = 1, no = 3), v2 = c(a = 1, b = 2, c = 3))
val_labels(df)

## ------------------------------------------------------------------------
v <- labelled(c(1,2,2,2,3,9,1,3,2,NA), c(yes = 1, no = 3, "don't know" = 8, refused = 9), c(FALSE, FALSE, TRUE, TRUE))
v

## ------------------------------------------------------------------------
missing_val(v)
missing_val(v) <- 9
v
missing_val(v) <- NULL
v
missing_val(v) <- c(8, 9)
v

## ---- eval=FALSE---------------------------------------------------------
#  missing_val(v) <- c(7, 8, 9)

## ------------------------------------------------------------------------
missing_val(v, force = FALSE) <- c(7, 8, 9)
v
missing_val(v, force = TRUE) <- c(7, 8, 9)
v

## ------------------------------------------------------------------------
missing_val(v)
val_label(v, 7) <- NULL
missing_val(v)

## ------------------------------------------------------------------------
v <- c(1,2,2,2,3,9,1,3,2,NA)
val_label(v, 1) <- "yes"
val_label(v, 3) <- "no"
val_label(v, 9) <- "refused"
val_label(v, 2) <- "maybe"
val_label(v, 8) <- "don't know"
v

## ------------------------------------------------------------------------
sort_val_labels(v)
sort_val_labels(v, decreasing = TRUE)

## ------------------------------------------------------------------------
sort_val_labels(v, according_to = "l")

## ------------------------------------------------------------------------
v <- labelled(c(1,2,2,2,3,9,1,3,2,NA), c(yes = 1, no = 3, "don't know" = 8, refused = 9), c(FALSE, FALSE, TRUE, TRUE))
v
missing_to_na(v)

## ------------------------------------------------------------------------
nolabel_to_na(v)

## ------------------------------------------------------------------------
size <- labelled(c(1.88, 1.62, 1.78, 99, 1.91), c("not measured" = 99))
size

## ------------------------------------------------------------------------
val_labels_to_na(size)

## ------------------------------------------------------------------------
v <- labelled(c(1,2,2,2,3,9,1,3,2,NA), c(yes = 1, no = 3, "don't know" = 8, refused = 9), c(FALSE, FALSE, TRUE, TRUE))
v
to_factor(v)

## ------------------------------------------------------------------------
to_factor(v, levels = "v")
to_factor(v, levels = "p")

## ------------------------------------------------------------------------
to_factor(v, ordered = TRUE)

## ------------------------------------------------------------------------
to_factor(v, missing_to_na = TRUE)
to_factor(missing_to_na(v))

## ------------------------------------------------------------------------
to_factor(v, sort_levels = "n")
to_factor(v, sort_levels = "v")
to_factor(v, sort_levels = "l")

## ------------------------------------------------------------------------
f <- factor(1:3, labels = c("a", "b", "c"))
to_labelled(f)

## ------------------------------------------------------------------------
v
to_labelled(to_factor(v))

## ---- eval=FALSE---------------------------------------------------------
#    # from foreign
#    library(foreign)
#    df <- to_labelled(read.spss(
#      "file.sav",
#      to.data.frame = FALSE,
#      use.value.labels = FALSE,
#      use.missings = FALSE
#   ))
#   df <- to_labelled(read.dta(
#     "file.dta",
#     convert.factors = FALSE
#   ))
#  
#   # from memisc
#   library(memisc)
#   nes1948.por <- UnZip("anes/NES1948.ZIP", "NES1948.POR", package="memisc")
#   nes1948 <- spss.portable.file(nes1948.por)
#   df <- to_labelled(nes1948)
#   ds <- as.data.set(nes19480)
#   df <- to_labelled(ds)

