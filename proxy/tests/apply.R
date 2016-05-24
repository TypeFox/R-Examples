
## tests on apply C wrappers

library(proxy)

set.seed(20070630)

## matrix

f <- function(x, y) sum(x*y) / sqrt(sum(x*x)) / sqrt(sum(y*y))

x <- matrix(runif(20), 5, 4)
x
y <- matrix(runif(20), 5, 4)
y

.External("R_apply_dist_matrix", x, NULL, FALSE, f)
.External("R_apply_dist_matrix", x, x, FALSE, f)
.External("R_apply_dist_matrix", x, y, FALSE, f)
.External("R_apply_dist_matrix", x, y, TRUE, f)

# coerce

z <- y * 100
storage.mode(z) <- "integer"
.External("R_apply_dist_matrix", x, z, FALSE, f)
.External("R_apply_dist_matrix", z, x, FALSE, f)
.External("R_apply_dist_matrix", z, z, FALSE, f)
.External("R_apply_dist_matrix", z, NULL, FALSE, f)

## list

x <- unlist(apply(x, 1, list), recursive = FALSE)
x
y <- unlist(apply(y, 1, list), recursive = FALSE)

.External("R_apply_dist_list", x, NULL, FALSE, f)
.External("R_apply_dist_list", x, x, FALSE, f)
.External("R_apply_dist_list", x, y, FALSE, f)
.External("R_apply_dist_list", x, y, TRUE, f)

## logical matrix

f <- function(a, b, c, d, n)
    a / sqrt(a+b) / sqrt(a+c)

x <- t(sapply(x, ">", 0.5))
x
y <- t(sapply(y, ">", 0.5))

.External("R_apply_dist_binary_matrix", x, NULL, FALSE, f)
.External("R_apply_dist_binary_matrix", x, x, FALSE, f)
.External("R_apply_dist_binary_matrix", x, y, FALSE, f)
.External("R_apply_dist_binary_matrix", x, y, TRUE, f)

## data.frame

f <- function(x, y) sum(x*y) / sqrt(sum(x*x)) / sqrt(sum(y*y))

x <- data.frame(unlist(apply(x, 2, list), recursive = FALSE))
names(x) <- letters[1:4]
x
y <- data.frame(unlist(apply(y, 2, list), recursive = FALSE))
names(y) <- letters[1:4]

.External("R_apply_dist_data_frame", x, NULL, FALSE, f)
.External("R_apply_dist_data_frame", x, x, FALSE, f)
.External("R_apply_dist_data_frame", x, y, FALSE, f)
.External("R_apply_dist_data_frame", x, y, TRUE, f)

#

f <- function(x, y) {
    if (rownames(x) == 1 && rownames(y) == 1) {
        print(x)
        str(x)
        print(y)
    }
    sum(x == y) / length(x)
}

x <- data.frame(1:5, LETTERS[1:5])
x

y <- data.frame(1:6, LETTERS[c(1,1:5)], row.names = letters[1:6])
y

all.equal(x, y)
identical(attributes(x[[1]]), attributes(y[[1]]))
identical(attributes(x[[2]]), attributes(y[[2]]))

.External("R_apply_dist_data_frame", x, NULL, FALSE, f)
.External("R_apply_dist_data_frame", x, x, FALSE, f)
.External("R_apply_dist_data_frame", x, y, FALSE, f)
.External("R_apply_dist_data_frame", x, y[-6,], TRUE, f)


# test parameters

f <- function(x, y, p = 1)
    sum(x == y) / length(x) * p

.External("R_apply_dist_data_frame", x, y[-6,], TRUE, f)
.External("R_apply_dist_data_frame", x, y[-6,], TRUE, f, p = 2)

###
