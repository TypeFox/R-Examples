library(phyclust, quiet = TRUE)

### Examples to use ms().
set.seed(1234)

### An example from msdoc.pdf by Hudson, R.R.
tbs.matrix <- matrix(c(3.0, 3.5, 5.0, 8.5), nrow = 2)

### Generate an ancestral tree.
ret <- ms(nsam = 5, nreps = 2, opts = "-t tbs -r tbs 1000",
          tbs.matrix = tbs.matrix)
print(ret)
