# Case 1: define the reference set or lookup set within the function. Useful for small or binary
# sets like m(ale)/f(emale)
is_one_of("m", "m") == TRUE

is_one_of("m", c("f", "m")) == TRUE

is_one_of("y", c("f", "m")) == FALSE

is_one_of(c("b", "c", "d"), c("a", "b", "c", "d", "e")) == TRUE


# Case 2: use an external lookup table. The external lookup table must have at least one column
# called exactly 'VALUES'. May have also another one 'LABELS'. Useful for long lookup tables like
# list of countries.

# some preparation work for using a temporary directory
owd <- getwd()
td <- tempdir()
setwd(td)


VALUES <- LETTERS[1:10]
LABELS <- VALUES
db <- cbind(VALUES, LABELS)
db <- as.data.frame(db, stringsAsFactors = FALSE)
names(db) <- c("VALUES", "LABELS")
write.csv(db, "sample.csv", row.names = FALSE)

is_one_of("A", "sample.csv") == TRUE
is_one_of("Z", "sample.csv") == FALSE

# switching back to your working directory
setwd(owd)
 
