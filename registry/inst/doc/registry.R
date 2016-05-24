### R code from vignette source 'registry.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: registry.Rnw:56-58
###################################################
options(width = 80)
library("registry")


###################################################
### code chunk number 2: registry.Rnw:122-125
###################################################
library(registry)
R <- registry()
print(R)


###################################################
### code chunk number 3: registry.Rnw:141-142
###################################################
checkAge <- function(x) stopifnot(is.na(x) || x > 0 && x < 100)


###################################################
### code chunk number 4: registry.Rnw:147-148
###################################################
checkPhone <- function(x) stopifnot(!is.na(x$mobile) || !is.na(x$home))


###################################################
### code chunk number 5: registry.Rnw:153-155
###################################################
R <- registry(registry_class = "Addressbook", entry_class = "Address",
              validity_FUN = checkPhone)


###################################################
### code chunk number 6: registry.Rnw:158-164
###################################################
print.Addressbook <-
function(x, ...) {
    writeLines(sprintf("An address book with %i entries.\n", length(x)))
    invisible(x)
}
print(R)


###################################################
### code chunk number 7: registry.Rnw:170-172
###################################################
R$set_field("last", type = "character", is_key = TRUE, index_FUN = match_partial_ignorecase)
R$set_field("first", type = "character", is_key = TRUE, index_FUN = match_partial_ignorecase)


###################################################
### code chunk number 8: registry.Rnw:175-176
###################################################
R$set_field("address", type = "character")


###################################################
### code chunk number 9: registry.Rnw:181-183
###################################################
R$set_field("mobile", type = "character")
R$set_field("home", type = "character")


###################################################
### code chunk number 10: registry.Rnw:187-188
###################################################
R$set_field("age", type = "integer", validity_FUN = checkAge)


###################################################
### code chunk number 11: registry.Rnw:192-195
###################################################
R$set_field("type", type = "character",
            alternatives = c("Business", "Private"),
            default = "Business")


###################################################
### code chunk number 12: registry.Rnw:198-199
###################################################
R$get_field("type")


###################################################
### code chunk number 13: registry.Rnw:206-210
###################################################
R$set_entry(last = "Smith", first = "Mary", address = "Vienna",
            home = "734 43 34", type = "Private", age = 44L)
R$set_entry(last = "Smith", first = "Peter", address = "New York",
            mobile = "878 78 87")


###################################################
### code chunk number 14: registry.Rnw:213-215
###################################################
R$set_entry("Myers", "John", "Washington", "52 32 34", "898 89 99",
            33L, "Business")


###################################################
### code chunk number 15: registry.Rnw:218-223
###################################################
TRY <- function(expr) tryCatch(expr, error = print)
TRY(R$set_entry(last = "Smith", first = "Mary"))
TRY(R$set_entry(last = "Miller", first = "Henry"))
TRY(R$set_entry(last = "Miller", first = "Henry", age = 12.5))
TRY(R$set_entry(last = "Miller", first = "Henry", age = 999L))


###################################################
### code chunk number 16: registry.Rnw:226-227
###################################################
R$get_entry(last = "Smith", first = "mar")


###################################################
### code chunk number 17: registry.Rnw:231-234
###################################################
print.Address <- function(x) with(x,
    writeLines(sprintf("%s %s, %s; home: %s, mobile: %s; age: %i (%s)", first, last, address, home, mobile, age, type)))
R$get_entry(last = "Smith", first = "mar")


###################################################
### code chunk number 18: registry.Rnw:241-242
###################################################
R[["Myers"]]


###################################################
### code chunk number 19: registry.Rnw:246-248
###################################################
R$set_entry(last = "Frears", first = c("Joe", "Jonathan"),
            address = "Washington", home = "721 42 34")


###################################################
### code chunk number 20: registry.Rnw:251-252
###################################################
identical(R[["Frears", "Jonathan"]], R[["Frears", "Joe"]])


###################################################
### code chunk number 21: registry.Rnw:258-259
###################################################
R$get_entries("Smith")


###################################################
### code chunk number 22: registry.Rnw:262-263
###################################################
R$grep_entries("Priv")


###################################################
### code chunk number 23: registry.Rnw:266-268 (eval = FALSE)
###################################################
## R$get_entries()
## R[]


###################################################
### code chunk number 24: registry.Rnw:271-272
###################################################
summary(R)


###################################################
### code chunk number 25: registry.Rnw:276-279
###################################################
R[["Smith", "Peter"]]
R$modify_entry(last = "Smith", first = "Peter", age = 22L)
R[["Smith", "Peter"]]


###################################################
### code chunk number 26: registry.Rnw:282-284
###################################################
R$delete_entry(last = "Smith", first = "Peter")
R[["Smith", "Peter"]]


###################################################
### code chunk number 27: registry.Rnw:294-299
###################################################
R$seal_entries()
TRY(R$delete_entry("Smith", "Mary"))
R$set_entry(last = "Slater", first = "Christian", address = "Boston",
            mobile = "766 23 88")
R[["Slater"]]


###################################################
### code chunk number 28: registry.Rnw:302-307
###################################################
R$get_permissions()
R$restrict_permissions(delete_entries = FALSE)
TRY(R$delete_entry("Slater"))
R$modify_entry(last = "Slater", first = "Christian", age = 44L)
R[["Slater"]]


