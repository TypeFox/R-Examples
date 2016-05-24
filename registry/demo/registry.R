##########################
### registry test instances

library(registry)

.my_check_fun <- function(x) if (x$Z == 999 && x$New2 == 999) stop("No evil allowed!")

## create registry
R <- registry(entry_class = "simple.list",
              validity_FUN = .my_check_fun)
R

## set fields
R$set_field("names", type = "character", is_key = TRUE)
R$set_field("X", type = TRUE, is_mandatory = TRUE)
R$set_field("Y", type = "character")
R$set_field("Z", default = 123)
R$get_fields()

## add entries
R$set_entry(names = c("test", "bar"), X = TRUE, Y = "bla")
R$set_entry(names = "test2", X = FALSE, Y = "foo", Z = 99)
R$set_entry(names = "test3", X = FALSE, Y = "bar", Z = "chars")
R$get_entry("test")
R[["test2"]]
R[["test3"]]
R[["bar"]]

## add new field
R$set_field("New")
R$get_field("New")

## change entries
R$modify_entry(names = "test", New = 123)
R$modify_entry(names = "test2", New = "test")

## field check function (checks for strict positive values)
R$set_field("New2", type = "numeric", validity_FUN = function(x) stopifnot(x > 0))
R$set_entry(names = "test5", X = TRUE, New2 = 2)

## add field with fixed alternatives
R$set_field("New3", type = "character", alternatives = c("A", "B"))
R$get_field("New")
R$set_entry(names = "test6", X = TRUE, New3 = "A")

## print/summary = as.data.frame
R
summary(R)

## seal entries
R$seal_entries()
R$set_field("New4")
R$set_entry(names = "test7", X = TRUE, Y = "bla")
R$delete_entry("test7")
R$modify_entry(names = "test", New4 = "test")

## error cases:
TRY <- function(...) stopifnot(inherits(try(..., silent = TRUE), "try-error"))
TRY(R$set_field("bla", type = "character", default = 123))
TRY(R$set_entry("err1", Y = "bla"))
TRY(R$set_entry("err2", X = "bla"))
TRY(R$set_entry("err3", X = TRUE, New2 = -2))
TRY(R$set_entry("err4", X = TRUE, Z = 999, New2 = 999))
TRY(R$set_entry("err5", X = TRUE, New3 = "C"))
TRY(R$modify_entry("Bla", "New", 123))
TRY(R$modify_entry("X", "Bla", 123))
TRY(R$modify_entry("test","X",TRUE))
