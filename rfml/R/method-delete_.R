## -----------------------------------------------------------------------
## delete ml.data.frame objects
## -----------------------------------------------------------------------
# Currently not used!
setGeneric (
  "delete",
  def = function (x, ...) {
    rst <- standardGeneric("delete")
  },
  signature = "x")

## -----------------------------------------------------------------------

setMethod (
  "delete",
  signature (x = "ml.data.frame"),
  def = function (x, cascade = FALSE) {

  })

## -----------------------------------------------------------------------

setMethod (
  "delete",
  signature (x = "character"),
  def = function (x, conn.id = 1, is.temp = FALSE, cascade = FALSE) {

  })
