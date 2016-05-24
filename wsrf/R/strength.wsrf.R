strength <- function(object, ...) UseMethod("strength")

strength.wsrf <- function (object, ...) {

  object[[.STRENGTH_IDX]]

}