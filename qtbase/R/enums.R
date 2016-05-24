## convenience functions for bitwise operations on enums/flags

## A QtEnum is a enum value represented by an R integer vector.

## QFlags arise through '|' combination of enum values. Smoke
## communicates QFlags as plain unsigned integers. We could use an R
## double vector to avoid misinterpretation of the sign bit; however,
## QFlags properties are not happy when we pass them doubles instead
## of integers. Thus, we keep enums/flags as integers.

"|.QtEnum" <- function(x, y) {
  structure(packBits(intToBits(x) | intToBits(y), "integer"),
            class = "QtEnum")
}

"&.QtEnum" <- function(x, y) {
  structure(packBits(intToBits(x) & intToBits(y), "integer"),
            class = "QtEnum")
}

print.QtEnum <- function(x, ...) {
  cat("Enum value: ", names(x), " (", x, ")\n", sep = "")
}
