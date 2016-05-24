## Simple validating text entry

library(qtbase)

qsetClass("PositiveValidator", Qt$QValidator)

qsetMethod("validate", positiveValidator, function(input, pos) {
  val <- suppressWarnings(as.integer(input))
  if (!is.na(val)) {
    if (val > 0)
      Qt$QValidator$Acceptable
    else Qt$QValidator$Invalid
  } else {
    if (input == "")
      Qt$QValidator$Acceptable
    else Qt$QValidator$Invalid
  }
})

e <- Qt$QLineEdit()
v <- PositiveValidator(e)
e$setValidator(v)
e$show()
