## Version 4 UUIDs have the form xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx
## where x is any hexadecimal digit and y is one of 8, 9, A, or B
## e.g., f47ac10b-58cc-4372-a567-0e02b2c3d479
uuid <- function(uppercase=FALSE) {
  
  hex_digits <- c(as.character(0:9), letters[1:6])
  hex_digits <- if (uppercase) toupper(hex_digits) else hex_digits
  
  y_digits <- hex_digits[9:12]
  
  paste(
    paste0(sample(hex_digits, 8, replace=TRUE), collapse=''),
    paste0(sample(hex_digits, 4, replace=TRUE), collapse=''),
    paste0('4', paste0(sample(hex_digits, 3, replace=TRUE), collapse=''), collapse=''),
    paste0(sample(y_digits,1), paste0(sample(hex_digits, 3, replace=TRUE), collapse=''), collapse=''),
    paste0(sample(hex_digits, 12, replace=TRUE), collapse=''),
    sep='-')
}
