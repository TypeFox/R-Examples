nonce <-
function(length = 10) { 
    paste(sample(c(letters, LETTERS, 0:9), length, replace = TRUE), collapse = "") 
}
