getTokenizers <-
function()
    c("MC_tokenizer", "scan_tokenizer")

# http://www.cs.utexas.edu/users/dml/software/mc/
MC_tokenizer <-
NLP::Token_Tokenizer(function(x)
{
    x <- as.character(x)
    ASCII_letters <- "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    id <- sprintf("[%s]+", ASCII_letters)
    http <- sprintf("(http://%s(\\.%s)*)", id, id)
    email <- sprintf("(%s@%s(\\.%s)*)", id, id, id)
    http_or_email <- sprintf("%s|%s", http, email)

    c(unlist(regmatches(x, gregexpr(http_or_email, x))),
      unlist(strsplit(gsub(http_or_email, "", x),
                      sprintf("[^%s]", ASCII_letters))))
})

scan_tokenizer <-
NLP::Token_Tokenizer(function(x)
    scan(text = as.character(x), what = "character", quote = "", quiet = TRUE))
