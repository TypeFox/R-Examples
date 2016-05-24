cleanstring <-
# helper function for extracting variable names from model objects 
function(string)
{
    if(!is.character(string)) stop("'string' should be of type character.")
    s <- unlist(strsplit(string, split = ""))
    if (any(c("(", ")", "^") %in% s)){
        beg <- if ("(" %in% s)
            which(s == "(")[1] + 1
        else 1
        end <- if (any(c(")", "^") %in% s))
            which(s == ")" | s == "^")[1] - 1
        else length(s)
        s <- s[beg:end]
    }
    paste(s, collapse = "")
}