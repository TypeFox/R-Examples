###-----------------------------------------------------------------------------
###
###  .paste.row
###
###

.paste.row <- function (x, label = "", name = "", ins = 0, i = NULL, val = FALSE,
                        ...){
  .print.val <- function (x, range = TRUE, digits = getOption ("digits"),
                          max.print = 5, shorten.to = c (2,1)){
    if (is.list (x)){                   # also for data.frames 
       paste ("", "columns/entries", paste (names (x), collapse = ", "))
     } else {
       if (length (x) == 0)
         return ("")
       
       if (any (is.na (x)))
         text <- "+ NA"
       else
         text <- ""
       
       if (range) 
         x <- sort (unique (as.vector (x)))
       else
         x <- as.vector (x)
       
       if (length (x) > max.print){
         from <- format (head (x, shorten.to [1]), digits = digits, trim = TRUE)
         to <- format (tail (x, shorten.to [2]), digits = digits, trim = TRUE)
         
         text <- paste (paste (from, collapse = " "), "...",
                        paste (to, collapse = " "), text, collapse = " ")
       } else {
         text <- paste (paste (format (x, digits = digits, trim = TRUE),
                               collapse = " "),
                        text,  collapse = " ")
       }
       
       paste (if (range) " rng ", text, collapse = "")
       
     }
  }
  
  label <- paste (as.character (label), "", collapse = " ")
  
  paste (paste (rep (" ", ins), collapse = ""),
         if (!is.null (i)) paste(i, ". ", collapse = "", sep = ""),
         name,
         ": ",
         label,
         "[",
         paste (class (x), collapse = ", "),
         if (! is.null (dim (x)))
         paste (if (is.matrix (x) & all (class (x) != "matrix")) " matrix x " else
                if (is.array (x) & all (class (x) != "array") & all (class (x) != "matrix"))
                " array x ",
                paste (dim (x) [-1], collapse = " x ")
                , sep = ""),
         "]",
         if (val) .print.val (x, ...),
         sep ="", collapse = "")
}
