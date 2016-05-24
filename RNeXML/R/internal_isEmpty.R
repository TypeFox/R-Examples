
isEmpty <- function (obj) 
{
    if (!isS4(obj)) {
        if (length(obj) > 0) 
            FALSE
        else TRUE
    }
    else {
        if (identical(obj, new(class(obj)[1]))) 
            out <- TRUE
        else {
            empty <- sapply(slotNames(obj), function(s) {
                if (isS4(slot(obj, s))) 
                  isEmpty(slot(obj, s))
                else {
                  if (length(slot(obj, s)) == 0) 
                    TRUE
                  else if (length(slot(obj, s)) > 0) 
                    FALSE
                }
            })
            out <- !any(!empty)
        }
        out
    }
}

