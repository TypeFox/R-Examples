is.diffproc <-
function (obj) 
{
    if (inherits(obj, "diffproc") & is.list(obj) & (length(obj) == 
        4)) 
        if (all(lapply(obj, length) == 1) & all(unlist(lapply(obj, 
            is.character))) & identical(names(obj), c("mean", 
            "var", "tpdf", "tpdF"))) 
            return(all(unlist(lapply(obj[2:3], function(l) !inherits(try(D(parse(text = l), 
                "x"), silent = TRUE), "try-error")))) & all(unlist(lapply(obj[c(1, 
                4)], function(l) !inherits(try(parse(text = l), 
                silent = TRUE), "try-error")))))
    return(FALSE)
}
