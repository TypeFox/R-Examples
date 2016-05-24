# convert clustering result to label vector
setMethod("labels", signature(object="ExClust"),
    function(object, type="names")
    {
        if (type == "names")
        {
            if (length(names(object@idx)) == 0)
                stop("no names available, use other type")
            else
                out <- names(object@idx)
        }
        else if (type == "exemplars")
            out <- object@idx
        else if (type == "enum")
        {
            out <- array(dim=object@l)

            for (i in 1:length(object@exemplars))
                out[which(object@idx == object@exemplars[i])] <- i
        }
        else
            stop("type '", type, "' unknown")

        attributes(out) <- NULL

        out
    }
)
