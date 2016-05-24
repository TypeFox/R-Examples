## annotate() can use a single annotator or an annotator pipeline or
## something coercible to this, such as a list of annotators, and
## recursively calls the given annotators and merges annotations.

annotate <-
function(s, f, a = Annotation())
{
    s <- as.String(s)
    for(e in as.Annotator_Pipeline(f))
        a <- merge(a, e(s, a))

    a
}
