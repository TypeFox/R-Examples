# This function is used to initialize the "seqInfo" object in the
# .odfEnv environment.  This must be done before the Sweave phase.
# It used to be done as part of the preprocessing phase.  We could
# continue to do then, but it isn't logically part of preprocessing.
# It was just convenient, and perhaps confusing, to do it then.

getSeqInfo <- function(node)
{
   sfun <- function(s) xmlGetAttr(s, 'text:ref-name')
   refnames <- as.character(treeapply(node, 'text:sequence', sfun, onlyFirst=FALSE, rooted=FALSE))
   list(Table=refnames, Illustration=refnames)
}
