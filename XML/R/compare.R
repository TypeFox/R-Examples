
# For comparing if two XML documents are "similar" whatever that means.
# We look at the distribution of node names

summary.XMLInternalDocument =
function(object, ...)
{
  counts = sort(table(xpathSApply(object, "//*", xmlName, ...)), decreasing = TRUE)
  list(nameCounts = counts,
       numNodes = sum(counts))
}


compareXMLDocs =
function(a, b, ...)
{
 sa = summary(a, ...)
 sb = summary(b, ...)

 inAOnly = setdiff(names(sa$nameCounts), names(sb$nameCounts))
 inBOnly = setdiff(names(sb$nameCounts), names(sa$nameCounts))

 common.ids = intersect(names(sa$nameCounts), names(sb$nameCounts)) #  != sb$nameCounts[names(sa$nameCounts)
 diffs = sa$nameCounts[common.ids] -  sb$nameCounts[common.ids]
 diffs = diffs[diffs != 0]
 
 list(inA = sa$nameCounts[inAOnly],  inB = sb$nameCounts[inBOnly], countDiffs = diffs)
 #all.equal(sa, sb)
}
