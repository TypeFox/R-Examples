# Probability distribution and quantile functions
# These functions require the Paradigm4 extensions to SciDB.

setOldClass("phyper")
setGeneric("phyper", function(x, ...) stats::phyper(x,...))

setOldClass("qhyper")
setGeneric("qhyper", function(x, ...) stats::qhyper(x,...))

setMethod("phyper", signature(x="scidb_or_scidbdf"),
  function(x, q, m, n, k, new="p",`eval`=FALSE)
  {
    query = sprintf("hygecdf(%s,%s,%s,%s)",q,m,n,k)
    bind(x, new, query, `eval`=eval)
  })

setMethod("qhyper", signature(x="scidb_or_scidbdf"),
  function(x, p, m, n, k, new="q", `eval`=FALSE)
  {
    query = sprintf("ihygecdf(%s,%s,%s,%s)",p,m,n,k)
    bind(x, new, query, `eval`=eval)
  })

# Input
# a: SciDB array (scidb or scidbdf object)
# x: YES/YES count
# m,n,k: Marginal counts (see doc)
# alternative: {"two.sided", "greater", "less"}
scidb_fisher.test = function(a,x="x",m="m",n="n",k="k",alternative="two.sided", `eval`=FALSE)
{
  pvalname = make.unique_(a@attributes, "pval")
  oddsname = make.unique_(a@attributes, "estimate")
  query = sprintf("apply(%s, %s, fishertest_p_value(%s,%s,%s,%s,'%s'), %s, fishertest_odds_ratio(%s,%s,%s,%s))",
           a@name, pvalname, x,m,n,k,alternative, oddsname,x,m,n,k)
  .scidbeval(query,depend=list(a),`eval`=eval,gc=TRUE)
}
