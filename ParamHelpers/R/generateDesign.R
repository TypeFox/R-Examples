#FIXME: generateDesign will NOT work if there are dependencies
# over multiple levels of params and one only states the dependency only
#  wrt to the "last" param. also see daniels unit test.
#  it works as long all dependencies are stated, we need to at least document this

#FIXME: it really makes no sense to calculate the distance for params that are NA
# when we do the design and augment it right? think about what happens here


#' @title Generates a statistical design for a parameter set.
#'
#' @description
#' The following types of columns are created:
#' \tabular{ll}{
#'  numeric(vector)   \tab  \code{numeric}  \cr
#'  integer(vector)   \tab  \code{integer}  \cr
#'  discrete(vector)  \tab  \code{factor} (names of values = levels) \cr
#'  logical(vector)   \tab  \code{logical}
#' }
#' If you want to convert these, look at \code{\link[BBmisc]{convertDataFrameCols}}.
#' Dependent parameters whose constraints are unsatisfied generate \code{NA} entries in their
#' respective columns.
#' For discrete vectors the levels and their order will be preserved, even if not all levels are present.
#'
#' Currently only lhs designs are supported.
#'
#' The algorithm currently iterates the following steps:
#' \enumerate{
#'   \item{We create a space filling design for all parameters, disregarding \code{requires},
#'     a \code{trafo} or the forbidden region.}
#'   \item{Forbidden points are removed.}
#'   \item{Parameters are trafoed (maybe); dependent parameters whose constraints are unsatisfied
#'     are set to \code{NA} entries.}
#'   \item{Duplicated design points are removed. Duplicated points are not generated in a
#'    reasonable space-filling design, but the way discrete parameters and also parameter dependencies
#'    are handled make this possible.}
#'   \item{If we removed some points, we now try to augment the design in a space-filling way
#'     and iterate.}
#' }
#'
#' Note that augmenting currently is somewhat experimental as we simply generate missing points
#' via new calls to \code{\link[lhs]{randomLHS}}, but do not add points so they are maximally
#' far away from the already present ones. The reason is that the latter is quite hard to achieve
#' with complicated dependences and forbidden regions, if one wants to ensure that points actually
#' get added... But we are working on it.
#'
#' \code{generateDesign} will NOT work if there are dependencies over multiple levels of
#' parameters and the dependency is only given with respect to the \dQuote{previous} parameter.
#' A current workaround is to state all dependencies on all parameters involved.
#' (We are working on it.)
#'
#' @template arg_gendes_n
#' @template arg_parset
#' @param fun [\code{function}]\cr
#'   Function from package lhs.
#'   Possible are: \code{\link[lhs]{maximinLHS}}, \code{\link[lhs]{randomLHS}},
#'   \code{\link[lhs]{geneticLHS}}, \code{\link[lhs]{improvedLHS}}, \code{\link[lhs]{optAugmentLHS}},
#'   \code{\link[lhs]{optimumLHS}}
#'   Default is \code{\link[lhs]{randomLHS}}.
#' @param fun.args [\code{list}]\cr
#'   List of further arguments passed to \code{fun}.
#' @template arg_trafo
#' @param augment [\code{integer(1)}]\cr
#'   Duplicated values and forbidden regions in the parameter space can lead to the design
#'   becoming smaller than \code{n}. With this option it is possible to augment the design again
#'   to size \code{n}. It is not guaranteed that this always works (to full size)
#'   and \code{augment} specifies the number of tries to augment.
#'   If the the design is of size less than \code{n} after all tries, a warning is issued
#'   and the smaller design is returned.
#'   Default is 20.
#' @template ret_gendes_df
#' @export
#' @useDynLib ParamHelpers c_generateDesign c_trafo_and_set_dep_to_na
#' @examples
#' ps = makeParamSet(
#'   makeNumericParam("x1", lower = -2, upper = 1),
#'   makeIntegerParam("x2", lower = 10, upper = 20)
#' )
#' # random latin hypercube design with 5 samples:
#' generateDesign(5, ps)
#'
#' # with trafo
#' ps = makeParamSet(
#'   makeNumericParam("x", lower = -2, upper = 1),
#'   makeNumericVectorParam("y", len = 2, lower = 0, upper = 1, trafo = function(x) x/sum(x))
#' )
#' generateDesign(10, ps, trafo = TRUE)
generateDesign = function(n = 10L, par.set, fun, fun.args = list(), trafo = FALSE, augment = 20L) {

  n = asInt(n)
  z = doBasicGenDesignChecks(par.set)
  lower = z$lower
  upper = z$upper

  requirePackages("lhs", why = "generateDesign", default.method = "load")
  if (missing(fun))
    fun = lhs::randomLHS
  else
    assertFunction(fun)
  assertList(fun.args)
  assertFlag(trafo)
  augment = asInt(augment, lower = 0L)

  ### precompute some useful stuff
  pars = par.set$pars
  lens = getParamLengths(par.set)
  k = sum(lens)
  pids = getParamIds(par.set, repeated = TRUE, with.nr = TRUE)
  lower2 = setNames(rep(NA_real_, k), pids)
  lower2 = insert(lower2, lower)
  upper2 = setNames(rep(NA_real_, k), pids)
  upper2 = insert(upper2, upper)
  values = getParamSetValues(par.set)
  types.df = getParamTypes(par.set, df.cols = TRUE)
  types.int = convertTypesToCInts(types.df)
  types.df[types.df == "factor"] = "character"
  # ignore trafos if the user did not request transformed values
  trafos = if(trafo)
    lapply(pars, function(p) p$trafo)
  else
    replicate(length(pars), NULL, simplify = FALSE)
  par.requires = lapply(pars, function(p) p$requires)


  nmissing = n
  iter = 0
  # result objects
  res = data.frame()
  des = matrix(nrow = 0, ncol = k)
  repeat {
    ### get design, types converted, trafos, conditionals set to NA
    # create new design or augment if we already have some points
    newdes = if (nmissing == n)
      do.call(fun, insert(list(n = nmissing, k = k), fun.args))
    else
      lhs::randomLHS(nmissing, k = k)
    # preallocate result for C
    newres = makeDataFrame(nmissing, k, col.types = types.df)
    newres = .Call(c_generateDesign, newdes, newres, types.int, lower2, upper2, values)
    colnames(newres) = pids
    # check each row if forbidden, then remove
    if (hasForbidden(par.set)) {
      #FIXME: this is pretty slow, but correct
      fb = rowSapply(newres, isForbidden, par.set = par.set)
      newres = newres[!fb, , drop = FALSE]
      newdes = newdes[!fb, , drop = FALSE]
    }
    newres = .Call(c_trafo_and_set_dep_to_na, newres, types.int, names(pars), lens, trafos, par.requires, new.env())
    # add to result (design matrix and data.frame)
    des = rbind(des, newdes)
    res = rbind(res, newres)
    # remove duplicates
    to.remove = duplicated(res)
    des = des[!to.remove, , drop = FALSE]
    res = res[!to.remove, , drop = FALSE]
    nmissing = n - nrow(res)

    # enough points or augment tries? we are done!
    iter = iter + 1L
    if (nmissing == 0L || iter >= augment)
      break
  }

  if (nrow(res) < n)
    warningf("generateDesign could only produce %i points instead of %i!", nrow(res), n)

  colnames(res) = pids
  res = fixDesignFactors(res, par.set)
  attr(res, "trafo") = trafo
  return(res)
}
