.makeObj.tm <- function() {
  tm <- list()
  class(tm) <- "tm"
  tm
}

##' Check whether an object is of type tm (tree model)
##' @title Tree Models
##' @param x an object to be tested
##' @return TRUE if an object has class "tm", FALSE otherwise
##' @keywords tm
##' @export
##' @author Melissa J. Hubisz
is.tm <- function(x) {
  return(class(x)=="tm")
}


.rootLeaf.from.pointer.tm <- function(x, tree) {
  if (is.null(x$externalPtr))
    stop(".rootLeaf.from.pointer.tm expects list with externalPtr")
  id <- .Call.rphast("rph_tm_rootLeaf", x$externalPtr)
  if (is.null(id) || id == -1) return(NULL)
  .Call.rphast("rph_tree_nodeName", tree, id)
}


# NOTE: tm's are only stored as external pointers with RPHAST functions,
# never by the user.  So from.pointer and as.pointer are internal functions.
# Treemodels stored as external pointers are NOT protected.
# from.pointer and as.pointer do not call freeall.rphast

from.pointer.one.lsmodel.tm <- function(x, i) {
  rv <- list()
  rv[["defn"]] <- .Call.rphast("rph_tm_altmodel_def", x$externalPtr, i)
  rv[["subst.mod"]] <- .Call.rphast("rph_tm_altmodel_substMod", x$externalPtr, i)
  val <- .Call.rphast("rph_tm_altmodel_backgd", x$externalPtr, i)
  if (!is.null(val)) rv[["backgd"]] <- val
  rv[["rate.matrix"]] <- rphast.simplify.list(.Call.rphast("rph_tm_altmodel_rateMatrix",
                                                    x$externalPtr, i))
  val <- .Call.rphast("rph_tm_altmodel_sel", x$externalPtr, i)
  if (!is.null(val))
    rv[["selection"]] <- val[1]
  val <- .Call.rphast("rph_tm_altmodel_bgc", x$externalPtr, i)
  if (!is.null(val))
    rv[["bgc"]] <- val[1]
  rv
}

from.pointer.lsmodel.tm <- function(x) {
  numModel <- .Call.rphast("rph_tm_num_altmodel", x$externalPtr)
  if (numModel == 0L) return(NULL)
  rv <- list()
  for (i in 1:numModel) {
    rv[[i]] <- from.pointer.one.lsmodel.tm(x, i)
    if (numModel == 1L) return(rv[[i]])
  }
  rv
}


from.pointer.tm <- function(x) {
  if (is.null(x$externalPtr))
    stop("from.pointer.tm expects list with externalPtr")
  tm <- .makeObj.tm()
  tm$alphabet <- .Call.rphast("rph_tm_alphabet", x$externalPtr)
  tm$backgd <- .Call.rphast("rph_tm_backgd", x$externalPtr)
  tm$rate.matrix <- rphast.simplify.list(.Call.rphast("rph_tm_rateMatrix", x$externalPtr))
  tm$subst.mod <- .Call.rphast("rph_tm_substMod", x$externalPtr)
  tm$likelihood <- .Call.rphast("rph_tm_likelihood", x$externalPtr)
  tm[["alpha"]]<- .Call.rphast("rph_tm_alpha", x$externalPtr)
  tm$nratecats <- .Call.rphast("rph_tm_nratecats", x$externalPtr)
  tm$rate.consts <- .Call.rphast("rph_tm_rK", x$externalPtr)
  tm$rate.weights <- .Call.rphast("rph_tm_freqK", x$externalPtr)
  tm$tree <- fix.semicolon.tree(.Call.rphast("rph_tm_tree", x$externalPtr))
  tm$root.leaf <- .rootLeaf.from.pointer.tm(x, tm$tree)
  selection <- .Call.rphast("rph_tm_selection", x$externalPtr)
  if (selection[1]==1) tm$selection <- selection[2]
  site.model <- .Call.rphast("rph_tm_site_model", x$externalPtr)
  if (site.model) tm$site.model <- site.model
  lsmodel <- from.pointer.lsmodel.tm(x)
  if (!is.null(lsmodel)) tm$ls.model <- lsmodel
  tm
}


as.pointer.tm <- function(x) {
  rv <- list()
  rv$externalPtr <- .Call.rphast("rph_tm_new",
                                 x$tree,
                                 x$alphabet,
                                 x$backgd,
                                 x$rate.matrix,
                                 x$subst.mod,
                                 x$likelihood,
                                 x[["alpha"]],
                                 x$nratecats,
                                 x$rate.consts,
                                 x$rate.weights,
                                 x$root.leaf,
                                 x$selection,
                                 x$site.model)
  if (!is.null(x$ls.model)) {
    if (!is.null(x$ls.model$defn)) {
      numModel <- 1L
    } else numModel <- length(x$ls.model)
    for (i in 1:numModel) {
      if (numModel == 1L && !is.null(x$ls.model$defn)) {
        lsmod <- x$ls.model
      } else lsmod <- x$ls.model[[i]]

      .Call.rphast("rph_tm_add_alt_mod", rv$externalPtr, lsmod[["defn"]])
      .Call.rphast("rph_tm_altmod_set_subst_mod", rv$externalPtr, i,
                   lsmod[["subst.mod"]])
      .Call.rphast("rph_tm_altmod_set_backgd", rv$externalPtr, i,
                   lsmod[["backgd"]])
      .Call.rphast("rph_tm_altmod_set_ratematrix", rv$externalPtr, i,
                   lsmod[["rate.matrix"]])
      .Call.rphast("rph_tm_altmod_set_sel_bgc", rv$externalPtr, i,
                   lsmod[["selection"]], lsmod[["bgc"]])
    }
  }
  rv
}


# should we check to make sure it has all the necessary elements and no extra ones?
as.tm.list <- function(l) {
  class(l) <- "tm"
  l
}

##' Read a tree model from a file
##' @param filename The file containing a tree model
##' @title Read a Tree Model
##' @return An object of class "tm"
##' @seealso \code{\link{tm}}
##' @keywords tm
##' @export
##' @example inst/examples/read-tm.R
##' @author Melissa J. Hubisz and Adam Siepel
read.tm <- function(filename) {
  check.arg(filename, "filename", "character", null.OK=FALSE)
  x <- list()
  x$externalPtr <- .Call.rphast("rph_tm_read", filename)
  tm <- from.pointer.tm(x)
  tm
}


##' Write a tree model to a file (or to the terminal)
##' @title Wrting Tree Models
##' @param tm An object of class "tm"
##' @param filename The filename to write to (use NULL for output to terminal)
##' @param append Whether to append the tree to the end of the file
##' (if FALSE, overwrites file).  Not used if filename is \code{NULL}
##' @seealso \code{\link{tm}}
##' @keywords tm
##' @export
##' @example inst/examples/write-tm.R
##' @author Melissa J. Hubisz and Adam Siepel
write.tm <- function(tm, filename=NULL, append=FALSE) {
  check.arg(filename, "filename", "character", null.OK=TRUE)
  check.arg(append, "append", "logical", null.OK=TRUE)
  tm <- as.pointer.tm(tm)
  invisible(.Call.rphast("rph_tm_print", tm$externalPtr, filename, append))
}


##' Tree model summary
##' @title Tree Model Summary
##' @param object An object of class tm
##' @param ... Parameters to be passed to/from other functions
##' @seealso \code{\link{tm}}
##' @export
##' @export summary.tm
##' @keywords tm
##' @method summary tm
##' @example inst/examples/summary-tm.R
##' @author Melissa J. Hubisz and Adam Siepel
summary.tm <- function(object, ...) {
  write.tm(object, NULL)
}


##' Coerce a tree model into a list
##' @title Tree Model to List
##' @param x an object of class tm
##' @param ... arguments to be passed to/from other functions
##' @seealso \code{\link{tm}}
##' @export
##' @export as.list.tm
##' @method as.list tm
##' @example inst/examples/as-list-tm.R
##' @author Melissa J. Hubisz and Adam Siepel
as.list.tm <- function(x, ...) {
  class(x) <- "list"
  x
}

##' Print a tree model
##' @title Printing Tree Models
##' @param x An object of class \code{tm}.
##' @param aslist Logical.  If \code{TRUE}, print the tree model as a list
##' rather than in tree model format.
##' @param ... arguments to be passed to/from other functions
##' @seealso \code{\link{tm}}
##' @export
##' @export print.tm
##' @keywords tm
##' @method print tm
##' @example inst/examples/print-tm.R
##' @author Melissa J. Hubisz and Adam Siepel
print.tm <- function(x, aslist=FALSE, ...) {
  if (aslist) print(as.list(x), ...)
  else summary.tm(x, ...)
}


##' Check if a string represents a phast substitution model
##' @title Check Substitution Model Strings
##' @param mod A vector of character strings representing substitution
##' model names to be tested
##' @return A vector of logical values indicating whether each string
##' represents a defined substitution model
##' @export
##' @example inst/examples/is-subst-mod-tm.R
##' @author Melissa J. Hubisz and Adam Siepel
is.subst.mod.tm <- function(mod) {
  result <- logical(length(mod))
  for (i in 1:length(mod)) 
    result[i] <- .Call.rphast("rph_subst_mods_is_valid_string", mod[i])
  result
}

##' List all valid substitution models
##' @title List PHAST Substitution Models
##' @return a character vector with the names of all valid substitution
##' models
##' @export
##' @example inst/examples/subst-mods.R
##' @author Melissa J. Hubisz and Adam Siepel
subst.mods <- function() {
  .Call.rphast("rph_subst_mods_list_all", NULL)
}
    


##' Make a new tree model
##'
##' Tree models represent a substitution process along a phylogenetic
##' tree.  They are stored as a list, with components defined by the
##' arguments to this function.
##' @title Tree Models
##' @param tree A character string representing a phylogenetic tree in
##' newick foromat
##' @param subst.mod A character string giving a valid substitution mod.
##' See \code{\link{subst.mods}}.
##' @param rate.matrix A square matrix representing the rate of substitution
##' from one state to the next.
##' @param backgd A numeric vector giving the equilibrium frequencies for
##' each state.
##' @param alphabet A character vector containing all valid states, given
##' in the order they are represented in rate.matrix and backgd.  Defaults
##' to "ACGT"
##' @param nratecats The number of rate categories in the model.  Defaults
##' to 1.
##' @param alpha  If nratecats > 1, weight for each category is computed using
##' a gamma distribution with shape parameter alpha.
##' @param rate.consts The rate for each rate category.  NULL if only
##' one category.
##' @param rate.weights Vector of numeric of length nratecats, determining
##' the weight of each rate category.  Must sum to 1 (will be normalized
##' otherwise).  May be defined implicitly by alpha.
##' @param selection If not NULL, then this is a numeric value giving the
##' selection parameter for this model.  If NULL then there is no selection
##' in the model.  If selection==0.0, means that selection has no effect
##' in the current model, but is part of the model, and by default the
##' selection parameter will be optimized by phyloFit.  The rate matrix
##' is assumed to already be scaled by the selection parameter, if provided.
##' @param root.leaf Usually NULL, but if set to the name of a leaf
##' node in the tree, the tree will be re-rooted at this leaf node.
##' @param likelihood an optional value giving the log likelihood of this
##' model for some alignment.
##' @keywords tm
##' @return An object of class \code{tm} representing a phylogenetic model.
##' @export
##' @example inst/examples/tm.R
##' @author Melissa J. Hubisz and Adam Siepel
tm <- function(tree, subst.mod, rate.matrix=NULL, backgd=NULL,
               alphabet="ACGT", nratecats=1, 
               alpha=0.0, rate.consts=NULL, rate.weights=NULL,
               selection=NULL,
               root.leaf=NULL, likelihood=NULL) {
  check.arg(tree, "tree", "character", null.OK=FALSE)
  check.arg(subst.mod, "subst.mod", "character", null.OK=FALSE)
  check.arg(rate.matrix, "rate.matrix", "numeric", null.OK=TRUE,
            min.length=NULL, max.length=NULL)  #will check later
  check.arg(backgd, "backgd", "numeric", null.OK=TRUE,
            min.length=NULL, max.length=NULL) # will check later
  check.arg(alphabet, "alphabet", "character", null.OK=FALSE)
  check.arg(nratecats, "nratecats", "integer", null.OK=FALSE)
  check.arg(likelihood, "likelihood", "numeric", null.OK=TRUE)
  check.arg(alpha, "alpha", "numeric", null.OK=FALSE)
  check.arg(rate.consts, "rate.consts", "numeric", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(selection, "selection", "numeric", null.OK=TRUE,
            min.length=1L, max.length=1L)
  check.arg(root.leaf, "root.leaf", "character", null.OK=TRUE)

  if (!is.subst.mod.tm(subst.mod))
    stop("invalid subst mod ", subst.mod)
  matsize <- NULL
  if (!is.null(rate.matrix)) {
    matsize <- as.integer(sqrt(length(rate.matrix)))
    if (matsize*matsize != length(rate.matrix))
      stop("rate.matrix should be square matrix")
    if (!is.matrix(rate.matrix))
      rate.matrix <- matrix(rate.matrix, nrow=matsize, ncol=matsize, byrow=TRUE)
  }
  if (!is.null(backgd)) {
    if (!is.null(matsize)) {
      if (matsize != length(backgd))
        stop("length of backgd should match dimensions of rate.matrix")
    }
  }
  if (!is.null(rate.consts)) {
    if (nratecats != length(rate.consts))
      stop("rate.consts should have length equal to nratecats")
  }
  if (!is.null(rate.weights)) {
    if (is.null(rate.consts))
      stop("rate.weights needs to be specified if rate.consts is specified")
    if (nratecats != length(rate.weights))
      stop("rate.weights should have length equal to nratecats")
  }
  
  if (!is.null(root.leaf)) {
    if ( ! (.Call.rphast("rph_tree_isNode", tree, root.leaf))) {
      stop("tree has no node named ", root.leaf)
    }
  }

  tm <- .makeObj.tm()
  tm$alphabet <- alphabet
  tm$backgd <- backgd
  tm$rate.matrix <- rate.matrix
  tm$subst.mod <- subst.mod
  tm[["alpha"]] <- alpha
  tm$likelihood <- likelihood
  tm$nratecats <- nratecats
  tm$rate.consts <- rate.consts
  tm$rate.weights <- rate.weights
  tm$tree <- fix.semicolon.tree(tree)
  tm$selection <- selection
  tm$root.leaf <- root.leaf
  tm
}


##' Adjust tree model background frequencies while maintaining reversibility
##' @param tm An object of type \code{tm}
##' @param new.backgd A numeric vector of length 4 giving the background
##' frequencies of A,C,G,T
##' @param gc (Alternative to new.backgd) A numeric value giving the GC
##' content, which is used to calculate new background frequencies.
##' Assumes freq(C)==freq(G)==gc/2, and freq(A)==freq(T)==(1-gc)/2
##' @note Currently only works with models of order 0, without lineage-
##' specific models, and which use the default alphabet "ACGT".
##' @export
##' @example inst/examples/mod-backgd-tm.R
##' @author Melissa J. Hubisz and Adam Siepel
mod.backgd.tm <- function(tm, new.backgd=NULL, gc=NULL) {
  if (is.null(new.backgd)) {
    if (is.null(gc)) stop("either new.backgd or gc should be provided")
    if (length(gc) != 1L || gc < 0 || gc > 1) stop("invalid value of gc")
    new.backgd <- c((1-gc)/2, gc/2, gc/2, (1-gc)/2)
  }
  if (length(new.backgd) != 4 ||
      sum(new.backgd < 0) > 0 ||
      sum(new.backgd > 1) > 0)
    stop("invalid background frequencies")
  tm <- as.pointer.tm(tm)
  .Call.rphast("rph_tm_mod_freqs", tm$externalPtr, new.backgd)
  tm <- from.pointer.tm(tm)
  tm
}


##' BGC+selection factor
##' @param x The cumulative effect of bgc and selection.  If bgc and sel are
##' population-scaled parameters describing biased gene conversion and selection,
##' then x should be sel+bgc for strong mutations, sel-bgc for weak mutations, a
##' and sel for neutral mutations.
##' @return x/(1-e^(-x)), the factor to scale the rate matrix entry by.
##' @author Melissa J. Hubisz
##' @export
bgc.sel.factor <- function(x) {
  ifelse(abs(x) < 1.0e-6,
         1 + x/2 + x*x/12,
         x/(1-exp(-x)))
}



##' Apply bgc+selection parameters to a matrix
##' @param m A transition matrix
##' @param bgc The bgc (biased gene conversion) parameter, population-scaled.
##' @param sel The selection parameter (population-scaled)
##' @param alphabet The alphabet used for nucleotide states
##' @return A matrix with bgc+sel applied.  This matrix may no longer be
##' time reversible.
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
apply.bgc.sel <- function(m, bgc=0, sel=0, alphabet="ACGT") {
  if (is.null(bgc)) bgc <- 0
  if (is.null(sel)) sel <- 0
  rphast.simplify.list(.Call.rphast("rph_tm_apply_selection_bgc",
                                    as.matrix(m), alphabet, sel, bgc))
}


##' Unapply bgc+selection parameters from a matrix
##' @param m A transition matrix
##' @param bgc The bgc parameter which was used to calculate m
##' @param sel The selection parameter which was used to calculate m
##' @param alphabet The alphabet used for nucleotide states
##' @return A matrix reflecting m before bgc and sel were applied.
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
unapply.bgc.sel <- function(m, bgc=0, sel=0, alphabet="ACGT") {
  rphast.simplify.list(.Call.rphast("rph_tm_unapply_selection_bgc",
                                    as.matrix(m), alphabet, sel, bgc))
}


##' Add a lineage-specific model
##'
##' Lineage-specific models allow a different substitution model to be
##' defined on a specified set of branches.  An entirely different
##' substitution model can be used, as long as it is of the same order
##' and has the same number of states as the model used in the rest of the
##' tree.  Or, if the same substitution model is used, certain parameters
##' can be optimized separately from the main model, whereas others are
##' shared with the main model.
##'
##' A lineage-specific model is stored as a list with the following
##' elements: defn, rate.matrix, and optional elements backgd,
##' selection, bgc.
##'
##' defn is a character string which defines the model
##' in a way that phast can parse; it is a colon-delimited string with
##' 2 or 3 elements.  The first element indicates which branches the model
##' applies to, the second indicates which substitution model to use or
##' which parameters to optimize if the same substitution model is used
##' (and also may impose boundaries on these parameters).  The optional
##' third element is a list of parameters which will not be optimized by
##' phyloFit.
##'
##' backgd is the initial set of equilibrium frequencies for this model; if 
##' not present, then the equilibrium frequencies will be shared with the
##' main model.
##'
##' selection and bgc are optional parameters for the model with
##' biased gene conversion and selection.  If they are not provided this
##' model is not used.  Note that selection is defined relative to 
##' selection in the main model, if x$selection is not NULL (so the
##' total selection in the lineage-specific model is the sum of
##' the selection value in the main and lineage-specific model.
##'
##' A tree model can have multiple lineage-specific models; if a later
##' model applies to the same branch as an earlier model, then the later one
##' overrides it.
##'
##' All lineage-specific models are stored in the ls.model element
##' of the \code{tm} object.
##' @param x An object of type \code{tm}
##' @param branch If the lineage-specific model applies to a single branch,
##' it can be named here using the name of the node which descendes from the
##' branch.  See \code{\link{name.ancestors}} for naming internal nodes.
##' @param label (Alternative to branch).  The label which identifies the
##' branch(es) which this lineage-specific model should apply to.  Labels
##' are denoted in a tree with a pound sign and label following the node.
##' See \code{\link{label.branches}} and \code{\link{label.subtree}} to
##' add a label to a tree.
##' @param category An integer indicating which category/categories to
##' apply the lineage-specific model.  This only works if x$nratecats > 1.
##' A value of 0 or NULL implies all categories.  Otherwise this can be an
##' integer (or vector of integers) from 1..x$nratecats.
##' @param subst.mod A character string indicating the substitution model to
##' be used for the lineage-specific model.  If \code{NULL}, use the same
##' model as the rest of the tree.  See \code{\link{subst.mods}} for a
##' list of possible substituion models.
##' @param separate.params (Only applies if subst.mod is the same as main model)
##' A vector of character strings
##' indicating which parameters to estimate separately.  Possible values
##' are "kappa", "sel", "bgc", "gap_param", "backgd", and "ratematrix".
##' If backgd, selection, or bgc are provided as arguments, they are
##' automatically considered separate parameters and do not need to be
##' explicitly listed here.  "ratematrix" implies all
##' parameters describing the substitution model (but does not include
##' backgd, sel, or bgc).  Boundaries can be
##' optionally appended to parameter names with brackets, ie, "kappa[1,10]"
##' will set boundaries for kappa between 1 and 10 (see "Parameter
##' boundaries" section of \code{\link{phyloFit}}).
##' If subst.mod is different from the main model, then no parameters are
##' shared with main model.  However the equilibrium frequencies can be
##' shared by setting backgd to NULL.
##' @param const.params A character vector indicating which parameters to
##' hold constant at their initial values, rather than being optimized
##' upon a call to phyloFit.  Possible values are the same as for
##' separate.params, although no boundaries can be given here.
##' @param backgd The initial equilibrium frequencies to use for this
##' model.  If \code{NULL}, use the same as in the main model.
##' @param selection The selection parameter (from the sel+bgc model), 
##' relative to selection in the main model.
##' @param bgc The bgc parameter (from the sel+bgc model).
##' @return An object of type \code{tm}, identical to the input model but
##' with a new lineage-specific model added on.  This lineage-specific model
##' is not validated by this function.
##' @author Melissa J. Hubisz
##' @export
add.ls.mod <- function(x,
                       branch=NULL, label=NULL,
                       category=0,
                       subst.mod=NULL, separate.params=NULL,
                       const.params=NULL,
                       backgd=NULL,
                       selection=NULL,
                       bgc=NULL) {

  if (is.null(branch) && is.null(label))
    stop("either branch or label must be provided")
  if (! (is.null(branch) || is.null(label)))
    stop("only one of branch/label should be provided")
  check.arg(branch, "branch", "character", null.OK=TRUE)
  check.arg(label, "label", "character", null.OK=TRUE)
  # C code will check that branch/label exist on tree
  if (is.null(branch)) branchArg <- label
  else branchArg <- branch

  category <- check.arg(category, "category", "integer", null.OK=TRUE,
                        min.length=1L, max.length=NULL)
  if (!is.null(category)) {
    nrates <- x$nratecats
    if (is.null(nrates)) nrates <- 1
    if (sum(category < 0 | category > nrates) >= 1L)
      stop("category should be in range [0..,", nrates,"]")
    if (sum(category==0)==0) {
      catstr <- paste(category-1, sep=",")
      if (!is.null(branch)) {
        branch <- sprintf("%s#%s", branch, catstr)
      } else label <- sprintf("%s#%s", label, catstr)
    }
  }

  if (is.null(subst.mod)) subst.mod <- x[["subst.mod"]]
  if (subst.mod != x$subst.mod && !is.null(separate.params))
    stop("separate.params can only be non-NULL when the substitution model is the same as the main tree")
  if (!is.subst.mod.tm(subst.mod))
    stop("unknown substitution mod ", subst.mod)

  if (subst.mod == x[["subst.mod"]]) {
    if (!is.null(selection) && sum(const.params=="sel")==0) {
      if (is.null(separate.params))
        separate.params <- "sel"
      else if (length(grep("^sel(\\[.*,.*\\]+){0,1}$", separate.params)) == 0L)
        separate.params <- c(separate.params, "sel")
    }
    if (!is.null(bgc) && sum(const.params=="bgc")==0) {
      if (is.null(separate.params))
        separate.params <- "bgc"
      else if (length(grep("^bgc(\\[.*,.*\\]+){0,1}$", separate.params))==0L)
        separate.params <- c(separate.params, "bgc")
    }
  }

  defn <- sprintf("%s:%s", ifelse(!is.null(branch), branch, label),
                  ifelse(is.null(separate.params), subst.mod,
                         paste(separate.params, collapse=",")))
  if (!is.null(const.params)) {
    defn <- sprintf("%s:%s", defn, paste(const.params, sep=",", collapse=","))
  }
  origPtr <- as.pointer.tm(x)
  .Call.rphast("rph_tm_add_alt_mod", origPtr$externalPtr, defn)
  newmod <- from.pointer.tm(origPtr)

  # apply selection and/or bgc to ls model if provided
  if (! (is.null(bgc) && is.null(selection))) {
    if (is.null(newmod$ls.model$defn)) {
      nummod <- length(newmod$ls.model)
      lsmod <- newmod$ls.model[[nummod]]
    } else lsmod <- newmod$ls.model

    # newmod$bgc should always be null
    if (!is.null(newmod$selection)) {
      currMat <- unapply.bgc.sel(newmod$rate.matrix,
                                 alphabet=newmod$alphabet,
                                 sel=newmod$selection)
    } else currMat <- newmod$rate.matrix

    if (is.null(newmod$selection)) total.sel <- 0 else total.sel <- newmod$selection
    if (!is.null(selection)) total.sel <- total.sel + selection
    lsmod$rate.matrix <- apply.bgc.sel(currMat,
                                        alphabet=newmod$alphabet,
                                        bgc=if (is.null(bgc)) 0 else bgc,
                                        sel=total.sel)
    if (!is.null(selection)) lsmod$selection <- selection
    if (!is.null(bgc)) lsmod$bgc <- bgc

    if (is.null(newmod$ls.model$defn)) {
      newmod$ls.model[[nummod]] <- lsmod
    } else newmod$ls.model <- lsmod
  }
  newmod
}



##' Set the rate matrix of a tree model using model-specific parameters.
##'
##' The params argument is a numeric vector with the parameters specific to the
##' model being used.  Here is the meaning of params for each model:
##' \itemize{
##' \item{"JC69","F81": These models have no parameters; params should be NULL.}
##' \item{"K80","HKY85","HKY_CODON": params should be a single value representing the
##' transition-transversion ratio (kappa).}
##' \item{"HKY85+Gap": params should be a numeric vector of length 2; the first
##' element represents the transition/transversion ratio (kappa), and the second
##' is the "gap parameter", the factor by which substitution rates are
##' multiplied if they involve an indel event.}
##' \item{"REV": params should be a numeric vector of length 6 (assuming a
##' model with 4 states).  With n states the vector should be of length
##' n*(n-1)/2.  The first parameter applies to the entry in the 1st row, 2nd
##' column; the next to the 1st row, 3rd column, etc until the end of the
##' first row; the next parameter applies to the 2nd row, 3rd column; etc.}
##' \item{"SSREV": params should be a numeric vector of length 4. Assuming
##' an alphabet "ACGT", the first parameter is the substitution rate from A->C,
##' C->A, T->G, and G->T.  The second is the rate from A->G, G->A, T->C,
##' and C->T,  The third is the rate from A->T and T->A, and the last is
##' C->G and G->C.}
##' \item{"UNREST": params should be a numeric vector of length n*n-n (where
##' n is the number of states.  params fills in the rate matrix starting at
##' the first row going across, skipping diagonals.}
##' \item{"R2": Parameters should be a numeric vector of length 48.  There
##' are 16 states in order AA, AC, AG, AT, CA, ..., TT (assuming alphabet
##' order ACGT).  Parameters are filled in starting row 1, column 2, going
##' across and then down, filling in the matrix above the diagonal, and
##' reflecting into the matrix below the diagonal.  Only cells which
##' represent a substitution which requires exactly one mutation are filled
##' in; cells requiring greater than 1 mutation have rate 0.}
##' \item{"U2": Similar to R2, except parameters are a numeric vector of
##' length 96.  Parameters are filled in starting row 1, column to, going
##' across and then down, filling entire matrix (rather than refelcting
##' across the diagonal)}
##' \item{"R2S": Similar to R2, but with strand symmetry.  params should be a
##' vector of length 24.  Parameters are filled in a similar fashion as R2,
##' except the same parameter applies to substitutions which are strand
##' symmetric.}
##' \item{"U2S": Similar to U2, but with strand symmetry.  params should be a
##' vector of length 48.  Parameters are filled in a similar fashion as U2,
##' except the same parameter applies to substitutions which are strand
##' symmetric.}
##' \item{"R3","R3S","U3","U3S": Similar to R2, R2S, U2, and U2S, except there
##' are 64 states instead of 16.  params should be numeric vector of length
##' 288, 148, 576, or 288, for models R3, R3S, U3, and U3S respectively.}}
##' @param x An object of type \code{tm}.
##' @param params Parameters specific to the substitution model.  Should be a
##' numeric vector of length appropriate for the model.  See details below.
##' @param scale A logical value.  If \code{TRUE}, scale the matrix so that the
##' expected number of mutations per unit time is one per base pair.
##' @return An object of type \code{tm} with a rate matrix set according
##' to params.
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
set.rate.matrix.tm <- function(x, params=NULL, scale=TRUE) {
  x <- as.pointer.tm(x)
  .Call.rphast("rph_tree_model_set_matrix", x$externalPtr, params, scale)
  from.pointer.tm(x)
}

##' Get the parameters describing a rate matrix
##' @param x An object of type \code{tm}.
##' @return A numeric vector of parameters which can be used to describe
##' the transition matrix under the substitution model indicated in x.
##' May be NULL for certain models which have no parameters (JC69, F81).
##' The meaning of the parameters is described in
##' \code{\link{set.rate.matrix.tm}}.
##' @note The params returned may not describe the rate matrix passed
##' in, if the rate matrix does not follow the model indicated in x.
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
get.rate.matrix.params.tm <- function(x) {
  x <- as.pointer.tm(x)
  .Call.rphast("rph_tree_model_get_rate_matrix_params", x$externalPtr)
}


##' Make a bubble plot of a transition matrix
##' @param x A square matrix representing a continuous-time markov model;
##' rows should sum to zero, with negative values cross the diagonal
##' @param eq.freq A numeric vector giving the equilibrium frequencies
##' of each state.  If provided, the equilibrium frequencies will be
##' plotted along the bottom.
##' @param max.cex A scaling factor which determines the size of the
##' largest circle
##' @param eq.freq.max.cex A scaling factor which determines the size
##' of the largest circle in the equilibrium frequencies.
##' @param alphabet A character vector representing the state names for
##' each row/column of the matrix.  Can either be a vector of size
##' \code{nrow(m)} or a single character string with \code{nrow(m)}
##' characters.  Can also be \code{NULL} for no row/column labels.
##' @param col If \code{NULL}, all circles will be drawn in black.
##' Otherwise, col can be a matrix of the same dimension of \code{m},
##' each entry should indicate the color used for the corresponding
##' cell in the transition matrix.
##' @param eq.freq.col (Only applicable when eq.freq provided).  Should
##' be vector of same length as eq.freq, though values will be recycled.
##' Values in the vector indicate colors to draw the equilibrium frequency
##' bubbles.
##' @param filled If \code{TRUE}, plot filled circles.
##' @param add If \code{TRUE}, add to the existing plot.  Otherwise create
##' a new plot.
##' @param ... Further arguments to be passed to \code{plot}.
##' @author Melissa J. Hubisz
##' @method plot rate.matrix
##' @export plot.rate.matrix
##' @export
plot.rate.matrix <- function(x, eq.freq=NULL, max.cex=10.0,
                             eq.freq.max.cex=5.0,
                             alphabet=NULL, col=NULL,
                             eq.freq.col=NULL, filled=TRUE,
                             add=FALSE, ...) {
  check.arg(x, "x", "matrix", min.length=NULL, max.length=NULL)
  if (nrow(x) != ncol(x)) stop("x should be square matrix")
  m <- x
  check.arg(eq.freq, "eq.freq", min.length=nrow(m), max.length=nrow(m),
            null.OK=TRUE)
  if (!is.null(eq.freq)) {
    if (!is.null(eq.freq.col)) {
      eq.freq.col <- rep(eq.freq.col, length.out=length(eq.freq))
    } else eq.freq.col <- "black"
  }
  if (!is.null(col)) {
    if (nrow(col) != nrow(m) || ncol(col) != ncol(m))
      stop("col should have same dimensions as m")
    col <- as.numeric(col)
  } else {
    col <- "black"
  }
    
  check.arg(max.cex, "max.cex", "numeric", null.OK=FALSE)
  pch <- 21
  if (filled) {
    bg <- col
    bg.eq.freq <- eq.freq.col
  }else {
    bg <- NULL
    bg.eq.freq <- NULL
  }
  
  x <- sapply(1:ncol(m), rep, nrow(m))
  y <- matrix(rep(nrow(m):1, ncol(m)), nrow=nrow(m))
  diag(m) <- NA
  m <- sqrt(m)
  m <- m*max.cex/max(m, na.rm=TRUE)

  if (!add) {
    par(mar=c(2,4,4,2))
    xlim <- c(0.5, ncol(m)+0.5)
    ylim <- c(0.5, nrow(m)+0.5)
    if (!is.null(eq.freq)) ylim[1] <- -1
    plot(as.numeric(x), as.numeric(y), cex=as.numeric(m), xlab="", ylab="", xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, bty="n", col=col, bg=bg, pch=pch, ...)
  } else {
    points(as.numeric(x), as.numeric(y), cex=as.numeric(m), col=col, bg=bg, pch=pch, ...)
  }

  if (!is.null(alphabet)) {
    if (length(alphabet) == 1L && nchar(alphabet) == nrow(m)) {
      newalphabet <- character(nrow(m))
      for (i in 1:nrow(m)) newalphabet[i] <- substr(alphabet, i, i)
      alphabet <- newalphabet
    }
    if (length(alphabet) != nrow(m))
      stop("alphabet should be length 1 with nchar==nrow(m) or length nrow(m)")
    for (i in 1:length(alphabet)) {
      mtext(alphabet[i], side=c(3,2), at=c(i, nrow(m)-i+1), las=1, cex=1.2)
    }
  }
  
  if (!is.null(eq.freq)) {
    eq.freq <- sqrt(eq.freq)
    eq.freq <- eq.freq*eq.freq.max.cex/max(eq.freq, na.rm=TRUE)
    if (is.null(eq.freq.col)) eq.freq.col <- "black"
    points(1:nrow(m), rep(-0.5, length.out=nrow(m)),
           cex=eq.freq, col=eq.freq.col, bg=bg.eq.freq, pch=pch)
    mtext("eq.freq", side=2, at=-0.5, las=1)
  }
  
  invisible(NULL)
}

##' Make a bubble plot of the transition matrix for a tree model.
##' @param x An object of type \code{tm}.
##' @param show.eq.freq If \code{TRUE}, show bubbles representing equilibrium
##' frequencies along the bottom of plot.
##' @param max.cex A scaling factor which determines the size of the
##' largest circle
##' @param eq.freq.max.cex A scaling factor which determines the size
##' of the largest circle in the equilibrium frequencies.
##' @param alphabet A character vector representing the state names for
##' each row/column of the matrix.  Can either be a vector of size
##' \code{nrow(m)} or a single character string with \code{nrow(m)}
##' characters.  Can also be \code{NULL} for no row/column labels.
##' @param col If \code{NULL}, all circles will be drawn in black.
##' Otherwise, col can be a matrix of the same dimension of \code{m},
##' each entry should indicate the color used for the corresponding
##' cell in the transition matrix.
##' @param eq.freq.col Should
##' be vector of same length as eq.freq, though values will be recycled.
##' Values in the vector indicate colors to draw the equilibrium frequency
##' bubbles.
##' @param filled If \code{TRUE}, plot filled circles.
##' @param add If \code{TRUE}, add to the existing plot.  Otherwise create
##' a new plot.
##' @param ... Further arguments to be passed to \code{plot}.
##' @author Melissa J. Hubisz
##' @example inst/examples/plot-tm.R
##' @method plot tm
##' @export plot.tm
##' @export
plot.tm <- function(x, show.eq.freq=TRUE, max.cex=10.0,
                    eq.freq.max.cex=5.0,
                    alphabet=NULL, col=NULL,
                    eq.freq.col=NULL, filled=TRUE,
                    add=FALSE, ...) {
  plot.rate.matrix(x$rate.matrix,
                   if (show.eq.freq) eq.freq=x$backgd else NULL,
                   max.cex=max.cex,
                   eq.freq.max.cex=eq.freq.max.cex,
                   x[["alphabet"]],
                   col=col,
                   eq.freq.col=eq.freq.col,
                   filled=filled,
                   add=add, ...)
}

  
##' Make a bubble plot of a lineage-specific transition matrix of a
##' tree model.
##' @param x An object of type \code{tm}.
##' @param i An integer identifying which element of tm[["ls.model"]] to
##' plot.
##' @param show.eq.freq If \code{TRUE}, show bubbles representing equilibrium
##' frequencies along the bottom of plot.
##' @param max.cex A scaling factor which determines the size of the
##' largest circle
##' @param eq.freq.max.cex A scaling factor which determines the size
##' of the largest circle in the equilibrium frequencies.
##' @param alphabet A character vector representing the state names for
##' each row/column of the matrix.  Can either be a vector of size
##' \code{nrow(m)} or a single character string with \code{nrow(m)}
##' characters.  Can also be \code{NULL} for no row/column labels.
##' @param col If \code{NULL}, all circles will be drawn in black.
##' Otherwise, col can be a matrix of the same dimension of \code{m},
##' each entry should indicate the color used for the corresponding
##' cell in the transition matrix.
##' @param eq.freq.col Should
##' be vector of same length as eq.freq, though values will be recycled.
##' Values in the vector indicate colors to draw the equilibrium frequency
##' bubbles.
##' @param filled If \code{TRUE}, plot filled circles.
##' @param add If \code{TRUE}, add to the existing plot.  Otherwise create
##' a new plot.
##' @param ... Further arguments to be passed to \code{plot}.
##' @author Melissa J. Hubisz
##' @example inst/examples/plot-lsmodel-tm.R
##' @method plot lsmodel.tm
##' @export plot.lsmodel.tm
##' @export
plot.lsmodel.tm <- function(x, i=1, show.eq.freq=TRUE, max.cex=10.0,
                             eq.freq.max.cex=5.0,
                             alphabet=NULL, col=NULL,
                             eq.freq.col=NULL, filled=TRUE,
                             add=FALSE, ...) {
  if (!is.null(x$ls.model$defn)) {
    if (i != 1) stop("only one ls.model in x, but i=", i)
    lsmod <- x$ls.model
  } else {
    if (length(x$ls.model) < i)
    stop("not enough ls.models (%i)", length(x$ls.model))
    lsmod <- x$ls.model[[i]]
  }
  if (show.eq.freq) {
    if (is.null(lsmod$backgd)) {
      eq.freq <- x$backgd
    } else eq.freq <- lsmod$backgd
  } else eq.freq <- NULL
  plot.rate.matrix(lsmod$rate.matrix,
                   eq.freq,
                   max.cex=max.cex,
                   eq.freq.max.cex=eq.freq.max.cex,
                   x[["alphabet"]],
                   col=col,
                   eq.freq.col=eq.freq.col,
                   filled=filled,
                   add=add, ...)
}

##' Set up a tree model for branch site selection analysis
##' @param mod an object of type \code{tm}
##' @param foreground a character string giving a tree branch name or label
##' identifying foreground branches
##' @param bgc If \code{TRUE}, then use 8 categories of sites; four with bgc
##' in the foreground and four without.
##' @param altModel If \code{TRUE}, then optimize the foreground positive
##' selection parameter (constrained > 0).  Otherwise hold constant at 0.
##' @param init.sel.neg Initial value for negative selection parameter
##' @param init.sel.pos Initial value for positive selection paramter
##' @param init.bgc Initial value for bgc parameter (Ignored if bgc==FALSE)
##' @param init.weights Numeric vector of length three giving the initial
##' weight parameters.  The first two values determine the relative
##' frequencies of negatively, neutral, and positively selected sites.  The
##' last parameter determines the frequency of sites affected by bgc, and is
##' ignored if bgc==FALSE. All values should be >= 0.
##' @return An object of type \code{tm} which can be used as the init.mod
##' argument to phyloFit to perform the branch-site test.
##' @author Melissa J. Hubisz
##' @export
setup.branch.site.tm <- function(mod, foreground, bgc=FALSE, altModel=TRUE,
                                 init.sel.neg=0, init.sel.pos=0,
                                 init.bgc=0, init.weights=NULL) {
  mod <- as.pointer.tm(mod)
  if (!is.null(init.weights)) {
    if (length(init.weights) != 3L && bgc)
      stop("init.weights should be length 3 if bgc==TRUE")
    if (length(init.weights) != 3L && length(init.weights) != 2L && !bgc)
      stop("init.weights should be length 2 or 3 if bgc==FALSE")
    if (sum(init.weights < 0) > 0L)
      stop("all values in init.weights should be >=0")
  }
  .Call.rphast("rph_tm_setup_site_model", mod$externalPtr, foreground, bgc,
               altModel, init.sel.neg, init.sel.pos, init.bgc, init.weights)
  from.pointer.tm(mod)
}
