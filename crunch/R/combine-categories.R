##' Combine categories or responses
##'
##' @param variable Categorical, Categorical Array, or Multiple Response
##' variable
##' @param combinations list of named lists containing (1) "categories":
##' category ids or names for categorical types, or for multiple response,
##' "responses": subvariable names, aliases, or positional indices; (2) a
##' "name" for the new category or response; and (3) optionally, other category
##' ("missing", "numeric_value") or subvariable ("alias", "description")
##' attributes. If \code{combinations} is omitted, the resulting variable will
##' essentially be a copy (but see \code{link{copy}} for a more natural way to
##' do that, if desired).
##' @param ... Additional variable metadata for the new derived variable
##' @return A \code{\link{VariableDefinition}} that will create the new
##' comined-category or -response derived variable. Categories/responses not
##' referenced in \code{combinations} will be appended to the end of the
##' combinations.
##' @examples
##' \dontrun{
##' ds$fav_pet2 <- combine(ds$fav_pet, name="Pets (combined)",
##'     combinations=list(list(name="Mammals", categories=c("Cat", "Dog")),
##'                       list(name="Reptiles", categories=c("Snake", "Lizard"))))
##' ds$pets_owned2 <- combine(ds$allpets, name="Pets owned (collapsed)",
##'     combinations=list(list(name="Mammals", responses=c("Cat", "Dog"))))
##' }
##' @export
combine <- function (variable, combinations=list(), ...) {
    ## Validate inputs
    if (!(type(variable) %in% c("categorical", "categorical_array", "multiple_response"))) {
        halt("Cannot combine ", dQuote(name(variable)), ": must be type ",
            "categorical, categorical_array, or multiple_response")
    }
    if (!is.list(combinations) || !all(vapply(combinations, is.list, logical(1)))) {
        halt("'combinations' must be a list of combination specifications")
    }

    ## Get basic variable metadata
    newvar <- copyVariableReferences(variable)
    newvar$alias <- NULL ## Let server specify, or specify in ..., or on <-
    newvar <- updateList(newvar, list(...))
    newvar$type <- NULL ## Type is function of the derivation

    ## Construct expr
    if (type(variable) == "multiple_response") {
        combs <- combResps(subvariables(variable), combinations)
        newvar$expr <- zfunc("combine_responses",
            zcl(variable), list(value=combs))
        ## Give default name based on number of responses
        if (identical(newvar$name, name(variable))) {
            nvalidresps <- length(newvar$expr$args[[2]]$value)
            newvar$name <- paste0(newvar$name, " (", nvalidresps,
                ifelse(nvalidresps == 1, " response)", " responses)"))
        }
    } else {
        combs <- combCats(categories(variable), combinations)
        newvar$expr <- zfunc("combine_categories",
            zcl(variable), list(value=combs))
        ## Give default name based on number of categories
        if (identical(newvar$name, name(variable))) {
            nvalidcats <- length(Filter(Negate(function (x) isTRUE(x$missing)),
                newvar$expr$args[[2]]$value))
            newvar$name <- paste0(newvar$name, " (", nvalidcats,
                ifelse(nvalidcats == 1, " category)", " categories)"))
        }
    }
    class(newvar) <- "VariableDefinition"
    return(newvar)
}

combCats <- function (cats, combs) {
    ## Validate combinations
    if (!all(vapply(combs,
        function (x) all(c("name", "categories") %in% names(x)),
        logical(1)))) {

        halt("'combinations' must be a list of combination specifications. ",
            "See '?combine'.")
    }

    defaultCat <- list(missing=FALSE, numeric_value=NULL)
    ## Convert category names to ids
    ## Update each comb with default
    combs <- lapply(combs, function (x) {
        if (is.character(x$categories)) {
            x$categories <- n2i(x$categories, cats)
        }
        if (!is.numeric(x$categories)) {
            halt("Combinations must reference 'categories' by name or id")
        }
        x$combined_ids <- I(x$categories)
        x$categories <- NULL
        return(updateList(defaultCat, x))
    })

    ## Validate that they're all unique and nonmissing
    idsToCombine <- unlist(lapply(combs, vget("combined_ids")))
    badids <- setdiff(idsToCombine, ids(cats))
    if (length(badids)) {
        badnames <- vapply(
            Filter(function (x) any(x$combined_ids %in% badids), combs),
            vget("name"),
            character(1))
        halt(ifelse(length(badnames) == 1, "Combination ", "Combinations "),
            serialPaste(dQuote(badnames)),
            ifelse(length(badnames) == 1, " references", " reference"),
            ifelse(length(badids) == 1,
                " category with id ", " categories with ids "),
            serialPaste(badids),
            ifelse(length(badids) == 1,
                ", which does not exist", ", which do not exist"))
    }

    dupids <- duplicated(idsToCombine)
    if (any(dupids)) {
        dupnames <- i2n(idsToCombine[dupids], cats)
        halt(ifelse(length(dupnames) == 1, "Category ", "Categories "),
            serialPaste(dQuote(dupnames)),
            " referenced in multiple combinations")
    }

    ## Give valid ids to new combinations
    usedIds <- setdiff(ids(cats), idsToCombine)
    newIds <- setdiff(seq_len(length(usedIds) + length(combs)),
        usedIds)[seq_len(length(combs))]
    combs <- mapply(function (comb, i) {
        comb$id <- i
        return(comb)
    }, combs, newIds, SIMPLIFY=FALSE, USE.NAMES=FALSE)

    ## Append unreferenced cats to end
    oldCats <- lapply(cats@.Data[ids(cats) %in% usedIds], function (x) {
        x$combined_ids <- I(id(x))
        return(x)
    })
    combs <- c(combs, oldCats)

    ## One more validation
    newnames <- vapply(combs, vget("name"), character(1))
    dupnames <- duplicated(newnames)
    if (any(dupnames)) {
        halt("Duplicate category name given: ",
            serialPaste(dQuote(unique(newnames[dupnames]))))
    }
    return(combs)
}

combResps <- function (subvars, combs) {
    ## Validate combinations
    if (!all(vapply(combs,
        function (x) all(c("name", "responses") %in% names(x)),
        logical(1)))) {

        halt("'combinations' must be a list of combination specifications. ",
            "See '?combine'.")
    }

    ## Convert response names/aliases to urls
    subnames <- names(subvars)
    subaliases <- aliases(subvars)
    suburls <- urls(subvars)
    combs <- lapply(combs, function (x) {
        if (!is.character(x$responses)) {
            halt("Combinations must reference 'responses' by name or alias")
        }
        matches <- match(x$responses, subnames)
        if (any(is.na(matches))) {
            ## Try aliases instead
            matches <- match(x$responses, subaliases)
        }
        if (any(is.na(matches))) {
            halt("Response ", dQuote(x$name),
                " does not reference valid subvariables")
        }
        x$combined_ids <- I(suburls[matches])
        x$responses <- NULL
        return(x)
    })

    ## Append unreferenced subvars to end
    subvarsToCombine <- unlist(lapply(combs, vget("combined_ids")))
    oldSubvars <- lapply(setdiff(suburls, subvarsToCombine),
        function (u) {
            return(list(name=index(subvars)[[u]]$name,
                combined_ids=I(u)))
        })
    combs <- c(combs, oldSubvars)

    ## One more validation
    newnames <- vapply(combs, vget("name"), character(1))
    dupnames <- duplicated(newnames)
    if (any(dupnames)) {
        halt("Duplicate response name given: ",
            serialPaste(dQuote(unique(newnames[dupnames]))))
    }

    return(combs)
}
