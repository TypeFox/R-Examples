#' Create AGS from a mutation matrix
#'
#' Imports a TAB-delimited file with mutations.
#' This function creates a new list of AGSs from a table listing point (or otherwise defined quantitative) mutations. Such a matrix M typically has size Ngenes x Nsamples, so that the current function returns a list of \code{length=ncol(M)}. AGSs for each of the Nsamples are created as simple lists of all mutated genes G in a given sample S, i.e. any value X in the matrix M that satisfies condition \code{!is.na(X)} would be treated as a mutation. Eventual mutation types / categories are ignored. Wild type positions in the original TAB-delimited file should be filled with NAs.

#' @param MUT Matrix of size Ngenes x Nsamples (the both Ns can be arbitrarily defined, depending on the screen scale).
#' @param col.mask To include only columns with IDs that contain the specified mask. This parameter is aware of regular expression syntax, i.e. uses \code{grep(..., fixed = FALSE)}.
#' @param namesFromColumn Number of the column (if any) that contains the gene/protein names. Note that it is only necessary if the latter are NOT the unique rownames of the matrix. This could be sometimes useful for processing redundant gene  profiles with one-to-many mapping etc. Otherwise (i.e. the default), rownames shall contain gene IDs.
#' @param Lowercase render gene/protein IDs lower-case (Default:1)

#'
#' @examples
#' data("tcga.gbm",package="NEArender")
#' dim(tcga.gbm)
#' ags.list <- mutations2ags(tcga.gbm, col.mask="[-.]01$")
#' length(ags.list)
#' length(unique(unlist(ags.list)))
#' @export



mutations2ags  <- function(MUT, col.mask=NA, namesFromColumn=NA, Lowercase = 1) {
if (is.null(MUT)) {stop("Not enough parameters...");}
mgs.list <- NULL;

if (is.na(namesFromColumn)) {
m1 <- MUT;
} else {
m1 <- MUT[,(namesFromColumn+1):ncol(MUT)];
}
if (!is.na(col.mask)) {m1 <- m1[,colnames(m1)[grep(col.mask,colnames(m1))]];}
mgs.list <- apply(m1, 2, function (x) tolower(names(x))[which(!is.na(x) )]);
return(mgs.list);
}
