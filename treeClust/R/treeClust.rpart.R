treeClust.rpart <- function (i, dfx, d.num, control, rcontrol)
{
#
# Build a single tree from data "dfx" with column i (numeric) as the
# response.
#
    outlist <- list (DevRat = 0, Size = 1)
    if (length (unique (dfx[,i])) == 1)
        return (outlist)
#
# Check for NAs in the response. We will use this later.
#
    if (any (is.na (dfx[,i])))
        response.had.NAs <- TRUE
    else
        response.had.NAs <- FALSE
#
# "rcontrol" elements: use "information" for categorical responses.
#
    if (is.element (class (dfx[,i]), c("factor", "character")))
        rcontrol$parms <- list (split = "information")
#
# LOOP: Build the tree. If it has one leaf, drop this tree and quit.
# For some reason, the "-" notation in the formula doesn't seem to work
# consistently. So if we need to remove a column from the formula, we 
# remove it from the data.
#
    resp.name <- names (dfx)[i]
    rcontrol$formula <- eval (parse (text = paste (resp.name, "~ .")))
    cols.to.drop <- NULL
    while (1) {
        if (length (cols.to.drop) == 0) {
            rcontrol$data <- dfx
        } else {
            rcontrol$data <- dfx[,-cols.to.drop]
        }
        mytree <- do.call (rpart::rpart, rcontrol)
        if (nrow (mytree$frame) == 1)
             return (outlist)

#
# Extract the CP table. Compute the cutoff value based on the one-se
# rule (or whatever value serule has). Then find the smallest row whose
# xerror is smaller than that value. If that's row 1, drop this tree.
#
        cptbl <- mytree$cptable
        min.cp.dex <- which (cptbl[,"xerror"] == min(cptbl[,"xerror"]))[1]
        serule.value <-         cptbl[min.cp.dex,"xerror"] + 
               control$serule * cptbl[min.cp.dex,"xstd"]
        best.row <- min(which (cptbl[,"xerror"] <= serule.value))
        if(best.row == 1.)
            return (outlist)
#
# Prune to the CP value in this row, except so as to avoid rounding error,
# prune to something a little bigger -- say, halfway between the CP in
# this row and the one above.
#
        mytree <- rpart::prune.rpart (mytree, 
                  cp = (cptbl[best.row, "CP"] + cptbl[best.row - 1,"CP"])/2)
#
# The thing named "dev" really is the deviance for a regression tree,
# but not for a classification tree. Those we have to compute ourselves.
#
        if (is.factor (dfx[,i]))
            devs <- rp.deviance (mytree)
        else
            devs <- mytree$frame$dev
        orig.dev <- devs[1]
        new.dev <- sum (devs[mytree$frame$var == "<leaf>"])

        outlist$DevRat <- (orig.dev - new.dev)/orig.dev
        outlist$Size <- sum (mytree$frame[,"var"] == "<leaf>")
#
# If the DevRat does not exceed the threshold, we're done.
#
        if (outlist$DevRat <= control$DevRatThreshold)
            break
#
# Otherwise, remove from the set of predictors the one that was chosen at 
# the root. If we're removed them all, something weird is going on, but quit.
#
        cols.to.drop <- c(cols.to.drop, 
                          which (names (dfx) == mytree$frame[1,"var"]))
        outlist <- list (DevRat = 0, Size = 1) # reset this
        if (length (cols.to.drop) >= ncol (dfx) - 1) 
            return (outlist)

} # end of "while(1)" loop

#
# Save the leaf membership values. If this response had NAs, though, we
# will need to generate the full set, first.
#
    if (response.had.NAs) {
        mytree$where.orig <- mytree$where # for debugging
        mytree$where <- rpart.predict.leaves (mytree, dfx, type = "where")
# "where" without missings
        outlist$leaf.where <- factor (mytree$where)
    }
    else 
        outlist$leaf.where <- factor (mytree$where)
#
# Save the tree if asked.
# 
    if (control$keep.trees)
        outlist$tree <- mytree
    return (outlist)

}

