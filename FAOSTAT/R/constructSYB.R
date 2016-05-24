##' Construct/Creat new variable.
##'
##' A function used to construct new variables from existing variables.
##'
##' Currently two types of construction are supported, either share or
##' growth rate computation.
##'
##' Share can be a share of total or share of another variable depending
##' on whether an additional variable is supplied or not.
##'
##' @param data The data frame containing the raw variable
##' @param origVar1 The variable name to be used in construction, refer
##' to Details for more information and useage.
##' @param origVar2 The variable name to be used in construction, refer
##' to Details for more information and useage.
##' @param newVarName The name assigned to the new variable, if missing
##' then .SC/.SH/.GR/.CH will be appended depending on the type of
##' construction
##' @param constructType The type of construction, refer to Details
##' for more information.
##' @param grFreq The frequency for the growth rate to be computed.
##' @param grType The method for the growth to be calculated, currently
##' supports least squares and geometric.
##' @param baseYear The base year to be used for constructing index.
##' @export
##' @return A data frame containing both the original data frame and the
##' processed data and also a list indicating whether the construction
##' passed or failed.

constructSYB = function(data, origVar1, origVar2, newVarName = NA,
    constructType = c("share", "growth", "change", "index"),
    grFreq = 1, grType = c("ls", "geo"), baseYear = 2000){
    ## The length of the variables must be the same
    checkDim = try(data.frame(origVar1, origVar2, newVarName,
        constructType, grFreq, baseYear))
    if(inherits(checkDim, "try-error"))
        stop("length of original variables are not the same")
    n = length(origVar1)
    result = data.frame(newVar = newVarName,
        Success = logical(length(newVarName)),
        Reason = character(length(newVarName)),
        stringsAsFactors = FALSE)
    printLab(paste("Constructing new variables (", n, " in Total)", sep = ""))
    for(i in 1:n){
        cat(paste("(", i, "): ", sep = ""))
        if(origVar1[i] %in% colnames(data) &&
           (origVar2[i] %in% colnames(data) || is.na(origVar2[i]))){
            switch(constructType[i],
                   share = {tmp = try(shConstruct(data = data,
                                totVar = origVar2[i],
                                shareVar = origVar1[i],
                                newVarName = newVarName[i]), silent = TRUE)},
                   growth = {tmp = try(grConstruct(data = data,
                                 origVar = origVar1[i],
                                 newVarName = newVarName[i],
                                 type = grType[i], n = grFreq[i]),
                                 silent = TRUE)},
                   change = {tmp = try(chConstruct(data = data,
                                 origVar = origVar1[i],
                                 newVarName = newVarName[i],
                                 n = grFreq[i]), silent = TRUE)},
                   index = {tmp = try(indConstruct(data = data,
                                origVar = origVar1[i],
                                newVarName = newVarName[i],
                                baseYear = baseYear), silent = TRUE)}
                   )
            ## Sometimes nan and Inf are a result of divisible by zero, in
            ## this case we will replace them with missing value for unknown
            ## information. Since Inf will usually cause problem in
            ## computation.
            print(str(tmp))
            tmp[(is.nan(tmp[, 3]) | (tmp[, 3] == Inf)) &
                !is.na(tmp[, 3]), 3] = NA
            if(!(inherits(tmp, "try-error"))){
                cat(paste("PASS: ", newVarName[i],
                          " sucessfully constructed\n", sep = ""))
                result[i, "Success"] = TRUE
                result[i, "Reason"] = "Construction Successful"
                data = merge(data, tmp, by = c("FAOST_CODE", "Year"))
            } else {
                cat(paste("FAIL: ", newVarName[i],
                          ", Check error message\n", sep = ""))
                result[i, "Success"] = FALSE
                result[i, "Reason"] = attr(tmp, "condition")$message
            }
        } else {
            cat(paste("FAIL: ",  origVar1[i], " or ",
                      origVar2[i], " not found in data\n", sep = ""))
            result[i, "Success"] = FALSE
            result[i, "Reason"] = "Variable not found in data"
        }
    }
    cat(paste("\nNumber of variables successfully computed: ",
              sum(result$Success), " out of ", NROW(result), "\n", sep = ""))
    list(data = data, result = result)
}
