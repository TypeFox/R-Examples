`callAllSol` <- function(k, chart) {
    return(.Call("allSol", k, chart, PACKAGE="QCA"))
}


`callFindSubsets` <- function(x, noflevels, mbase, maximum) {
    return(.Call("findSubsets", x, noflevels, mbase, maximum, PACKAGE="QCA"))
}


`callRemoveRedundants` <- function(expressions, noflevels, mbase) {
    return(.Call("removeRedundants", expressions, noflevels, mbase, PACKAGE="QCA"))
}


`callSolveChart` <- function(combos, chart) {
    return(.Call("solveChart", combos, chart, PACKAGE="QCA"))
}


`callSuperSubsetMem` <- function(conditions, noflevels, mbase, fuzzycc, outcome, relation) {
    return(.Call("superSubsetMem", conditions, noflevels, mbase, fuzzycc, outcome, relation, PACKAGE="QCA"))
}


`callSuperSubset` <- function(conditions, nk, fuzzycc, outcome, relation) {
    return(.Call("superSubset", conditions, nk, fuzzycc, outcome, relation, PACKAGE="QCA"))
}


`callTruthTable` <- function(conditions, tt, fuzzycc, outcome) {
    return(.Call("truthTable", conditions, tt, fuzzycc, outcome, PACKAGE="QCA"))
}






