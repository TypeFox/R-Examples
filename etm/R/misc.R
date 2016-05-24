### Some useful miscellaneous functions ###

tra_ill <- function(state.names = c("0", "1", "2")) {

    if (length(state.names) != 3)
        stop("An illness-death model has 3 states")
    
    tra <- matrix(FALSE, ncol = 3, nrow = 3,
                  dimnames = list(state.names, state.names))
    tra[1, 2:3] <- TRUE
    tra[2, 3] <- TRUE
    tra
}

tra_ill_comp <- function(nComp = 2,
                         state.names = as.character(seq(0, nComp + 1, 1))) {

    if (nComp == 1)
        stop("No competing risks. Use 'tra_ill' instead")

    nstates <- length(state.names)
    if (length(state.names) != nComp + 2)
        stop(paste("Something is wrong with 'state.names'. The specified multistate model has ",
                   nComp + 2L, " states", sep = "")) 
    
    tra <- matrix(FALSE, nstates, nstates,
                  dimnames = list(state.names, state.names))
    tra[1, 2:nstates] <- TRUE
    tra[2, 3:nstates] <- TRUE
    tra
}

tra_comp <- function(nComp = 2,
                     state.names = as.character(seq(0, nComp))) {

    if (nComp == 1)
        stop("That's not a competing risks model. Use 'tra_surv' instead")
    nstates <- length(state.names)
    if (nstates != nComp + 1L)
        stop(paste("Something is wrong with 'state.names'. The specified multistate model has ",
                   nComp + 1L, " states", sep = "")) 
    
    tra <- matrix(FALSE, nstates, nstates,
                  dimnames = list(state.names, state.names))
    tra[1, 2:nstates] <- TRUE
    tra
}

tra_surv <- function(state.names = c("0", "1")) {

    if (length(state.names) != 2)
        stop("Survival model has 2 states")
    
    tra <- matrix(FALSE, ncol = 2, nrow = 2,
                  dimnames = list(state.names, state.names))
    tra[1, 2] <- TRUE
    tra
}

