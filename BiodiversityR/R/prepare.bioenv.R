`prepare.bioenv` <-
function(env, as.numeric = c()) {
    env2 <- env
    if (is.null(as.numeric) == F) {
        for (i in 1:length(as.numeric)) {
            if(any(names(env) == as.numeric[i])) { 
                env2[, as.numeric[i]] <- as.numeric(env[, as.numeric[i]])
            }
        }
    }
    vars <- names(env2)
    for (i in 1:length(vars)) {
        focal.var <- which(names(env2)==vars[i])
        if (is.numeric(env2[, focal.var]) == F) {env2 <- env2[, -focal.var]}
    }
    return(env2)
}

