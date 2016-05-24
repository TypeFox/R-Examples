# Functions to create the variable name vector and degrees of freedom vector
# Includes a function to calculate the lev for each treatment/interaction
# Each function uses recursion in order to calculate the required vector

# Variable naming function
var.names <- function(num.ivs, ivs) {
    if(num.ivs == 0) {
        return()
    }
    if(length(ivs) == 2) {
        out.vars <- c(ivs, paste(ivs, collapse = "*"))
        return(out.vars)
    }
    else {
        all.variables <- sapply(1:length(ivs), FUN = function(x) {
            combn(ivs,x)
        }
        )
        interactions <- sapply(1:length(all.variables[[num.ivs]][1, ]), FUN = function(x) {
            c(paste(all.variables[[num.ivs]][,x], collapse = "*"))
        }
        )
        return(c(var.names(num.ivs - 1, ivs), interactions))
    }
}

# DF vector function
df.vector <- function(num.ivs, df.from.lev) {
    if(num.ivs == 0) {
        return()
    }
    if(length(df.from.lev) == 2) {
        out.df <- c(df.from.lev, prod(df.from.lev))
        return(out.df)
    }
    else {
        all.df <- sapply(1:length(df.from.lev), FUN = function(x) {
            combn(df.from.lev,x)
        }
        )
        interactions.df <- sapply(1:length(all.df[[num.ivs]][1,]), FUN = function(x) {
            c(prod(all.df[[num.ivs]][,x]))
        }
        )
        return(c(df.vector(num.ivs - 1, df.from.lev), interactions.df))
    }
}

# lev vector
lev.vector <- function(num.ivs, effect.lev) {
    if(num.ivs == 0) {
        return()
    }
    if(length(effect.lev) == 2) {
        out.df <- c(effect.lev, prod(effect.lev))
        return(out.df)
    }
    else {
        all.lev <- sapply(1:length(effect.lev), FUN = function(x) {
            combn(effect.lev,x)
        }
        )
        interactions.lev <- sapply(1:length(all.lev[[num.ivs]][1,]), FUN = function(x) {
            c(prod(all.lev[[num.ivs]][,x]))
        }
        )
        return(c(lev.vector(num.ivs - 1, effect.lev), interactions.lev))
    }
}
