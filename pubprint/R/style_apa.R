#############################################################################
# style_apa.R
#############################################################################

#' @include style_default.R
#' @include utils.R
NULL

#########
# initialise function
#########
style.apa.init <- function()
{
    return(list(character = style.default.character,
                numeric = style.default.numeric,
                p.value = style.apa.p.value,
                df = style.apa.df,
                t.test = style.apa.t.test,
                fisher = style.apa.fisher,
                pearson = style.apa.pearson,
                chisquared = style.apa.chisquared,
                shapiro = style.apa.shapiro,
                ks.test = style.apa.ks.test,
                summary.aovlist = style.apa.summary.aovlist,
                summary.lm = style.apa.summary.lm))
}

#########
# working functions
#########

# p value function
# x: all items will be returned as p value
# return: character vector with appropriate items
style.apa.p.value <- function(x)
{
    ret <- c(out.value(unlist(x),
                       name = "p",
                       nsmall = 3,
                       replace0 = TRUE,
                       leading0 = FALSE))

    return(ret)
}

# df function
# x: all items will be returned as df
# return: character vector with appropriate items
style.apa.df <- function(x)
{
    ret <- c(out.value(unlist(x),
                       drop0trailing = TRUE))

    return(ret)
}

# t-Test function
# x: first item is t.test (htest class) and second item a option d value (numeric)
# return: character vector with appropriate items
style.apa.t.test <- function(x,
                             print.d = TRUE,
                             print.estimate = TRUE,
                             estimate.names = names(x[[1]]$estimate))
{
    ret <- c()

    # estimates
    if (print.estimate)
        ret <- out.value(unname(x[[1]]$estimate), name = out.names(estimate.names))

    # t value
    ret <- c(ret,
             out.value(x[[1]]$statistic[[1]],
                       name = out.names(names(x[[1]]$statistic)),
                       inbracket = style.apa.df(x[[1]]$parameter)))

    # p value
    ret <- c(ret,
             style.apa.p.value(x[[1]]$p.value))

    # optional d value
    if (print.d && 1 < length(x))
        ret <- c(ret,
                 out.value(x[[2]],
                           name = "d"))

    return(ret)
}

# fisher.test function
# x: first item is fisher test (htest)
# return: character vector with appropriate items
style.apa.fisher <- function(x,
                             print.estimate = TRUE,
                             estimate.names = names(x[[1]]$estimate))
{
    ret <- c()

    # estimates
    if (print.estimate && !is.null(x[[1]]$estimate))
        ret <- out.value(unname(x[[1]]$estimate), name = out.names(estimate.names))

    # p value
    ret <- c(ret,
             style.apa.p.value(x[[1]]$p.value))

    return(ret)
}

# Pearson product-moment correlation
# x: first item is cor.test (htest class)
# return: character vector with appropriate items
style.apa.pearson <- function(x)
{
    # estimates
    ret <- out.value(unname(x[[1]]$estimate), name = out.names(names(x[[1]]$estimate)))

    # t value
    ret <- c(ret,
             out.value(x[[1]]$statistic[[1]],
                       name = out.names(names(x[[1]]$statistic)),
                       inbracket = style.apa.df(x[[1]]$parameter)))

    # p value
    ret <- c(ret,
             style.apa.p.value(x[[1]]$p.value))

    return(ret)
}

# chi squared test function
# x: first item is chi.squared (htest class)
# return: character vector with appropriate items
style.apa.chisquared <- function(x,
                                 chi.n.name = "N")
{
    ret <- c()
    myinb <- out.concat(out.value(x[[1]]$parameter,
                                  nsmall = 0),
                        out.value(sum(x[[1]]$observed),
                                  name = chi.n.name,
                                  nsmall = 0))

    # chi value
    ret <- c(ret,
             out.value(x[[1]]$statistic[[1]],
                       name = out.names(names(x[[1]]$statistic)),
                       inbracket = myinb))

    # p value
    ret <- c(ret,
             style.apa.p.value(x[[1]]$p.value))


    return(ret)
}

# Shapiro wilk test
# x: first item is shapiro.test (htest class)
# return: character vector with appropriate items
style.apa.shapiro <- function(x)
{
    ret <- c()

    # W value
    ret <- c(ret,
             out.value(x[[1]]$statistic[[1]],
                       name = out.names(names(x[[1]]$statistic))))

    # p value
    ret <- c(ret,
             style.apa.p.value(x[[1]]$p.value))

    return(ret)
}

# bartlett/friedman test function
# x: first item is bartlett or friedman test (htest class)
# return: character vector with appropriate items
style.apa.bartlett <- function(x)
{
    ret <- c()

    # chi squared value
    ret <- c(ret,
             out.value(x[[1]]$statistic[[1]],
                       name = out.names(names(x[[1]]$statistic)),
                       inbracket = x[[1]]$parameter))

    # p value
    ret <- c(ret,
             style.apa.p.value(x[[1]]$p.value))

    return(ret)
}

# Kolmogorow smirnow test
# x: first item is ks.test (htest class)
# return: character vector with appropriate items
style.apa.ks.test <- function(x)
{
    ret <- c()

    # D^- value
    ret <- c(ret,
             out.value(x[[1]]$statistic[[1]],
                       name = out.names(names(x[[1]]$statistic))))

    # p value
    ret <- c(ret,
             style.apa.p.value(x[[1]]$p.value))

    return(ret)
}

# ANOVA summary.aovlist function
# x: first item is summary.aovlist
# effect: string for which effect has to be returned
# return: character vector with appropriate items
style.apa.summary.aovlist <- function(x,
                                      effect = 1)
{
    ret <- c()

    # F value
    ret <- c(ret,
             out.value(x[[1]][[1]][[4]][effect], # 4 is F value
                       name = "F",
                       inbracket = out.concat(c(x[[1]][[1]][[1]][effect],
                                                tail(x[[1]][[1]][[1]], 1)))))

    # p value
    ret <- c(ret,
             style.apa.p.value(x[[1]][[1]][[5]][effect]))

    return(ret)
}

# regression analysis (summary.lm function)
# x: first item is summary.lm
# result: either "model", "equation" or a numeric
# r.squared: return normal or adjusted r.squared?
# respname: response name
# coeff.name: coefficient names (first character is name of intercept)
# return: character vector with appropriate items
style.apa.summary.lm <- function(x,
                                 result = "model",
                                 r.squared = c("normal", "adjusted"),
                                 respname = out.above("y", "^"),
                                 coeff.name = names(x[[1]]$coefficients[, 1]))
{
    r.squared = match.arg(r.squared)

    ret <- c()

    if ("model" == result)
    {
        # R squared
        ret <- c(ret,
                 out.value(switch(EXPR = r.squared,
                                  "normal" = x[[1]]$r.squared,
                                  "adjusted" = x[[1]]$adj.r.squared),
                           name = paste0("R", out.superscript("2"))))
        
        # F value
        ret <- c(ret,
                 out.value(x[[1]]$fstatistic[[1]],
                           name = "F",
                           inbracket = out.concat(x[[1]]$fstatistic[[2]],
                                                  x[[1]]$fstatistic[[3]])))

        # p value
        ret <- c(ret,
                 style.apa.p.value(pf(x[[1]]$fstatistic[[1]],
                                      x[[1]]$fstatistic[[2]],
                                      x[[1]]$fstatistic[[3]],
                                      lower.tail = FALSE)))
    }
    else if ("equation" == result)
    {
        ret <- c(ret,
                 out.value(x[[1]]$coefficients[1, 1]))

        # TODO: support here return of beta values as well?
        for (i in 2:length(x[[1]]$coefficients[, 1]))
            ret <- c(ret,
                     out.concat(out.value(x[[1]]$coefficients[i, 1]),
                                name = coeff.name[i],
                                sep = "*"))

        ret <- paste(ifelse(ret < 0,
                            ret,
                            paste0("+", ret)),
                     collapse = "")

        ret <- paste(respname,
                     "=",
                     ret,
                     collapse = "")
    }
    else if ("numeric" == class(result))
    {
        # b value
        # TODO: maybe beta would be better?
        ret <- c(ret,
                 out.value(x[[1]]$coefficients[result, 1],
                           name = "b"))

        # t vailue
        ret <- c(ret,
                 out.value(x[[1]]$coefficients[result, 3],
                           name = "t",
                           inbracket = x[[1]]$df[2]))

        # p value
        ret <- c(ret,
                 style.apa.p.value(x[[1]]$coefficients[result, 4]))
    }
    else
        stop("Argument \"result\" must either be \"model\", \"equation\" or a numeric.")

    return(ret)
}
