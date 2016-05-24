Diff <- function(player1, player2, formula = NULL, id = "..", data = NULL,
                 separate.ability = NULL, refcat = NULL, contrasts = NULL,
                 subset = NULL) {
    player.one <- player1[[id]]
    player.two <- player2[[id]]

    if (!is.factor(player.one) || !is.factor(player.two) ||
        !identical(levels(player.one), levels(player.two)))
        stop("'player1$", id, "' and 'player2$", id,
             "' must be factors with the same levels")
    if (!identical(attr(player.one, "contrasts"),
                   attr(player.two, "contrasts")))
        stop("'player1$", id, "' and 'player2$", id,
             "' must have the same contrasts attribute")
    if(is.null(formula)) formula <- reformulate(id)

    players <- levels(player.one)
    nplayers <- nlevels(player.one)
    ncontests <- length(player.one)
    D <- matrix(nrow = ncontests, ncol = nplayers)
    D <- col(D) == as.numeric(player.one)
    D <- D - (col(D) == as.numeric(player.two))
    colnames(D) <- paste(id, players, sep = "")

    fixed <- nobars(formula)
    X <- offset <- missing <- term.labels <- NULL
    saturated <- FALSE
    sep <- list()
    empty <- is.null(fixed) || is.empty.model(mt <- terms(fixed))
    if (!empty) {
        factors <- attr(mt, "factors")
        term.labels <- as.character(colnames(factors))
        vars <- rownames(factors)
        indexed <- grep("[[][^],]+[],]", vars)
        if (length(indexed)) { #set NAs to zero
            indices <- gsub("[^[]*[[]([^],]+)[],].*", "\\1", vars[indexed])
            vars <-  gsub("[[][^]]*[]]", "", vars[indexed])
            ## assumes no overlap, e.g. no age[..]:judge.gender[judge]
            grp <- split(vars, indices)
            for (ind in names(grp)) {
                vars <- model.frame(terms(reformulate(grp[[ind]])),
                                    data = data, na.action = na.pass)
                lev <- levels(eval(as.name(ind), c(player1, data)))
                as.sep <- rowSums(is.na(vars)) | lev %in% separate.ability
                if (any(as.sep)) {
                    sep[[ind]] <- as.sep
                    vars[sep[[ind]], ] <- lapply(vars, function(x)
                                                 max(levels(x)[1], 0))
                    colnames(vars) <- gsub(".*[$[],? ?\"?([^]\"]*).*", "\\1",
                                           grp[[ind]])
                    labels <- gsub("([^[$]*)[[$].*", "\\1", grp[[ind]])
                    for (lab in intersect(labels, grp[[ind]]))
                        data[lab] <- vars[lab]
                    for (lab in setdiff(labels, grp[[ind]]))
                        data[[lab]] <- vars[, labels == lab, drop = FALSE]
                }
            }
            if (length(sep)) {
                fixed <- reformulate(c(names(sep), attr(mt, "term.labels"),
                                       rownames(attr(mt, "factors"))[attr(mt, "offset")]))
                mt <- terms(fixed)
            }
        }

        idterm <- id %in% rownames(attr(mt, "factors"))
        mf1 <- model.frame(mt, data = c(player1, data), na.action = na.pass)
        if (nrow(mf1) != ncontests)
            stop("Predictor variables are not of the correct length --",
                 "they probably need indexing in 'formula'.")
        mf2 <- model.frame(mt, data = c(player2, data), na.action = na.pass)
        if (idterm){
            if (!is.null(refcat)) {
                mf1[[id]] <- relevel(mf1[[id]], refcat)
                mf2[[id]] <- relevel(mf2[[id]], refcat)
                if (!missing(contrasts)) contrasts[[id]] <- "contr.treatment"
            } else {
                ## 'else' defined by contrasts arg/contrasts attr of id factor
                ## leave refcat NULL
                if (is.null(contrasts))
                    contrasts[[id]] <- attr(player.one, "contrasts")
            }
        }
        offset <- model.offset(mf1)
        if (!is.null(offset)) offset <- offset - model.offset(mf2)

        if (length(sep)){ #create separate effect factor
            recode <- function(x, keep){
                lev <- levels(x)
                ext <- make.unique(c(lev[keep], "nosep"))[sum(keep) + 1]
                levels(x)[!keep]  <- ext
                relevel(x, ref = ext)
            }
            for (ind in names(grp)) {
                mf1[ind] <- recode(mf1[[ind]], sep[[ind]])
                mf2[ind] <- recode(mf2[[ind]], sep[[ind]])
            }
        }

        X1 <- model.matrix(fixed, mf1, contrasts = contrasts)
        X2 <- model.matrix(fixed, mf2, contrasts = contrasts)
        X <- X1 - X2
        ## will need to check for saturation in each set of indexed var
        ## - however as only allowing (1|..) just consider player id for now

        saturated <- qr(na.omit(X))$rank == qr(na.omit(cbind(D, X)))$rank && !idterm
        if (all(X[,1] == 0)) X <- X[, -1, drop = FALSE]
        attr(X, "assign") <- attr(X1, "assign")[-1]
    }

    random <- findbars(formula[[2]])
    if (!is.null(random)) {
        if (!is.list(random)) random <- list(random)
        if (length(random) > 1 ||
            random[[1]] != parse(text = paste("1|", id, sep = ""))[[1]])
            stop("Currently '(1 | ", id, ")' is the only random effects",
                 "structure allowed.")
        random <- D
    }
    else if (!empty && (!idterm & !saturated))
        warning("Ability modelled by predictors but no random effects",
                call. = FALSE)

    if (length(sep)) {
        attr(X, "assign") <- attr(X, "assign") - 1
        if (!is.null(random))
            random <- D[,!sep[[id]], drop = FALSE]
    }
    list(X = X, random = random, offset = offset,
         term.labels = term.labels, refcat = refcat, contrasts = contrasts,
         saturated = saturated)
}

