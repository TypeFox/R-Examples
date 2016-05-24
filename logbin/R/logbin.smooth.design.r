logbin.smooth.design <- function(interpret, allref, design.knots, design.param)
{
    p.env <- environment(interpret$full.formula)
    formula.lhs <- Reduce(paste, deparse(interpret$full.formula[[2L]]))
    formula.rhs <- Reduce(paste, deparse(interpret$full.formula[[3L]]))
    monotonic.new <- allref$monotonic
    data.new <- allref$data
    nobs <- nrow(data.new)
    smoothterms <- names(interpret$smooth.spec)
    knots <- list()
    for (smth in smoothterms) {
        smthlabel <- interpret$smooth.spec[[smth]]$termlabel
        smthtype <- class(interpret$smooth.spec[[smth]])
        if (attr(allref$terms,"dataClasses")[smth] != "numeric")
            stop(gettextf("error in %s: variable %s must be numeric",
                    smthlabel, smth), domain = NA)
        x <- allref$data[[smth]]
        if (smthtype == "Iso.smooth") {
            unique.x <- sort(unique(x[!is.na(x)]))
            knots[[smth]] <- unique.x
            x.new <- matrix(0, nrow = nobs, ncol = length(unique.x)-1)
            for (i in 2:length(unique.x))
                x.new[, i-1] <- as.numeric(x < unique.x[i])
            colnames(x.new) <- paste(smthlabel, 2L:length(unique.x), sep = "")
        } else if (smthtype == "B.smooth") {
            min.x <- min(x, na.rm = TRUE)
            max.x <- max(x, na.rm = TRUE)
            knots.in <- interpret$smooth.spec[[smth]]$knots
            if (!is.null(knots.in)) {
                if (any(knots.in <= min.x) || any(knots.in >= max.x))
                    stop(gettextf("error in %s: specified knots must be in the interior of the range of %s",
                        smthlabel, smth), domain = NA)
                knots[[smth]] <- c(rep(min.x,3), knots.in, rep(max.x,3))
            } else {
                num.knots <- as.numeric(design.knots[smth])
                unique.x <- sort(unique(x))
                if(num.knots > length(unique.x)-2)
                    stop(gettextf("error in %s: number of knots greater than number of unique values of %s",
                        smthlabel, smth), domain = NA)
                knots.int <- quantile(unique.x, probs = seq(0, 1, len = num.knots+2), names = FALSE)
                knots[[smth]] <- c(rep(min.x,2), knots.int, rep(max.x, 2))
            }
			B <- matrix(NA, nrow = nobs, ncol = length(knots[[smth]])-3)
            B[!is.na(x),] <- splines::splineDesign(knots[[smth]], x, ord = 3)
            colnames(B) <- paste(smthlabel, seq_len(ncol(B)), sep = "")
            ref <- allref$allref[[smth]][[as.numeric(design.param[smth])]]
            if (length(ref)==1)
                x.new <- B[, -ref, drop=FALSE] 
            else {
                x.temp <- B[, ref][, -ncol(B), drop=FALSE]
                x.new <- x.temp
                for (j in 2L:ncol(x.temp)) 
                    x.new[,j] <- rowSums(x.temp[, 1L:j])
            }
        }
        nvar.new <- NCOL(x.new)
        varnames.new <- colnames(x.new)
        data.new <- cbind(data.new, -x.new)
        formula.rhs <- gsub(smthlabel, paste(paste("`",varnames.new,"`",sep=""), 
									collapse=" + "), formula.rhs, fixed = TRUE)
        pos <- which(names(monotonic.new) == smth)
        mono.add <- rep(TRUE, nvar.new)
        names(mono.add) <- varnames.new
        monotonic.new <- append(monotonic.new[-pos], mono.add, after = pos-1)
    }
    formula.new <- as.formula(paste(formula.lhs,"~",formula.rhs), env = p.env)
    list(formula = formula.new, data = data.new, monotonic = monotonic.new, knots = knots)
}