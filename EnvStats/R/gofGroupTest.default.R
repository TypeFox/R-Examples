gofGroupTest.default <-
function (object, group, test = "sw", distribution = "norm", 
    est.arg.list = NULL, n.classes = NULL, cut.points = NULL, 
    param.list = NULL, estimate.params = ifelse(is.null(param.list), 
        TRUE, FALSE), n.param.est = NULL, correct = NULL, digits = .Options$digits, 
    exact = NULL, ws.method = "normal scores", data.name = NULL, 
    group.name = NULL, parent.of.data = NULL, subset.expression = NULL, 
    ...) 
{
    if (is.null(data.name)) 
        data.name <- deparse(substitute(object))
    if (is.null(group.name)) 
        group.name <- deparse(substitute(group))
    y <- as.vector(unlist(object))
    if (!is.numeric(y)) 
        stop("All elements of 'object' must be numeric")
    if (!is.factor(group)) 
        group <- factor(unlist(group))
    if (any(!is.finite(group))) 
        stop("NA's/Inf's not allowed in 'group'.")
    if (length(y) != length(group)) 
        stop("'y' and 'group' must have the same length.")
    y.list <- split(y, group)
    n.grps <- length(y.list)
    if (n.grps == 1) {
        warning("Only one group supplied, so the function 'gofTest' was called.")
        ret.list <- gofTest(y = unlist(y.list), test = test, 
            distribution = distribution, est.arg.list = est.arg.list, 
            n.classes = n.classes, cut.points = cut.points, param.list = param.list, 
            estimate.params = estimate.params, n.param.est = n.param.est, 
            correct = correct, digits = digits, exact = exact, 
            ws.method = ws.method, data.name = data.name)
    }
    else {
        c.names <- names(y.list)
        if (identical(c.names, as.character(1:n.grps))) {
            c.names <- paste(group.name, c.names, sep = ".")
        }
        dum.list <- lapply(y.list, gofTest, test = test, distribution = distribution, 
            est.arg.list = est.arg.list, n.classes = n.classes, 
            cut.points = cut.points, param.list = param.list, 
            estimate.params = estimate.params, n.param.est = n.param.est, 
            correct = correct, digits = digits, exact = exact, 
            ws.method = ws.method)
        distribution <- dum.list[[1]]$distribution
        dist.abb <- dum.list[[1]]$dist.abb
        sample.size <- sapply(dum.list, function(x) x$sample.size)
        names(sample.size) <- c.names
        p.value.vec <- sapply(dum.list, function(x) x$p.value)
        names(p.value.vec) <- c.names
        bad.obs.vec <- sapply(dum.list, function(x) x$bad.obs)
        names(bad.obs.vec) <- c.names
        group.gof.list <- gofTest(p.value.vec, test = "ws", ws.method = ws.method)
        group.statistic <- group.gof.list$statistic
        group.sample.size <- group.gof.list$sample.size
        group.parameters <- group.gof.list$parameters
        group.p.value <- group.gof.list$p.value
        group.alternative <- paste("At least one group", "does not come from a", 
            paste(distribution, "Distribution."), sep = paste("\n", 
                space(33), sep = ""))
        group.method <- group.gof.list$method
        group.bad.obs <- group.gof.list$bad.obs
        names(group.bad.obs) <- "Group.bad.obs"
        group.scores <- group.gof.list$scores
        ret.list <- list(distribution = distribution, dist.abb = dist.abb, 
            statistic = group.statistic, sample.size = sample.size, 
            parameters = group.parameters, p.value = c(p.value.vec, 
                group.p.value), alternative = group.alternative, 
            method = group.method, data.name = data.name, grouping.variable = group.name)
        if (!is.null(parent.of.data)) 
            ret.list$parent.of.data <- parent.of.data
        if (!is.null(subset.expression)) 
            ret.list$subset.expression <- subset.expression
        ret.list <- c(ret.list, list(bad.obs = c(bad.obs.vec, 
            group.bad.obs), n.groups = n.grps, group.names = c.names, 
            group.scores = group.scores))
        oldClass(ret.list) <- "gofGroup"
    }
    ret.list
}
