#' @export
stemplot <-
function (x, ...) 
{
    UseMethod("stemplot", x)
}


#' @export
stemplot.data.frame <-
function (x, ...) 
{
    doIt <- which(sapply(x, is.numeric))
    for (column in doIt) {
        cat(paste(rep("-", 25), sep = ""))
        cat(paste("\n", "Stemplot of ", names(x)[column], "\n\n", 
            sep = ""))
        print(stemplot.numeric(x[, column], ...))
    }
}


#' @export
stemplot.numeric <-
function (x, unit, split, Min, Max, rule.line = c("Dixon", "Velleman", 
    "Sturges"), style = c("bare", "Tukey"), trim.outliers = TRUE, 
    depths = TRUE, reverse.negative.leaves = TRUE, ...) 
{
    data <- x
    rule.line <- match.arg(rule.line)
    style <- match.arg(style)
    n <- length(data <- sort(data))
    row.max <- floor(c(Dixon = 10 * log(n, 10), Velleman = 2 * 
        sqrt(n), Sturges = 1 + log(n, 2)))[rule.line]
    stats <- boxplot(data, plot = FALSE)
    if (missing(Min)) 
        Min <- if (trim.outliers) 
            stats$stats[1, 1]
        else min(data, na.rm = TRUE)
    if (missing(Max)) 
        Max <- if (trim.outliers) 
            stats$stats[5, 1]
        else max(data, na.rm = TRUE)
    if (Min == Max) {
        Min <- min(data, na.rm = TRUE)
        Max <- max(data, na.rm = TRUE)
    }
    spannweite.red <- Max - Min
    zeilen.intervall.laenge <- spannweite.red/row.max
    factor <- if (missing(unit)) 
        10^ceiling(log(zeilen.intervall.laenge, 10))
    else 10^round(log(unit * 10, 10))
    z <- zeilen.intervall.laenge/factor
    delta.tick <- c(0.2, 0.2, 0.5, 1)[sum(z > c(0, 0.1, 0.2, 
        0.5))]
    if (missing(split)) 
        split <- round(1/delta.tick)
    else delta.tick <- 1/split
    data.tr <- data/factor
    Min.tr <- Min/factor
    Max.tr <- Max/factor
    spannweite.red <- Max.tr - Min.tr
    sk.min <- floor(Min.tr)
    sk.max <- ceiling(Max.tr)
    skala <- seq(sk.min, sk.max, by = delta.tick)
    if (sk.min < 0) 
        skala <- c(sk.min - delta.tick, skala)
    if (sk.max < 0) 
        skala <- skala[-length(skala)]
    lo.limit <- if (trim.outliers) 
        skala[1]
    else -Inf
    lo.log <- if (skala[1] < 0) 
        data.tr <= lo.limit
    else data.tr < lo.limit
    n.sk <- length(skala)
    hi.limit <- if (trim.outliers) 
        skala[n.sk] + delta.tick
    else Inf
    hi.log <- if (skala[n.sk] >= 0) 
        data.tr >= hi.limit
    else data.tr > hi.limit
    n.lower.extr.values <- sum(lo.log)
    n.upper.extr.values <- sum(hi.log)
    if (0 < n.lower.extr.values) {
        lower.line <- paste("LO:", paste(data[lo.log], collapse = " "))
    }
    if (0 < n.upper.extr.values) {
        upper.line <- paste("HI:", paste(data[hi.log], collapse = " "))
    }
    data.tr.red <- data.tr[(!lo.log) & (!hi.log)]
    stem <- ifelse(data.tr.red < 0, ceiling(data.tr.red), floor(data.tr.red))
    leaf <- floor(abs(data.tr.red * 10 - stem * 10))
    class.of.data.tr <- unlist(c(sapply(data.tr.red[data.tr.red < 
        0], function(x, sk) length(sk) - sum(-sk <= -x), skala), 
        sapply(data.tr.red[data.tr.red >= 0], function(x, sk) sum(sk <= 
            x), skala)))
    class.of.data.tr <- c(1:length(skala), class.of.data.tr)
    class.negative <- skala < 0
    leaf.grouped <- split(c(rep(-1, length(skala)), leaf), class.of.data.tr)
    leaf.grouped <- lapply(leaf.grouped, function(x) {
        sort(x[-1])
    })
    if (reverse.negative.leaves) {
        for (i in seq(class.negative)) if (class.negative[i]) 
            leaf.grouped[[i]] <- rev(leaf.grouped[[i]])
    }
    leaf.grouped.ch <- paste("|", unlist(lapply(leaf.grouped, 
        paste, collapse = "")))
    class.neg.zero <- floor(skala) == -1
    line.names <- skala
    line.names[class.negative] <- line.names[class.negative] + 
        1
    line.names <- as.character(floor(line.names))
    line.names[class.neg.zero] <- "-0"
    if (style == "Tukey") {
        switch(as.character(split), `1` = {
        }, `2` = {
            h <- round(2 * (skala%%1))
            line.names <- paste(line.names, ifelse(skala < 0, 
                c(".", "*")[1 + h], c("*", ".")[1 + h]), sep = "")
        }, `5` = {
            h <- round(5 * (skala%%1))
            line.names[h > 0 & h < 4] <- ""
            line.names <- paste(line.names, ifelse(skala < 0, 
                c(".", "s", "f", "t", "*")[1 + h], c("*", "t", 
                  "f", "s", ".")[1 + h]), sep = "")
        })
    }
    ragged.left <- function(ch.lines) {
        max.n <- max(n.lines <- nchar(ch.lines))
        h <- paste(rep(" ", max.n), collapse = "")
        ch.lines <- paste(substring(h, 1, 1 + max.n - n.lines), 
            ch.lines)
        ch.lines
    }
    line.names <- ragged.left(line.names)
    n.class <- unlist(lapply(leaf.grouped, length))
    select <- (cumsum(n.class) > 0) & rev((cumsum(rev(n.class)) > 
        0))
    depth <- cumsum(n.class) + n.lower.extr.values
    depth.rev <- rev(cumsum(rev(n.class)) + n.upper.extr.values)
    uplow <- depth >= depth.rev
    pos.median <- which(uplow)[1] + (-1:0)
    h <- abs(depth[pos.median] - depth.rev[pos.median])
    pos.median <- pos.median[1] + (h[1] > h[2])
    depth[uplow] <- depth.rev[uplow]
    depth <- paste(depth, "")
    depth[pos.median] <- paste("(", n.class[pos.median], ")", 
        sep = "")
    depth[n.class == 0] <- " "
    depth <- if (depths) 
        ragged.left(depth)
    else ""
    info <- c(paste("1 | 2: represents", 1.2 * factor), paste(" leaf unit:", 
        factor/10), paste("         n:", n), "")
    stem <- paste(depth, line.names, leaf.grouped.ch)
    if ((style != "Tukey") || (split != 5) || (sum(select) > 
        4)) 
        stem <- stem[select]
    if (exists("lower.line")) 
        stem <- c(lower = lower.line, stem)
    if (exists("upper.line")) 
        stem <- c(stem, upper = upper.line)
    result <- list(info = info, stem = stem)
    class(result) <- "stemplot"
    result
}
