plot.nn <-
function (x, rep = NULL, x.entry = NULL, x.out = NULL, radius = 0.15, 
    arrow.length = 0.2, intercept = TRUE, intercept.factor = 0.4, 
    information = TRUE, information.pos = 0.1, col.entry.synapse = "black", 
    col.entry = "black", col.hidden = "black", col.hidden.synapse = "black", 
    col.out = "black", col.out.synapse = "black", col.intercept = "blue", 
    fontsize = 12, dimension = 6, show.weights = TRUE, file = NULL, 
    ...) 
{
    net <- x
    if (is.null(net$weights)) 
        stop("weights were not calculated")
    if (!is.null(file) && !is.character(file)) 
        stop("'file' must be a string")
    if (is.null(rep)) {
        for (i in 1:length(net$weights)) {
            if (!is.null(file)) 
                file.rep <- paste(file, ".", i, sep = "")
            else file.rep <- NULL
            dev.new()
            plot.nn(net, rep = i, x.entry, x.out, radius, arrow.length, 
                intercept, intercept.factor, information, information.pos, 
                col.entry.synapse, col.entry, col.hidden, col.hidden.synapse, 
                col.out, col.out.synapse, col.intercept, fontsize, 
                dimension, show.weights, file.rep, ...)
        }
    }
    else {
        if (is.character(file) && file.exists(file)) 
            stop(sprintf("%s already exists", sQuote(file)))
        result.matrix <- t(net$result.matrix)
        if (rep == "best") 
            rep <- as.integer(which.min(result.matrix[, "error"]))
        if (rep > length(net$weights)) 
            stop("'rep' does not exist")
        weights <- net$weights[[rep]]
        if (is.null(x.entry)) 
            x.entry <- 0.5 - (arrow.length/2) * length(weights)
        if (is.null(x.out)) 
            x.out <- 0.5 + (arrow.length/2) * length(weights)
        width <- max(x.out - x.entry + 0.2, 0.8) * 8
        radius <- radius/dimension
        entry.label <- net$model.list$variables
        out.label <- net$model.list$response
        neuron.count <- array(0, length(weights) + 1)
        neuron.count[1] <- nrow(weights[[1]]) - 1
        neuron.count[2] <- ncol(weights[[1]])
        x.position <- array(0, length(weights) + 1)
        x.position[1] <- x.entry
        x.position[length(weights) + 1] <- x.out
        if (length(weights) > 1) 
            for (i in 2:length(weights)) {
                neuron.count[i + 1] <- ncol(weights[[i]])
                x.position[i] <- x.entry + (i - 1) * (x.out - 
                  x.entry)/length(weights)
            }
        y.step <- 1/(neuron.count + 1)
        y.position <- array(0, length(weights) + 1)
        y.intercept <- 1 - 2 * radius
        information.pos <- min(min(y.step) - 0.1, 0.2)
        if (length(entry.label) != neuron.count[1]) {
            if (length(entry.label) < neuron.count[1]) {
                tmp <- NULL
                for (i in 1:(neuron.count[1] - length(entry.label))) {
                  tmp <- c(tmp, "no name")
                }
                entry.label <- c(entry.label, tmp)
            }
        }
        if (length(out.label) != neuron.count[length(neuron.count)]) {
            if (length(out.label) < neuron.count[length(neuron.count)]) {
                tmp <- NULL
                for (i in 1:(neuron.count[length(neuron.count)] - 
                  length(out.label))) {
                  tmp <- c(tmp, "no name")
                }
                out.label <- c(out.label, tmp)
            }
        }
        grid.newpage()
        for (k in 1:length(weights)) {
            for (i in 1:neuron.count[k]) {
                y.position[k] <- y.position[k] + y.step[k]
                y.tmp <- 0
                for (j in 1:neuron.count[k + 1]) {
                  y.tmp <- y.tmp + y.step[k + 1]
                  result <- calculate.delta(c(x.position[k], 
                    x.position[k + 1]), c(y.position[k], y.tmp), 
                    radius)
                  x <- c(x.position[k], x.position[k + 1] - result[1])
                  y <- c(y.position[k], y.tmp + result[2])
                  grid.lines(x = x, y = y, arrow = arrow(length = unit(0.15, 
                    "cm"), type = "closed"), gp = gpar(fill = col.hidden.synapse, 
                    col = col.hidden.synapse, ...))
                  if (show.weights) 
                    draw.text(label = weights[[k]][neuron.count[k] - 
                      i + 2, neuron.count[k + 1] - j + 1], x = c(x.position[k], 
                      x.position[k + 1]), y = c(y.position[k], 
                      y.tmp), xy.null = 1.25 * result, color = col.hidden.synapse, 
                      fontsize = fontsize - 2, ...)
                }
                if (k == 1) {
                  grid.lines(x = c((x.position[1] - arrow.length), 
                    x.position[1] - radius), y = y.position[k], 
                    arrow = arrow(length = unit(0.15, "cm"), 
                      type = "closed"), gp = gpar(fill = col.entry.synapse, 
                      col = col.entry.synapse, ...))
                  draw.text(label = entry.label[(neuron.count[1] + 
                    1) - i], x = c((x.position - arrow.length), 
                    x.position[1] - radius), y = c(y.position[k], 
                    y.position[k]), xy.null = c(0, 0), color = col.entry.synapse, 
                    fontsize = fontsize, ...)
                  grid.circle(x = x.position[k], y = y.position[k], 
                    r = radius, gp = gpar(fill = "white", col = col.entry, 
                      ...))
                }
                else {
                  grid.circle(x = x.position[k], y = y.position[k], 
                    r = radius, gp = gpar(fill = "white", col = col.hidden, 
                      ...))
                }
            }
        }
        out <- length(neuron.count)
        for (i in 1:neuron.count[out]) {
            y.position[out] <- y.position[out] + y.step[out]
            grid.lines(x = c(x.position[out] + radius, x.position[out] + 
                arrow.length), y = y.position[out], arrow = arrow(length = unit(0.15, 
                "cm"), type = "closed"), gp = gpar(fill = col.out.synapse, 
                col = col.out.synapse, ...))
            draw.text(label = out.label[(neuron.count[out] + 
                1) - i], x = c((x.position[out] + radius), x.position[out] + 
                arrow.length), y = c(y.position[out], y.position[out]), 
                xy.null = c(0, 0), color = col.out.synapse, fontsize = fontsize, 
                ...)
            grid.circle(x = x.position[out], y = y.position[out], 
                r = radius, gp = gpar(fill = "white", col = col.out, 
                  ...))
        }
        if (intercept) {
            for (k in 1:length(weights)) {
                y.tmp <- 0
                x.intercept <- (x.position[k + 1] - x.position[k]) * 
                  intercept.factor + x.position[k]
                for (i in 1:neuron.count[k + 1]) {
                  y.tmp <- y.tmp + y.step[k + 1]
                  result <- calculate.delta(c(x.intercept, x.position[k + 
                    1]), c(y.intercept, y.tmp), radius)
                  x <- c(x.intercept, x.position[k + 1] - result[1])
                  y <- c(y.intercept, y.tmp + result[2])
                  grid.lines(x = x, y = y, arrow = arrow(length = unit(0.15, 
                    "cm"), type = "closed"), gp = gpar(fill = col.intercept, 
                    col = col.intercept, ...))
                  xy.null <- cbind(x.position[k + 1] - x.intercept - 
                    2 * result[1], -(y.tmp - y.intercept + 2 * 
                    result[2]))
                  if (show.weights) 
                    draw.text(label = weights[[k]][1, neuron.count[k + 
                      1] - i + 1], x = c(x.intercept, x.position[k + 
                      1]), y = c(y.intercept, y.tmp), xy.null = xy.null, 
                      color = col.intercept, alignment = c("right", 
                        "bottom"), fontsize = fontsize - 2, ...)
                }
                grid.circle(x = x.intercept, y = y.intercept, 
                  r = radius, gp = gpar(fill = "white", col = col.intercept, 
                    ...))
                grid.text(1, x = x.intercept, y = y.intercept, 
                  gp = gpar(col = col.intercept, ...))
            }
        }
        if (information) 
            grid.text(paste("Error: ", round(result.matrix[rep, 
                "error"], 6), "   Steps: ", result.matrix[rep, 
                "steps"], sep = ""), x = 0.5, y = information.pos, 
                just = "bottom", gp = gpar(fontsize = fontsize + 
                  2, ...))
        if (!is.null(file)) {
            weight.plot <- recordPlot()
            save(weight.plot, file = file)
        }
    }
}
calculate.delta <-
function (x, y, r) 
{
    delta.x <- x[2] - x[1]
    delta.y <- y[2] - y[1]
    x.null <- r/sqrt(delta.x^2 + delta.y^2) * delta.x
    if (y[1] < y[2]) 
        y.null <- -sqrt(r^2 - x.null^2)
    else if (y[1] > y[2]) 
        y.null <- sqrt(r^2 - x.null^2)
    else y.null <- 0
    c(x.null, y.null)
}
draw.text <-
function (label, x, y, xy.null = c(0, 0), color, alignment = c("left", 
    "bottom"), ...) 
{
    x.label <- x[1] + xy.null[1]
    y.label <- y[1] - xy.null[2]
    x.delta <- x[2] - x[1]
    y.delta <- y[2] - y[1]
    angle = atan(y.delta/x.delta) * (180/pi)
    if (angle < 0) 
        angle <- angle + 0
    else if (angle > 0) 
        angle <- angle - 0
    if (is.numeric(label)) 
        label <- round(label, 5)
    vp <- viewport(x = x.label, y = y.label, width = 0, height = , 
        angle = angle, name = "vp1", just = alignment)
    grid.text(label, x = 0, y = unit(0.75, "mm"), just = alignment, 
        gp = gpar(col = color, ...), vp = vp)
}
