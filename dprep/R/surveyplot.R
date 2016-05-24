surveyplot <-
function (datos, dataname = "", orderon = 0, class = 0, obs = rep(0, 
    0), lwd = 1) 
{
    if (orderon == 0) 
        datos = datos[order(datos[, ncol(datos)], decreasing = FALSE), 
            ]
    data = datos
    classes = as.numeric(factor(datos[, ncol(datos)]))
    r1 = dim(datos)[1]
    c1 = dim(datos)[2] - 1
    colors = classes + 1
    r = r1
    c = c1
    cnames = colnames(datos)[1:c]
    rnames = rownames(datos)
    graphtitle = paste("Survey Plot for ", dataname)
    if (orderon != 0) {
        neworder = order(data[, orderon], decreasing = TRUE)
        data = data[neworder, ]
        classes = classes[neworder]
        if (class == 0) 
            colors = colors[neworder]
    }
    data <- apply(data[, 1:c], 2, function(data) (data - min(data))/(max(data) - 
        min(data)))
    colnames(data) = cnames
    if (class != 0) {
        data = data[(classes == class), ]
        c = ncol(data)
        r = nrow(data)
        rnames = rownames(data)
        colors = 4
    }
    width1 = seq(0.01, (c) * (0.04), by = 0.04)
    x = rep(1:c, r)
    x = x[order(x, decreasing = FALSE)]
    x = matrix(x, ncol = c, byrow = FALSE)
    for (i in 1:nrow(x)) x[i, 2:ncol(x)] = x[i, 2:ncol(x)] + 
        width1[1:c - 1]
    temp = seq(from = 0, to = 1, length = r + 1)
    y = rep(temp, c)
    y = matrix(y, ncol = c, byrow = FALSE)
    y = y[-1, ]
    op = par(bg = gray(0.8), xaxs = "i", yaxs = "i", yaxp = c(0, 
        1, r + 1), las = 2, cex.axis = 0.75)
    plot(x, y, xlim = c(1, (c + 2)), axes = FALSE, type = "n", 
        xlab = "Features", ylab = "", main = graphtitle)
    axislabels = c(colnames(data), " ")
    width2 = seq(1, c + 1)
    width1 = c(0, width1)
    width2 = width2 + width1
    for (i in width2) lines(c(i, i), c(0, 1), col = "dark gray")
    segments(x, y, (x + data), y, col = colors, lwd = lwd)
    axis(1, at = width2, labels = axislabels, pos = 0)
    if (length(obs) != 0) {
        old.obs = obs
        obs = rep(0, 0)
        if (orderon != 0) {
            if (class == 0) 
                for (i in 1:length(old.obs)) obs = c(obs, which(neworder == 
                  old.obs[i]))
            else for (i in old.obs) obs = c(obs, which(rnames == 
                i))
        }
        else for (i in old.obs) obs = c(obs, which(rnames == 
            i))
        axis(2, at = y[obs, 1], old.obs, tick = FALSE, cex.axis = 0.5)
        segments(x[obs, 1:c], y[obs, 1:c], x[obs, 1:c] + data[obs, 
            1:c], y[obs, 1:c], col = "dark blue", lwd = lwd + 
            1)
    }
    par(op)
}
