smartBoxPlot <-
  function (the.split, groupNames, groupOrder = groupNames, groupColors,
            ylim., pointcol = "blue", col = "grey", round = TRUE, squeeze = 1,
            na.rm = TRUE,
            boxwex = 0.5,
            cex.lab = 1.5, cex.axis = 1.5, staplewex = 1,
            ...)
  {
    remove.missing <-
      function (x)
      {
        if (is.factor(x) || is.character(x))
          return(x[!is.element(x, c("", " ", "NA", NA))])
        else return(x[!is.na(x)])
      }
    boxplotdataJitter <-
      function (y, round = TRUE, nCirclesPerInch = 11.48936, print = FALSE)
      {
        nCirclesPerInch = 11.48936
        linesPerYAxis = nCirclesPerInch * par("pin")[2]
        if (print)
          catn("linesPerYAxis ", linesPerYAxis)
        if (length(y) == 0) {
          return(cbind(x = numeric(0), y = numeric(0)))
        }
        y = sort(y)
        if (round & (max(y) > min(y)))
          ygroup = round(linesPerYAxis * (y - min(y))/(max(y) -
                                                         min(y)))
        else ygroup = y
        ty = table(ygroup)
        maxEqual = max(ty)
        if (print)
          catn("maxEqual ", maxEqual)
        x = unlist(lapply(ty, function(n) cumsum((0:(n - 1)) * (-1)^(1:n))))
        xUserPerInch = (par("usr")[2] - par("usr")[1])/par("pin")[1]
        circleWidthUser = xUserPerInch/nCirclesPerInch
        if (print)
          catn("circleWidthUser ", circleWidthUser)
        xoffset = min(circleWidthUser, 2/3/maxEqual)
        if (print)
          catn("xoffset ", xoffset)
        x = x * xoffset
        return(data.frame(x = x, y = y))
      }

    if (na.rm)
      the.split = lapply(the.split, remove.missing)
    if (missing(groupNames))
      groupNames = names(the.split)
    if (missing(ylim.))
      ylim. = range(unlist(the.split))
    boxplotdata = lapply(the.split, boxplotdataJitter, round = round)
    for (element in 1:length(boxplotdata)) rownames(boxplotdata[[element]]) = NULL
    names(boxplotdata) = groupNames
    maxcircles = max(sapply(boxplotdata, function(group) max(rle(group[,
                                                                       2])$lengths)))
    for (group in groupNames) {
      boxplotdata[[group]][, "x"] = boxplotdata[[group]][,
                                                         "x"]/squeeze
      boxplotdata[[group]][, "x"] = boxplotdata[[group]][,
                                                         "x"] + (which(group == groupNames))
    }
    boxplotdataMatrix = boxplotdata[[1]]
    for (group in 2:length(boxplotdata)) boxplotdataMatrix = rbind(boxplotdataMatrix,
                                                                   boxplotdata[[group]])
    if (missing(groupColors))
      boxcol = "black"
    else {
      print(groupColors)
      if (length(groupColors) == length(groupNames)) {
        boxcol = groupColors
        pointsPerBox = unlist(sapply(boxplotdata, nrow))
        catn("pointsPerBox = ", pointsPerBox)
        pointcol = rep(groupColors, times = pointsPerBox)
      }
      else if (length(groupColors) == 1) {
        pointcol = groupColors
        boxcol = groupColors
      }
      else warning("groupColors ignored;  length = " %&% length(groupColors) %&%
                     " with " %&% length(groupNames) %&% "groups")
    }
    # print(pointcol)
    plot(boxplotdataMatrix, col = pointcol, axes = F, xlab = "",
         ylab = "", xlim = c(1/2, length(boxplotdata) + 1/2),
         ylim = ylim.)
    boxplot(the.split, boxwex = boxwex,
            xlim = c(1 - 0.4, length(the.split) + 0.4),
            ylim = ylim., border = boxcol, names = groupNames,
            cex.lab = cex.lab, cex.axis = cex.axis, staplewex = staplewex, add = TRUE,
            boxcol=boxcol, ...)
    if(drawDensities) {
      densitySplit = lapply(the.split, density)
      for(group in seq(along=densitySplit)) {
        lines(densitySplit[[group]]$y, densitySplit[[group]]$x)
      }
    }
    return(invisible(boxplotdata))
  }
