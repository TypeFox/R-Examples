.alphaOrth = function(k, p = 0, cc, cs) {
    alpha = sqrt((2^(k - p) * (2 * k + cs))/(2 * (2^(k - p) + cc)))
    return(alpha)
}
.alphaRot = function(k, p = 0) {
    alpha = (2^(k - p))^0.25
    return(alpha)
}
.centerPoints = function(nrow, k) {
    if (!is.numeric(nrow))
        stop("nrow must be numeric")
    if (!is.numeric(k))
        stop("ncol must be numeric")
    mat = as.data.frame(matrix(rep(0, nrow * k), nrow, ncol = k))
    return(mat)
}
.starFrame = function(k, alpha) {
    if (!is.numeric(k))
        stop("k must be numeric")
    if (!is.numeric(alpha))
        stop("alpha must be numeric")
    .starFrame = as.data.frame(matrix(0, nrow = k * 2, ncol = k))
    for (j in 1:k) {
        for (i in (2 * (j - 1) + 1):(2 * (j - 1) + 2)) {
            .starFrame[i, j] = ((-1)^i) * alpha
        }
    }
    return(.starFrame)
}
rsmDesign = function(k = 3, p = 0, alpha = "rotatable", blocks = 1, cc = 1, cs = 1, fp = 1,
    sp = 1, faceCentered = FALSE) {
    DB = FALSE
    if (blocks > 2^(k - 1) + 1)
        stop("Blocking not possible")
    if (alpha == "rotatable")
        alpha = .alphaRot(k, p)
    if (alpha == "orthogonal")
        alpha = .alphaOrth(k, p, cc = cc, cs = cs)
    if (alpha == "both") {
        found = FALSE
        for (i in seq(along = .rsmOrth)) {
            if (DB) {
                print(k)
                print(blocks)
            }
            if (.rsmOrth[[i]]$k == k)
                if (.rsmOrth[[i]]$blocks == blocks)
                  if (.rsmOrth[[i]]$p == p) {
                    cc = .rsmOrth[[i]]$cc
                    cs = .rsmOrth[[i]]$cs
                    p = .rsmOrth[[i]]$p
                    alpha = .alphaOrth(k, p, cc, cs)
                    found = TRUE
                    break
                  }
        }
        if (!found) {
            return("no design available")
        }
    }
    if (DB) {
        print("Values")
        print(k)
        print(alpha)
        print(cc)
        print(cs)
        print(blocks)
    }
    fdo = facDesign(k = k, p = p, replicates = fp)                              ###      
    if (cc > 0) {
        temp = as.data.frame(matrix(0, nrow = cc, ncol = ncol(cube(fdo))))
        names(temp) = names(cube(fdo))
        centerCube(fdo) = temp
        if (DB)
            print("centerCube")
    }
    if (DB)
        print("star not added")
    temp = .starFrame(k, alpha)
    starportion = data.frame()
    for (i in 1:sp) {
        starportion = rbind(temp, starportion)
    }
    names(starportion) = names(cube(fdo))
    star(fdo) = starportion
    if (DB)
        print("star added")
    if (cs > 0) {
        temp = as.data.frame(matrix(0, nrow = cs, ncol = ncol(cube(fdo))))
        names(temp) = names(cube(fdo))
        centerStar(fdo) = temp
    }
#    return(fdo)
    fdo = blocking(fdo, blocks)
    return(fdo)
}
.rsmOrth = vector(mode = "list", length = 7)
.rsmOrth[[1]] = list(k = 2, p = 0, col = 1, row = 2, blocks = 2, cc = 3, cs = 3)
.rsmOrth[[2]] = list(k = 2, p = 0, col = 1, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[3]] = list(k = 3, p = 0, col = 2, row = 3, blocks = 3, cc = 2, cs = 2)
.rsmOrth[[4]] = list(k = 3, p = 0, col = 2, row = 2, blocks = 2, cc = 2, cs = 2)
.rsmOrth[[5]] = list(k = 3, p = 0, col = 2, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[6]] = list(k = 4, p = 0, col = 3, row = 4, blocks = 5, cc = 2, cs = 2)
.rsmOrth[[7]] = list(k = 4, p = 0, col = 3, row = 3, blocks = 3, cc = 2, cs = 2)
.rsmOrth[[8]] = list(k = 4, p = 0, col = 3, row = 2, blocks = 2, cc = 2, cs = 2)
.rsmOrth[[9]] = list(k = 4, p = 0, col = 3, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[10]] = list(k = 5, p = 0, col = 4, row = 5, blocks = 9, cc = 2, cs = 4)
.rsmOrth[[11]] = list(k = 5, p = 0, col = 4, row = 4, blocks = 5, cc = 2, cs = 4)
.rsmOrth[[12]] = list(k = 5, p = 0, col = 4, row = 3, blocks = 3, cc = 2, cs = 4)
.rsmOrth[[13]] = list(k = 5, p = 0, col = 4, row = 2, blocks = 2, cc = 2, cs = 4)
.rsmOrth[[14]] = list(k = 5, p = 0, col = 4, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[15]] = list(k = 5, p = 1, col = 5, row = 4, blocks = 5, cc = 6, cs = 1)
.rsmOrth[[16]] = list(k = 5, p = 1, col = 5, row = 3, blocks = 3, cc = 6, cs = 1)
.rsmOrth[[17]] = list(k = 5, p = 1, col = 5, row = 2, blocks = 2, cc = 6, cs = 1)
.rsmOrth[[18]] = list(k = 5, p = 1, col = 5, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[19]] = list(k = 6, p = 0, col = 6, row = 6, blocks = 17, cc = 1, cs = 6)
.rsmOrth[[20]] = list(k = 6, p = 0, col = 6, row = 5, blocks = 9, cc = 1, cs = 6)
.rsmOrth[[21]] = list(k = 6, p = 0, col = 6, row = 4, blocks = 5, cc = 1, cs = 6)
.rsmOrth[[22]] = list(k = 6, p = 0, col = 6, row = 3, blocks = 3, cc = 1, cs = 6)
.rsmOrth[[23]] = list(k = 6, p = 0, col = 6, row = 2, blocks = 2, cc = 1, cs = 6)
.rsmOrth[[24]] = list(k = 6, p = 0, col = 6, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[25]] = list(k = 6, p = 1, col = 7, row = 5, blocks = 9, cc = 4, cs = 2)
.rsmOrth[[26]] = list(k = 6, p = 1, col = 7, row = 4, blocks = 5, cc = 4, cs = 2)
.rsmOrth[[27]] = list(k = 6, p = 1, col = 7, row = 3, blocks = 3, cc = 4, cs = 2)
.rsmOrth[[28]] = list(k = 6, p = 1, col = 7, row = 2, blocks = 2, cc = 4, cs = 2)
.rsmOrth[[29]] = list(k = 6, p = 1, col = 7, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[30]] = list(k = 7, p = 0, col = 8, row = 6, blocks = 17, cc = 1, cs = 11)
#.rsmOrth[[31]] = list(k = 7, p = 0, col = 8, row = 6, blocks = 17, cc = 1, cs = 11) ###
.rsmOrth[[31]] = list(k = 7, p = 0, col = 8, row = 5, blocks = 9, cc = 1, cs = 11)
.rsmOrth[[32]] = list(k = 7, p = 0, col = 8, row = 4, blocks = 5, cc = 1, cs = 11)
.rsmOrth[[33]] = list(k = 7, p = 0, col = 8, row = 3, blocks = 3, cc = 1, cs = 11)
.rsmOrth[[34]] = list(k = 7, p = 0, col = 8, row = 2, blocks = 2, cc = 1, cs = 11)
.rsmOrth[[35]] = list(k = 7, p = 0, col = 8, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[36]] = list(k = 7, p = 1, col = 9, row = 5, blocks = 9, cc = 1, cs = 4)
rsmChoose = function() {
    DB = FALSE
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    colFun = colorRampPalette(c("yellow", "red"), space = "rgb")
    colPalette = colFun(169)
    numRows = 6
    numCol = 9
    blockVals = c(1, 2, 3, 5, 9, 17)
    factorVals = c(2, 3, 4, 5, 5, 6, 6, 7, 7)
    rsmList = .rsmOrth
    plot.new()
    par(mfrow = c(6, 9))
    par(mar = c(0, 0, 0, 0))
    par(oma = c(4, 4, 4, 4))
    for (i in 1:6) for (j in 1:9) {
        par(mfg = c(i, j))
        plot(0, 0, xaxs = "i", yaxs = "i", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, type = "n",
            xlab = "", ylab = "", bg = "red", fg = "green")
        box()
    }
    for (i in seq(along = rsmList)) {
        temp = rsmList[[i]]
        par(mfg = c(temp$row, temp$col))
        par(mfg = c(temp$row, temp$col))
        plot(0, 0, xaxs = "i", yaxs = "i", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, type = "n",
            xlab = "", ylab = "", bg = "red", fg = "green")
        rect(0, 0, 1, 1, col = colPalette[2^((temp$k) - (temp$p))])
        text(0.1, 0.9, paste("N =", 2^((temp$k) - (temp$p)) + temp$cc * (temp$blocks - 1) + temp$cs +
            (temp$k + temp$p) * 2), adj = c(0, 1), cex = 1.25)
        text(0.1, 0.75, paste("k =", temp$k), adj = c(0, 1), cex = 1.25)
        text(0.1, 0.6, paste("p =", temp$p), adj = c(0, 1), cex = 1.25)
        text(0.1, 0.45, ".centerPoints", adj = c(0, 1), cex = 1.25)
        text(0.1, 0.3, paste("Cube:", temp$cc), adj = c(0, 1), cex = 1.25)
        text(0.1, 0.15, paste("Axial:", temp$cs), adj = c(0, 1), cex = 1.25)
        box()
    }
    x = 1/18 + (0:8) * 2/18
    mtext(factorVals, at = x, side = 3, line = 0.5, outer = TRUE)
    mtext("number of factors k", at = 0.5, side = 3, line = 2.5, outer = TRUE)
    x = 1/12 + (5:0) * 2/12
    mtext(blockVals, at = x, side = 2, line = 0.5, outer = TRUE, las = 2)
    mtext("number of blocks", at = 0.5, side = 2, line = 2.5, outer = TRUE)
    cat("\nChoose a response surface design by clicking into the appropriate field")
    cat("\nWaiting for your selection:")
    cat("\n\n")
    flush.console()
    if (DB)
        cat("TODO: standardize the locator part in respect to possible figure region")
    x = numeric(0)
    y = numeric(0)
    xyList = locator(1)                                                         ###
    print(xyList)
    x = ceiling(xyList$x + 8)
    y = ceiling(5 - xyList$y)
    if (DB) {
        print(paste("x:", x))
        print(paste("y:", y))
    }
    if (length(x) < 1)
        return(rsmDesign(k = 2, p = 0, blocks = 2, alpha = "both"))
    if (length(y) < 1)
        return(rsmDesign(k = 2, p = 0, blocks = 2, alpha = "both"))
#    if (!(x %in% factorVals) || !(y %in% blockVals))                           ###
#        return(rsmDesign(k = 2, p = 0, blocks = 2, alpha = "both"))            ###
    blocks = blockVals[y]
    k = factorVals[x]
    if (x==5 || x==7 || x==9 )                                                  ###
     p = 1                                                                      ###
    else                                                                        ###
     p = 0                                                                      ###
    if (DB) {
        print(paste("blocks:", blocks))
        print(paste("k:", k))
  #      print(rsmList)                                                         ###
    }
    for (i in seq(along = rsmList)) {
        if (rsmList[[i]]$k == k)
            if (rsmList[[i]]$blocks == blocks)
                if (rsmList[[i]]$p == p)                                        ###
                   return(rsmDesign(k = k, p = rsmList[[i]]$p, blocks = blocks, ###
                          alpha = "both", cc = rsmList[[i]]$cc,                 ###
                          cs = rsmList[[i]]$cs))                                ###
    }
    return(cat("\nno selection recognized\n"))
}