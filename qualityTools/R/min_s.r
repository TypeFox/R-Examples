fracChoose = function() {
    DB = FALSE
    genList = list(6 * 9)
    genList = list(c("C = AB"), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c("D = ABC"), c("D = AB", "E = AC"), c("D = AB", 
        "E = AC", "F = BC"), c("D = AB", "E = AC", "F = BC", "G = ABC"), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c("E = ABCD"), c("E = ABC", "F = BCD"), 
        c("E = ABC", "F = BCD", "G = ACD"), c("E = BCD", "F = ACD", "G = ABC", "H = ABD"), c("E = ABC", "F = BCD", "G = ACD", "H = ABD", "J = ABCD"), c("E = ABC", 
            "F = BCD", "G = ACD", "H = ABD", "J = ABCD", "K = AB"), c("E = ABC", "F = BCD", "G = ACD", "H = ABD", "J = ABCD", "K = AB", "L = AC"), c(NULL), 
        c(NULL), c(NULL), c("F = ABCDE"), c("F = ABCD", "G = ABDE"), c("F = ABC", "G = ABD", "H = BCDE"), c("F = BCDE", "G = ACDE", "H = ABDE", "J = ABCE"), 
        c("F = ABCD", "G = ABCE", "H = ABDE", "J = ACDE", "K = BCDE"), c("F = ABC", "G = BCD", "H = CDE", "J = ACD", "K = AEF", "L = ADEF"), c(NULL), c(NULL), 
        c(NULL), c(NULL), c("G = ABCDEF"), c("G = ABCD", "H = ABEF"), c("G = ABCD", "H = ACEF", "J = CDEF"), c("G = BCDF", "H = ACDF", "J = ABDE", "K = ABCE"), 
        c("G = CDE", "H = ABCD", "J = ABF", "K = BDEF", "L = ADEF"), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c("H = ABCDEFG"), c("H = ACDFG", "J = BCEFG"), 
        c("H = ABCG", "J = BCDE", "K = ACDF"), c("H = ABCG", "J = BCDE", "K = ACDF", "L = ABCDEFG"))
    resList = list(6 * 9)
    resList = list(c(3), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(4), c(3), c(3), c(3), c(NULL), c(NULL), c(NULL), 
        c(NULL), c(NULL), c(NULL), c(5), c(4), c(4), c(4), c(3), c(3), c(3), c(NULL), c(NULL), c(NULL), c(6), c(4), c(4), c(4), c(4), c(4), c(NULL), c(NULL), 
        c(NULL), c(NULL), c(7), c(5), c(4), c(4), c(4), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(8), c(6), c(5), c(5))
    facMat = matrix(rep(3:11, 6), ncol = 9, byrow = TRUE)
    runMat = matrix(c(rep(2^2, 9), rep(2^3, 9), rep(2^4, 9), rep(2^5, 9), rep(2^6, 9), rep(2^7, 9)), ncol = 9, byrow = TRUE)
    par(mfrow = c(6, 9))
    par(mar = c(0, 0, 0, 0))
    par(oma = c(4, 4, 4, 4))
    colList = vector(mode = "list")
    colList[3] = "red"
    colList[4] = "yellow"
    colList[5] = "green"
    colList[6] = "green"
    colList[7] = "green"
    colList[8] = "green"
    k = 3
    N = 2^2
    m = 0
    for (i in seq(along = genList)) {
        res = unlist(resList[[i]])
        plot(0, 0, xaxs = "i", yaxs = "i", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "", bg = "red", fg = "green")
        box()
        if (!is.null(res)) 
            rect(0, 0, 1, 1, col = colList[[res]])
        yPos = 0.04
        xPos = 0.6
        item = rev(genList[[i]])
        for (j in seq(along = item)) {
            text(xPos, yPos + (j - 1) * 0.125, item[j], adj = c(0, 0), cex = 0.8)
        }
        p = log2((2^k)/N)
        if (!is.null(res)) {
            romNum = as.character(as.roman(res))
            text(0.1, 0.9, do.call("expression", list(substitute(2[romNum]^(k - p), list(k = k, p = p, romNum = romNum)))), adj = c(0, 1), cex = 1.5)
        }
        k = k + 1
        if ((i%%9) == 0) {
            N = N * 2
            k = 3
        }
    }
    cat("\nChoose a fractional factorial design by clicking into the appropriate field")
    cat("\nWaiting for your selection:")
    cat("\n\n")
    flush.console()
    mtext("number of runs N", side = 2, line = 2.5, outer = TRUE)
    mtext("number of variables k", side = 3, line = 2.5, outer = TRUE)
    for (numFac in 1:9) {
        mtext(numFac + 2, at = c(-0.05 + (1/9) * numFac), outer = TRUE, line = 0.5)
    }
    for (k in 1:6) {
        mtext(2^(k + 1), at = (7 - k)/6 - (0.5 * (1/6)), side = 2, outer = TRUE, line = 0.5)
    }
    if (DB) 
        cat("TODO: standardize the locator part in respect to possible figure region")
    xyList = NULL
    xyList = try(locator(1), silent = TRUE)
    x = 1
    y = 1
    if (!is.null(xyList)) {
        x = ceiling(xyList$x + 8)
        y = ceiling(6 - xyList$y)
    }
    mat = matrix(1:54, ncol = 9, byrow = TRUE)
    fdo = NULL
    if (!(x %in% 1:ncol(mat)) || !(y %in% 1:nrow(mat))) 
        return(fracDesign(k = 3, gen = NULL, replicates = 1))
    else index = mat[y, x]
    k = facMat[y, x]
    generator = genList[[index]]
    N = runMat[y, x]
    if (!is.null(generator)) {
        fdo = try(do.call("fracDesign", list(k = k, gen = generator)), silent = TRUE)
    }
    if (N >= 2^k & is.null(generator)) {
        replicates = N/(2^k)
        fdo = try(fracDesign(k = k, gen = NULL, replicates = replicates), silent = TRUE)
    }
    if (class(fdo) == "facDesign") 
        return(fdo)
    else return(genList[[mat[y, x]]])
} 
