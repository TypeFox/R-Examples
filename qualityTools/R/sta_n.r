.nblock = function(fdo) {
    if (class(fdo) != "facDesign") 
        stop(paste(deparse(substitute(data)), "needs to be an object of class 'facDesign'"))
    return(length(unique(fdo@block[[1]])))
}
starDesign = function(k, p = 0, alpha = c("both", "rotatable", "orthogonal"), cs, cc, data) {
    DB = FALSE
    fdo = NULL
    csFrame = NULL
    ccFrame = NULL
    starFrame = NULL
    blocks = 1
    alpha = alpha[1]
    if (DB) 
        print(alpha)
    if (missing(cc)) 
        cc = 1
    if (missing(cs)) 
        cs = 1
    if (!missing(k)) {
        nameVec = LETTERS[1:k]
    }
    if (!missing(data)) {
        fdo = data
        k = ncol(cube(fdo))
        if (class(fdo) != "facDesign") {
            stop(paste(deparse(substitute(data)), "needs to be an object of class 'facDesign'"))
        }
        if (nrow(star(fdo)) > 0) 
            stop(paste("star portion of", deparse(substitute(data)), "not empty"))
        k = length(names(fdo))
        nameVec = names(names(fdo))
        cc = nrow(centerCube(fdo))
        p = ncol(cube(fdo)) - log(nrow(unique(cube(fdo))), 2)
        blocks = .nblock(fdo) + 1
    }
    if (is.numeric(alpha)) 
        a = alpha
    if (alpha == "rotatable") 
        a = .alphaRot(k, p)
    if (alpha == "orthogonal") 
        a = .alphaOrth(k, p, cc = cc, cs = cs)
    if (alpha == "both") {
        found = FALSE
        for (i in seq(along = .rsmOrth)) {
            if (DB) {
                print(k)
                print(blocks)
                print(p)
            }
            if (.rsmOrth[[i]]$k == k) 
                if (.rsmOrth[[i]]$blocks == blocks) 
                  if (.rsmOrth[[i]]$p == p) {
                    found = TRUE
                    cc = .rsmOrth[[i]]$cc
                    cs = .rsmOrth[[i]]$cs
                    p = .rsmOrth[[i]]$p
                    a = .alphaOrth(k, p, cc, cs)
                    break
                  }
        }
        if (!found) {
            return("no starDesign with approximate rotatability and orthogonality available")
        }
    }
    starFrame = .starFrame(k, alpha = a)
    names(starFrame) = nameVec
    if (DB) 
        print(starFrame)
    if (!missing(data)) 
        star(fdo) = starFrame
    if (DB) 
        print("starFrame added")
    if (cs > 0) {
        csFrame = as.data.frame(matrix(0, nrow = cs, ncol = k))
        names(csFrame) = nameVec
        if (!missing(data)) {
            centerStar(fdo) = csFrame
            if (DB) 
                print("csFrame added")
        }
    }
    if (cc > 0) {
        ccFrame = as.data.frame(matrix(0, nrow = cc, ncol = k))
        names(ccFrame) = nameVec
        if (!missing(data)) {
            centerCube(fdo) = ccFrame
            if (DB) 
                print("ccFrame added")
        }
    }
    if (!missing(data)) 
        return(fdo)
    else return(list(star = starFrame, centerStar = csFrame, centerCube = ccFrame))
} 
