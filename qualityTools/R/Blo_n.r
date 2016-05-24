.letToNum = function(charVec) {
    charVec = toupper(charVec)
    numVec = numeric(length = length(charVec))
    if (!is.character(charVec)) 
        stop(".letToNum: characterVector needs to be provided!")
    alphabet = LETTERS[1:26]
    for (i in seq(along = charVec)) {
        numVec[i] = match(charVec[i], alphabet)
    }
    return(numVec)
}
.letterIndex = function(char) {
    if (char %in% LETTERS[1:26]) {
        return((1:26)[LETTERS[1:26] == char])
    }
    stop("no valid LETTER specified!")
}
.isEven = function(x) {
    if (x%%2 > 0) 
        return(FALSE)
    return(TRUE)
}
.isOdd = function(x) {
    return(!.isEven(x))
}
.lociv = function(charVec) {
    lenVec = numeric(length = length(charVec))
    for (i in seq(along = charVec)) {
        lenVec[i] = length(strsplit(charVec[i], split = "")[[1]])
    }
    return(lenVec)
}
.confoundings = function(blockGenVec, lSet, DB = FALSE) {
    biVec = character(0)
    for (i in 2:length(blockGenVec)) {
        mat = combn(blockGenVec, i)
        temp = apply(mat, 2, strsplit, split = "")
        comb = lapply(temp, unlist)
        comb = lapply(comb, c, lSet)
        if (DB) {
            print("here")
            print(comb)
        }
        combFreq = sapply(comb, table)%%2
        combBool = !apply(combFreq, 2, as.logical)
        chars = row.names(combFreq)
        if (DB) 
            print(combBool)
        biTemp = character(0)
        for (j in 1:ncol(combBool)) {
            biTemp = c(biTemp, paste(chars[combBool[, j]], collapse = ""))
        }
        if (DB) 
            print(biTemp)
        biVec = c(biVec, biTemp)
    }
    return(c(blockGenVec, biVec))
}
.rsm = vector(mode = "list", length = 7)
.rsm[[1]] = list(k = 3, blocks = 2, gen = c("ABC"))
.rsm[[2]] = list(k = 3, blocks = 4, gen = c("AB", "AC"))
.rsm[[3]] = list(k = 4, blocks = 2, gen = c("ABCD"))
.rsm[[4]] = list(k = 4, blocks = 4, gen = c("ABC", "ACD"))
.rsm[[5]] = list(k = 4, blocks = 8, gen = c("AB", "BC", "CD"))
.rsm[[6]] = list(k = 5, blocks = 2, gen = c("ABCDE"))
.rsm[[7]] = list(k = 5, blocks = 4, gen = c("ABC", "CDE"))
.rsm[[8]] = list(k = 5, blocks = 8, gen = c("ABE", "BCE", "CDE"))
.rsm[[9]] = list(k = 5, blocks = 16, gen = c("AB", "AC", "CD", "DE"))
.rsm[[10]] = list(k = 6, blocks = 2, gen = c("ABCDEF"))
.rsm[[11]] = list(k = 6, blocks = 4, gen = c("ABCF", "CDEF"))
.rsm[[12]] = list(k = 6, blocks = 8, gen = c("ABEF", "ABCD", "ACE"))
.rsm[[13]] = list(k = 6, blocks = 16, gen = c("ABF", "ACF", "BDF", "DEF"))
.rsm[[14]] = list(k = 6, blocks = 32, gen = c("AB", "BC", "CD", "DE", "EF"))
.rsm[[15]] = list(k = 7, blocks = 2, gen = c("ABCDEFG"))
.rsm[[16]] = list(k = 7, blocks = 4, gen = c("ABCFG", "CDEFG"))
.rsm[[17]] = list(k = 7, blocks = 8, gen = c("ABC", "DEF", "AFG"))
.rsm[[18]] = list(k = 7, blocks = 16, gen = c("ABD", "EFG", "CDE", "ADG"))
.rsm[[19]] = list(k = 7, blocks = 32, gen = c("ABG", "BCG", "CDG", "DEG", "EFG"))
.rsm[[20]] = list(k = 7, blocks = 64, gen = c("AB", "BC", "CD", "DE", "EF", "FG"))
.blockInteractions = function(fdo, blocks = 2, useTable = "rsm") {
    DB = FALSE
    if (!(blocks %in% c(0, 1, 2, 4, 8, 16, 32, 64))) 
        stop("blocks needs to be a power of 2 up to 64!")
    gen = NULL
    if (blocks %in% c(0, 1)) {
        if (DB) 
            print("TODO: Return the Identity as generator")
        return(gen)
    }
    if (length(useTable) > 0) {
        if (!(nrow(unique(cube(fdo))) >= 2^.numFac(fdo))) 
            stop("no blocking of a fractional factorial Design --> block on replicates instead!")
        if (identical(useTable, "rsm")) {
            for (i in seq(along = .rsm)) {
                if (.rsm[[i]]$k == .numFac(fdo) & .rsm[[i]]$blocks == blocks) 
                  return(.rsm[[i]]$gen)
            }
        }
        return(gen)
    }
    bgaci = matrix(nrow = 0, ncol = blocks - 1)
    if (!is.numeric(blocks)) 
        stop("blocks must be an integer")
    numCol = log2(blocks)
    blockGen = character(3)
    lSet = names(names(fdo))
    sSet = vector(mode = "list")
    for (i in length(lSet):2) {
        sSet = c(sSet, combn(lSet, i, simplify = FALSE))
    }
    if (blocks == 2) {
        index = order(sapply(sSet, length), decreasing = TRUE)[1]
        sSet = sapply(sSet, paste, collapse = "")
        return(sSet[index])
    }
    sSet = sapply(sSet, paste, collapse = "")
    if (DB) 
        print(sSet)
    possGen = combn(sSet, numCol, simplify = FALSE)
    for (i in seq(along = possGen)) {
        blockGenVec = unlist(possGen[[i]])
        if (DB) 
            print(blockGenVec)
        if (DB) 
            print(.confoundings(blockGenVec, lSet))
        newRow = .confoundings(blockGenVec, lSet)
        if (!any(newRow %in% c(lSet, ""))) 
            bgaci = rbind(bgaci, .confoundings(blockGenVec, lSet))
    }
    mat = unique(t(apply(bgaci, 1, sort)))
    temp = t(apply(mat, 1, .lociv))
    temp = t(apply(temp, 1, sort))
    ref = temp[1, ]
    index = 1
    for (i in 1:nrow(temp)) {
        if (any((ref - temp[i, ]) < 0)) {
            ref = temp[i, ]
            index = i
        }
    }
    for (i in 1:nrow(temp)) {
        if (!(any(ref - temp[i, ] > 0) | any(ref - temp[i, ] < 0))) {
            index = c(index, i)
        }
    }
    temp = unique((mat[index, ]))
    cat("\nSuggested Effects for Blocking:")
    cat("\n")
    cat(temp[1, 1:numCol])
    cat("\n")
    cat("\nInteractions Confounded with blocks:")
    cat("\n")
    cat(unique(temp[1, ]))
    cat("\n")
    cat("\n Alternate Effects for Blocking:")
    cat(temp[c(-1), 1:numCol])
    cat("\n")
    gen = temp[1, 1:numCol]
    return(gen)
}

#Reihenfolge muss noch auf Standard gesetzt werden
.blockGenCol = function(gen, fdo) {
    DB = FALSE
    blockVec = NULL
    .blockCol = NULL
    genList = gen
    genList = strsplit(genList, split = "")
    .fdo = fdo
    for (i in seq(along = genList)) {
        gen = genList[[i]]
        for (j in seq(along = gen)) {
            genTemp = .fdo[, gen[j]]
            if (j == 1) 
                blockVec = rep(1, length = length(genTemp))
            blockVec = blockVec * genTemp
            if (DB) 
                print(blockVec)
        }
        if (i == 1) 
            .blockCol = data.frame(B1 = blockVec)
        else .blockCol = cbind(.blockCol, blockVec)
    }
    names(.blockCol) = paste("B", 1:ncol(.blockCol), sep = "")
    return(.blockCol)
}
.blockCol = function(.blockGenCol) {
    DB = FALSE
    .blockCol = numeric(nrow(.blockGenCol))
    uniCol = unique(.blockGenCol)
    for (i in 1:nrow(uniCol)) {
        if (ncol(uniCol) == 1) 
            .blockCol[apply(t(as.data.frame(apply(.blockGenCol, 1, "==", uniCol[i, ]))), 2, all)] = i
        else .blockCol[apply(apply(.blockGenCol, 1, "==", uniCol[i, ]), 2, all)] = i
    }
    return(data.frame(Block = .blockCol))
}
randomize = function(fdo, random.seed, so = FALSE) {
    if (missing(random.seed)) 
        set.seed(93275938)
    else set.seed(random.seed)
    j = 1
    temp = runOrd(fdo)
    for (i in sort(unique(block(fdo)[, 1]))) {
        pos = !is.na(match(block(fdo)[, 1], i))
        count = sum(as.numeric(pos))
        if (so) {
            temp[pos, 1] = j:(j + (count - 1))
        }
        else {
            temp[pos, 1] = sample(j:(j + (count - 1)), count)
        }
        j = j + count
    }
    runOrd(fdo) = temp
    return(fdo)
}
blocking = function(fdo, blocks, BoR = FALSE, random.seed, useTable = "rsm", gen) {
    override = FALSE
    Block = data.frame(Block = rep(1, nrow(fdo))) #do not change
    block(fdo) = Block  #do not change
    fdo = randomize(fdo, so = TRUE)
    if (missing(random.seed)) {
        runif(1)
        random.seed = .Random.seed[sample(1:626, 1)]
    }
    if (missing(gen)) 
        gen = NULL
    if (blocks <= 1) {
        Block = data.frame(Block = rep(1, nrow(fdo)))
        block(fdo) = Block
        fdo = randomize(fdo, random.seed = random.seed)
        return(fdo)
    }
    if (nrow(star(fdo)) > 0 | nrow(centerStar(fdo)) > 0) {
        if (blocks == 2) {
            override = TRUE
            fdo = randomize(fdo, so = TRUE)
            numB1 = nrow(cube(fdo)) + nrow(centerCube(fdo))
            numB2 = nrow(fdo) - numB1
            block(fdo) = data.frame(Block = c(rep(1, numB1), rep(2, numB2)))
            #or by using standard order always
    
            blockGen(fdo) = data.frame(B1 = rep(NA, nrow(fdo)))
        }
        if (blocks %in% c(2, 3, 5, 9, 17)) 
            blocks = blocks - 1
        else stop("Blocking not possible")
    }
    else {
        if (!(blocks %in% c(1, 2, 4, 8, 16, 32, 64, 128))) 
            stop("Blocking not possible")
    }
    if (is.null(gen)) 
        gen = .blockInteractions(fdo, blocks, useTable)
    if (is.null(gen) & !override) {
        cat("\n")
        cat(paste("Blocking in", blocks, "blocks not possible!"))
        cat("\n")
        return(fdo)
    }
    if (!override) {
        .blockGenCol = .blockGenCol(gen, fdo)
        .blockCol = .blockCol(.blockGenCol)
        Block = .blockCol#[runOrd(fdo)[,1],] #TODO: fix this in block()
        BlockGenCol = .blockGenCol#[runOrd(fdo)[,1],] #TODO: fix this in block()
        #or by using standard order always
        block(fdo) = Block
        blockGen(fdo) = BlockGenCol
#       block(fdo) = .blockCol
#        blockGen(fdo) = .blockGenCol
    }
    numCC = nrow(centerCube(fdo))
    if (numCC > 0) {
        ccFrame = as.data.frame(matrix(0, nrow = numCC, ncol = ncol(cube(fdo))))
        names(ccFrame) = names(names(fdo))
        centerCube(fdo) = ccFrame
    }
    fdo = randomize(fdo, random.seed = random.seed)
    return(fdo)
} 
