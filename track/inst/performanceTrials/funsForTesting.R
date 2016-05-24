# Create lots of files with a variety of names, store about 10 different objects, save, stop & start & check

runSaveLoadTest <- function(what=c("create", "verify", "both"), dir="test1", nObjs=11, scale=7, repetitions=1, seed=1, simple=TRUE, clobber=FALSE, options=list(), verbose=1, sep.times=TRUE, obj.types=createTestObj(return.types=TRUE)) {
    if (env.is.tracked())
        stop("global env is currently tracked")
    if (is.numeric(obj.types))
        obj.types <- createTestObj(return.types=TRUE)[obj.types]
    if (any(is.na(obj.types)))
        stop("have NA obj.types")
    what <- match.arg(what)
    if ((what=="create" || what=="both") && file.exists(dir))
        if (clobber)
            unlink(dir, recursive=TRUE)
        else
            stop(dir, " already exists")
    track.start(dir, options=options)
    on.exit(track.stop())
    if (what=="verify") {
        test1opts <- get("test1opts", envir=globalenv(), inherits=FALSE)
        if (missing(seed))
            seed <- test1opts$seed
        if (missing(nObjs))
            nObjs <- test1opts$nObjs
        if (missing(scale))
            scale <- test1opts$scale
        if (missing(simple))
            simple <- test1opts$simple
        if (missing(obj.types))
            obj.types <- test1opts$obj.types
        obj.names <- test1opts$obj.names
    } else {
        set.seed(seed)
        obj.names <- character(nObjs)
        i <- 0
        while (i < nObjs) {
            n <- randObjName(simple=simple)
            if (!is.element(n, obj.names[seq(len=i)])) {
                i <- i+1
                obj.names[i] <- n
            }
        }
        assign("test1opts", list(seed=seed, nObjs=nObjs, scale=scale, simple=simple, obj.names=obj.names, obj.types=obj.types), envir=globalenv())
        track(list="test1opts")
    }
    start.time <- Sys.time()
    cat("Running '", switch(what, both="create", what), "' test w ", nObjs, " objs in dir '", dir, "', seed=", seed, ", scale=", scale, " & ", if (!simple) "non-", "simple names\n", sep="")
    if (verbose) {
        cat("track.options=\n", paste("  ", deparse(track.options(c("cache", "writeToDisk", "maintainSummary", "alwaysSaveSummary", "recordAccesses")), control=character(0)), "\n", sep=""), sep="")
        cat("Object types= ", paste(obj.types[seq(len=min(length(obj.types), nObjs))], collapse=", "), "\n", sep="")
    }
    total.obj.size <- 0
    total.access.time <- 0
    for (k in seq(len=repetitions)) {
        set.seed(seed)
        rep.start.time <- Sys.time()
        rep.obj.size <- 0
        access.time <- 0
        for (i in seq(len=nObjs)) {
            obj.gen <- createTestObj(obj.types[((i-1) %% length(obj.types))+1], scale=scale)
            if (is.element(what, c("create", "both"))) {
                track(list=obj.names[i])
                if (sep.times)
                    access.time <- access.time + system.time(assign(obj.names[i], value=obj.gen, pos=1))[3]
                else
                    assign(obj.names[i], value=obj.gen, pos=1)
            } else {
                if (sep.times)
                    access.time <- access.time + system.time(obj.old <- get(obj.names[i], inherits=FALSE, pos=1))[3]
                else
                    obj.old <- get(obj.names[i], inherits=FALSE, pos=1)
                same <- identical(obj.gen, obj.old)
                if (verbose>1 || !same)
                    cat("obj ", obj.names[i], " is ", if (same) "same" else "different", "\n", sep="")
            }
            rep.obj.size <- rep.obj.size + object.size(obj.gen)
        }
        rep.end.time <- Sys.time()
        rep.time <- as.numeric(rep.end.time - rep.start.time)
        total.obj.size <- total.obj.size + rep.obj.size
        total.access.time <- total.access.time + access.time
        if (repetitions > 1)
            if (sep.times)
                cat("Rep: Processed ", round(rep.obj.size/2^10), " Kb (", nObjs, " objects) in ",
                    round(access.time, digits=1), "/",
                    round(rep.time-access.time, digits=1),
                    " secs (",
                    round(rep.obj.size/2^10 / access.time), "/",
                    round(rep.obj.size/2^10 / (rep.time-access.time)),
                    " Kb/s)\n", sep="")
            else
                cat("Rep: Processed ", round(rep.obj.size/2^10), " Kb (", nObjs, " objects) in ",
                    round(rep.time, digits=1), " secs (",
                    round(rep.obj.size/2^10 / rep.time), " Kb/s)\n", sep="")
    }
    if (verbose>0)
        print(track.summary())
    on.exit()
    stop.time <- system.time(track.stop())[3]
    end.time <- Sys.time()
    total.time <- as.numeric(end.time - start.time)
    if (sep.times)
        cat("Total: Processed ", round(total.obj.size/2^10), " Kb (", nObjs, " objects) in ",
            round(total.access.time, digits=1), "/",
            round(total.time-total.access.time, digits=1),
            " secs (",
            round(total.obj.size/2^10 / total.access.time), "/", 
            round(total.obj.size/2^10 / (total.time-total.access.time)),
            " Kb/s)\n", sep="")
    else
        cat("Total: Processed ", round(total.obj.size/2^10), " Kb (", nObjs, " objects) in ",
            round(total.time, digits=1), " secs (",
            round(total.obj.size/2^10 / total.time), " Kb/s)\n", sep="")
    cat("Stop time= ", round(stop.time, digits=1), " secs (",
            round(total.obj.size/2^10 / stop.time), " Kb/s)\n", sep="")
    if (what=="both")
        runSaveLoadTest("verify", dir=dir, verbose=max(0, verbose-1))
    invisible(NULL)
}

randObjName <- function(len.exp.mean=8, simple=FALSE, max.len=NA) {
    # Create a random name for object.
    # If simple=TRUE, the name will follow the rules
    # of the trackObjs package for simple names.
    len <- ceiling(rexp(1, 1/len.exp.mean))
    if (!is.na(max.len))
        len <- min(len, max.len)
    if (simple) {
        len <- min(len, 55)
        repeat {
            nm <- paste(c(sample(letters, 1), sample(c(letters, 0:9, ".", "-"), rep=TRUE, size=len-1)), collapse="")
            if (trackObjs:::isSimpleName(nm))
                return(nm)
        }
    } else {
        return(paste(sample(c(letters, LETTERS, 0:9, ".", "_", "-", "&", "("), rep=TRUE, size=len), collapse=""))
    }
}

createTestObj <- function(code, scale, tiny=FALSE, return.types=FALSE) {
    # Create a test object with one of 11 different types (specfied
    # by 'code' as an integer or string).
    # Object size is approx 1kb * 2^scale.
    types <- c("vector", "named.vector", "matrix", "named.matrix", "array", "named.array", "list", "data.frame", "POSIXct", "factor", "named.factor")
    if (return.types)
        return(types)
    if (is.character(code)) {
        i <- pmatch(code, types)
        if (is.na(i))
            stop("code must be one of ", paste(types, collapse=", "))
        code <- types[i]
    }

    if (is.numeric(code))
        if (code < 0 || code > length(types))
            stop("code must be between 0 and ", length(types))
        else
            code <- types[code]
    if (code=="vector" || code=="named.vector") {
        # vector
        obj <- rnorm(max(1,round((if (tiny) 1 else (if (code=="vector") 128 else 24)) * 2^scale)))
        if (code=="named.vector")
            names(obj) <- make.names(rep(letters, length=length(obj)), unique=TRUE)
    } else if (code=="matrix" || code=="named.matrix") {
        # matrix
        d <- pmax(1,round(c(4,2) * (if (tiny) 1 else 4) * 2^(scale/2)))
        obj <- matrix(rnorm(prod(d)), nrow=d[1], ncol=d[2])
        if (code=="named.matrix")
            dimnames(obj) <- list(make.names(rep(letters, length=d[1], unique=TRUE)),
                                  make.names(rep(LETTERS, length=d[2], unique=TRUE)))
    } else if (code=="array" || code=="named.array") {
        # 3-d array
        d <- pmax(1,round(c(3,4,2) * (if (tiny) 1 else 1.75) * 2^(scale/3)))
        obj <- array(rnorm(prod(d)), dim=d)
        if (code=="named.array")
            dimnames(obj) <- list(make.names(rep(letters, length=d[1], unique=TRUE)),
                                  make.names(rep(LETTERS, length=d[2], unique=TRUE)),
                                  make.names(rep(paste("T", letters, sep=""), length=d[3], unique=TRUE)))
    } else if (code=="list") {
        # recursive named list (create as a tree with depth=scale-2, and terminal nodes containing vector, vector, matrix, 3-d array
        f <- function(scale, tiny) {
            if (scale <= 2)
                list(vector1=createTestObj(1,scale-2,tiny=tiny),
                     vector2=createTestObj(2,scale-2,tiny=tiny),
                     vector3=createTestObj(2,scale-2,tiny=tiny),
                     matrix=createTestObj(4,scale-2,tiny=tiny),
                     array=createTestObj(6,scale-2,tiny=tiny))
            else
                list(left=f(scale-1, tiny),
                     right=f(scale-1, tiny))
        }
        obj <- f(scale-if (.Machine$sizeof.pointer<8) 1 else 2, tiny=tiny)
    } else if (code=="data.frame") {
        obj <- data.frame(X=createTestObj(1,scale-2,tiny=tiny),
                          Y=createTestObj(1,scale-2,tiny=tiny),
                          Z=createTestObj(1,scale-2,tiny=tiny),
                          W=createTestObj(9,scale-2,tiny=tiny))
        # data frame
    } else if (code=="POSIXct") {
        # POSIXct
        obj <- strptime("1970-01-01", format="%Y-%m-%d") + seq(len=(if (tiny) 1 else 128) * 2^scale)
    } else if (code=="factor" || code=="named.factor") {
        # factor
        obj <- factor(LETTERS)[sample(26, max(1,round((if (tiny) 1 else (if (code=="factor") 256 else 26)) * 2^scale)), rep=TRUE)]
        if (code=="named.factor")
            names(obj) <- make.names(rep(letters, length=length(obj)), unique=TRUE)
    } else {
        stop("unknown code", code)
    }
    obj
}
