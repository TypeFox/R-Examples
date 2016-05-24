seqpcplot <- function(seqdata, group = NULL, weights = NULL,
                      cex = 1, lwd = 1/4, cpal = NULL, grid.scale = 1/5,
                      ltype = "unique", embedding = "most-frequent",
                      lorder = NULL , lcourse = "upwards",
                      filter = NULL, hide.col = "grey80",
                      alphabet = NULL, missing = "auto", order.align = "first",
                      title = NULL, xlab = NULL, ylab = NULL,
                      xaxis = TRUE, yaxis = TRUE, axes = "all",
                      xtlab = NULL, cex.plot = 1,
                      rows = NA, cols = NA, plot = TRUE, seed = NULL,
                      ...) {

  seqpcplot_private(seqdata = seqdata, group = group, weights = weights,
                    cex = cex, lwd = lwd, cpal = cpal,
                    grid.scale = grid.scale,
                    ltype = ltype, embedding = embedding,
                    lorder = lorder, lcourse = lcourse,
                    filter = filter, hide.col = hide.col,
                    alphabet = alphabet, missing = missing,
                    order.align = order.align,
                    title = title, xlab = xlab, ylab = ylab,
                    xaxis = xaxis, yaxis = yaxis,
                    axes = axes, xtlab = xtlab,
                    cex.plot = cex.plot,
                    rows = rows, cols = cols, plot = plot, seed = seed,
                    ...)

}

seqpcplot_private <- function(seqdata, weights = NULL, group,
                              cex = 1, lwd = 1/4, cpal = NULL,
                              grid.scale = 1/5, grid.fill = "grey95",
                              grid.lwd = 0.5, grid.shape = "default",
                              grid.border = "grey60", grid.col = "white",
                              border = NULL, border.lwd = 0,
                              lorder = NULL, lcourse = "upwards",
                              hide.col = "grey80", col.nobs = "black",
                              filter = NULL, alpha = 1,
                              ltype = "unique", embedding = "most-frequent",
                              sf.cex = 1, sf.cex.leaves = 1,
                              title = NULL, xlab = NULL, ylab = NULL,
                              xlim, ylim,
                              alphabet = NULL, alphabet.optim = FALSE,
                              missing = c("auto", "show", "hide"),
                              R = 1000, order.align = NULL, maxit = 300,
                              xtlab = xtlab,
                              xaxis = TRUE, yaxis = TRUE, axes = "all",
                              cex.plot = 1, rows = NA, cols = NA,
                              plot = TRUE, seed = NULL, add = FALSE,
                              verbose = FALSE, ...) {

  if (!"seqpcplot" %in% class(seqdata)) {

    ## Step 1: Check arguments ............................. #

    if (verbose) cat(" [>] check arguments\n")

    if (!(embedding %in% c("most-frequent", "uniformly"))) {
      stop("[!] invalid embedding input")
    }

    if (!(lcourse %in% c("upwards", "downwards"))) {
      stop("[!] invalid lcourse input")
    }

    if (!(order.align %in% c("first", "last", "time"))) {
      stop("[!] invalid lcourse input")
    }

    if (is.null(lorder)) {
      lorder <- ifelse(is.null(filter), "background", "foreground")
    }

    if (!(lorder %in% c("background", "foreground"))) {
      stop("[!] invalid lorder input")
    }
    
    if (!is.null(filter)) {
      filter <- construct.filter(x = filter)
    }

    missing <- match.arg(missing)
    
    mtext <- NULL
    if (is.list(filter) && # text below the title
        filter$type == "function" &&
        is.character(filter$value) &&
        filter$value %in% c("minfreq", "cumfreq")) {
      mtext <- "colored: "
    }

    ## set seed
    if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    if (!is.null(seed)) set.seed(seed)
    RNGstate <- .Random.seed

    ## Step 2: Check and prepare raw data .................. #

    if (!missing(seqdata)) {

      ## seqe format
      if (inherits(seqdata, "seqelist")) { # convert seqe data

        ##TMP <- TraMineR:::seqe2TSE(seqdata)
        TMP <- seqe2TSE(seqdata)
        id <- factor(TMP$id)
        if (nlevels(id) < length(seqdata)) {
          levels(id) <- c(levels(id), seq(nlevels(id), length(seqdata), 1))
        }
        x <- TMP$timestamp
        y <- TMP$event

        ## delete 'end' event
        if ("_end" %in% levels(y)) {
          subs <- y != "_end"
          id <- id[subs]
          x <- x[subs]
          y <- y[subs]
          y <- factor(y, levels = levels(y)[levels(y) != "_end"])
        }
        
        if (is.null(weights)) weights <- seqeweight(seqdata)

        ## STS format (STS or DSS representation)
      } else if (inherits(seqdata, "stslist")) {

        if (order.align %in% c("first", "last")) {
          seqdata <- seqdss(seqdata)
        }
        TMP <- as.data.frame(seqdata)
        
        TMP <- reshape(data = TMP, ids = 1:nrow(TMP), times = colnames(TMP), timevar = "time", varying = list(state = colnames(TMP)), v.names = "state", direction = "long")
        TMP$id <- factor(TMP$id)
        TMP <- TMP[TMP$state %in% c(attr(seqdata, "alphabet"), attr(seqdata, "nr")), ]

        TMP$state <- factor(TMP$state, levels = c(attr(seqdata, "alphabet"), attr(seqdata, "nr")), labels = c(attr(seqdata, "labels"), attr(seqdata, "nr")))

        TMP$time <- factor(TMP$time, levels = attr(seqdata, "names"), labels = attr(seqdata, "names"))
        levels(TMP$time) <- names(seqdata)
        TMP <- TMP[order(TMP$id, TMP$time, TMP$state),]
        id <- TMP$id
        x <- TMP$time
        y <- TMP$state

        if (missing %in% c("auto", "hide")) {
          SUBS <- y == attr(seqdata, "nr")
          if (missing == "hide" | missing == "auto" && !any(SUBS)) {
            id <- id[!SUBS]
            x <- x[!SUBS]
            y <- y[!SUBS]
            y <- factor(y, levels = levels(y)[levels(y) != attr(seqdata, "nr")])
          }
        }
        
        if (is.null(weights)) weights <- attr(seqdata, "weights")

        if (is.null(xlab)) {
          if (substr(attributes(seqdata)$names[1], 1, 2) == "ST") {
            xlab <- "Distinctive States"
          }
        }

        ## TSE format
      } else if (is.data.frame(seqdata)) {

        if (sum(c("id", "time", "event") %in% colnames(seqdata)) == 3) {
          id <- seqdata$id
          y <- seqdata$event
          x <- seqdata$time
        }

      } else {
        stop("[!] invalid seqdata argument")
      }

    } else {

      stop("[!] seqdata argument is missing")

    }

    if (sum(c(length(x), length(y)) != length(id)) != 0) {
      stop("[!] input vectors x, y and id have different lengths")
    }

    ## check the x variable (the time or sequence position)
    if (is.numeric(x)) {
      if (length(unique(x)) > 1) {
        xdiff <- diff(sort(unique(x)))

        if (all((xdiff / min(xdiff) - round(xdiff / min(xdiff),0)) != 0)) {
          x <- factor(x, levels = seq(min(x, na.rm = TRUE),
                           max(x, na.rm = TRUE), min(xdiff)))
        } else if (all(xdiff - round(xdiff) == 0)) {
          x <- factor(x, levels = seq(min(x), max(x), 1))
        } else if (order.align == "time") {
          warning("[!] Problems with distances between x axis positions. The x axis positions will not be illustrated adequately.")
          x <- factor(x)
        }

      } else { x <- factor(x) }  # pseudo case
    }
    if (order.align %in% c("first", "last")) {
      x <- unlist(tapply(X = as.integer(x), INDEX = list(id), FUN = ordering, align = order.align))
      x <- factor(x)
    }

    ## y (state or event categories) must be categorical
    if (!is.factor(y)) stop("[!] y is not categorical")
    if (!is.null(alphabet)) {
      if (sum(levels(y) %in% alphabet) != nlevels(y)) {
        stop("[!] incorrect alphabet")
      } else {
        y <- factor(y, levels = alphabet)
      }
    }

    ## convert id variable to a factor
    if (!is.factor(id)) id <- factor(id)

    ## construct/ convert group variable
    if (!is.null(group)) {
      if (!is.factor(group)) group <- factor(group)
      if (sum(is.na(group)) > 1) {
        group <- factor(group, levels = c(levels(group), "not available"))
        group[is.na(group)] <- "not available"
      }
      if (length(group) == length(id)) {
        group <- unique(data.frame(id, group))[, 2]
      }
      if (length(group) != nlevels(id)) {
        stop("[!] cannot link group and id vector")
      }
      if (nlevels(group) > length(unique(group))) {
        warning(paste("[!] erase empty group levels:", paste(setdiff(levels(group), unique(group)), collapse = ", ")))
        group <- factor(group, levels = levels(group)[levels(group) %in% unique(group)])
      }
    } else group <- factor(rep(1, nlevels(id)))
    names(group) <- levels(id)

    ## construct/ convert weights variable
    if (!is.null(weights)) {
      if (!is.numeric(weights)) stop("[!] weights are not numeric")
      if (min(weights, na.rm = TRUE) < 0) stop("[!] negative weights")
      if (length(weights) != nlevels(id) &&
          length(weights) == length(id)) {
        weights <- unique(data.frame(id, weights))[, 2]
      }
      if (length(weights) != nlevels(id) | (sum(is.na(weights))) > 0) {
        stop(paste("[!] cannot link weight and id vector"))
      }
    } else weights <- rep(1, nlevels(id))

    ## construct/ convert filter variable
    if (!is.null(filter)) {
      if (filter$type == "value") {
        if (!is.numeric(filter$value)) stop("[!] filter$value is not numeric")
        if (min(filter$value, na.rm = TRUE) < 0) stop("[!] negative filter$value values")
        if (length(filter$value) == length(id)) {
          filter$value <- unique(data.frame(id, filter$value))[, 2]
        }
        if (length(filter$value) != nlevels(id) | (sum(is.na(filter$value))) > 0) {
          stop(paste("[!] cannot link filter$value and id vector"))
        }
      }
      if (filter$type == "density") {
        if (!is.matrix(filter$value)) stop("[!] filter$value must be a matrix")
        if (ncol(as.matrix(filter$value)) != nlevels(id)) {
          stop(paste(("[!] cannot link dissimilarity matrix and id vector")))
        }
        if (is.null(filter$kernel)) filter$kernel <- "gaussian"
        if (is.null(filter$bw)) filter$bw <- 0.2
        if (is.character(filter$kernel)) {
          filter$kernel <- switch(filter$kernel, gaussian = dnorm, rectangular = rectangular, dnorm)
        }
        rownames(filter$value) <- colnames(filter$value) <- levels(id)
      }
    }

    ## remove NAs and omit.levels in x, y or id
    SUBS <- is.na(id) | is.na(x) | is.na(y)
    if (sum(SUBS) > 0) {
      warning(paste(" [!] found NA's in x, y or id. see entries ", paste(which(SUBS), collapse = ", ")))
      id <- id[!SUBS]
      x <- x[!SUBS]
      y <- y[!SUBS]
    }

    ## standardize weights
    wid.group.unscaled <- tapply(weights, group, sum)
    TMP1 <- table(group)[as.integer(group)]
    TMP2 <- tapply(weights, group, sum)[as.integer(group)]
    weights <- weights / TMP2 * TMP1
    names(weights) <- levels(id)

    ## id variable
    idlevs <- levels(id)
    idused <- intersect(levels(id), unique(id))
    nid.tot <- nlevels(id)
    id <- factor(id, levels = idused)
    weights.use <- weights[names(weights) %in% idused]
    nid.use <- nlevels(id)
    id <- as.integer(id)

    ## x variable
    if (is.null(xtlab)) xlevs <- levels(x) else xlevs <- rep(xtlab, length.out = nlevels(x))
    nx <- nlevels(x)
    x <- as.integer(x)

    ## y variable
    ylevs <- levels(y)
    ny <- nlevels(y)
    y <- as.integer(y)

    ## group variable
    nid.group.tot <- table(group)
    wid.group.tot <- tapply(weights, group, sum)
    group <- group[names(group) %in% idused]
    nid.group.use <- table(group)
    wid.group.use <- tapply(weights.use, group, sum)
    grouplev <- levels(group)
    ngroup <- nlevels(group)
    group <- as.integer(group)

    ## filter matrix
    if (!is.null(filter)) {
      if (is.matrix(filter$value)) {
        filter$value <- filter$value[rownames(filter$value) %in% idused, colnames(filter$value) %in% idused]
      }
    }

    ## order the data
    SUBS <- order(id, x, y)
    id <- id[SUBS]
    x <- x[SUBS]
    y <- y[SUBS]

    ## Step 3: Find curves to plot ......................... #

    if (verbose) cat(" [>] find plot trajectories\n")

    ## some needed variables
    x.list <- tapply(X=x,
                     INDEX=list(id),
                     FUN=function(x){return(x)})
    traj.string <- tapply(X = y, INDEX = list(id, x), paste,
                          collapse = ",")

    traj.string <- apply(X = traj.string, MARGIN = 1, FUN = createstring)
    traj.order <- order(nchar(traj.string), decreasing = TRUE)
    traj.group <- c()
    traj.id <- c()
    group.id <- c()
    group.name <- c()

    ## start algorithms
    if (ltype == "unique") {
      traj.group <- as.integer(factor(traj.string))
      traj.id <- 1:nid.use
      group.name <- 1:max(traj.group)
      group.id <- tapply(traj.id, traj.group, function(x){x[1]})
    } else if (ltype == "jitter") {
      traj.group <- 1:length(traj.string)
      traj.id <- 1:nid.use
      group.name <- 1:max(traj.group)
      group.id <- 1:nid.use
    } else if (ltype == "non-embeddable") {
      embedded <- c()
      while (sum(!(1:nid.use) %in% embedded) > 0) {
        i <- which(!traj.order %in% embedded)[1]
        gn <- ifelse(length(traj.group) == 0, 1, max(traj.group) + 1)
        group.id <- append(group.id, traj.order[i])
        group.name <- append(group.name, gn)

        ## check for embeddable trajectories
        if (is.null(embedded)) {
          SUBS <- which(unlist(lapply(X = traj.string, FUN = grepl, traj.string[traj.order[i]], fixed = TRUE)))
        } else {
          SUBS <- which(unlist(lapply(X = traj.string[-traj.order[1:(i - 1)]],FUN = grepl, traj.string[traj.order[i]], fixed = TRUE)))
          SUBS <- (1:nid.use)[-traj.order[1:(i - 1)]][SUBS]
        }
        traj.id <- append(traj.id, SUBS)
        traj.group <- append(traj.group, rep(gn, length(SUBS)))
        embedded <- unique(c(embedded, SUBS))
      }
    } else {
      stop("[!] invalid ltype argument")
    }

    ## erase difficult embeddings difficult to interprete in case
    ## of a group vector was entered
    if ((ltype == "non-embeddable" ) & (ngroup > 1)) {
      TMP <- table(traj.id, traj.group)
      uemb <- rowSums(TMP) == 1
      for (i in 1:ngroup) {
        ## for which non-embeddedable trajectories there are observed
        ## trajectories that were solely embedded into it
        uemb.group <-
          traj.group[traj.id %in% which(uemb & (group == i))]
        for (j in which((group == i) & (!uemb)))
          ## erase embeddings
          if (sum(uemb.group %in% traj.group[traj.id %in% j]) > 0) {
            SUBS2 <- (traj.id %in% j) & (!traj.group %in% uemb.group)
            traj.id <- traj.id[!SUBS2]
            traj.group <- traj.group[!SUBS2]
          }
      }
    }

    ## Step 4: Prepare point data ........................ #

    if (verbose) cat(" [>] prepare point data\n")

    ## setup pts object
    pts <- data.frame(id = id[id %in% group.id], x = x[id %in% group.id], y = y[id %in% group.id])
    pts$traj <- as.integer(factor(pts$id, levels = group.id, labels = group.name))
    ntraj <- max(pts$traj)

    ## if lines are demanded to run from the highest to lowest y category
    ## in case of simultaneous observations
    if (lcourse == "downwards") {
      pts <- pts[order(pts$traj, pts$x, factor(pts$y, levels = rev(sort(unique(pts$y))))),]
    }

    ## find multiple equal simultaneous observations
    ## and note their frequency
    TMP <- data.frame(table(x=pts$x, y=pts$y, traj = pts$traj))
    TMP <- TMP[TMP$Freq != 0,]
    pts <- merge(x = unique(pts), y = TMP, by = c("x", "y", "traj"), sort = FALSE)

    ## some useful extractors for later
    frqcols <- paste("n", 1:ngroup, sep = ".")

    ## calculate weighted trajectory shape group sizes
    pts[, frqcols] <- 1
    if (verbose) cat(" [>] determine point sizes\n")
    tab <- table(traj.id,traj.group)
    w <- 1/rowSums(tab)
    minx <- unlist(lapply(X = x.list, FUN = min))
    maxx <- unlist(lapply(X = x.list, FUN = max))

    ## option: trajectories that are embeddable into multiple
    ## non-embeddable trajectories are embedded into the most
    ## frequent non-embeddable trajectory
    if (ltype == "non-embeddable") {
      if (embedding=="most-frequent") {
        TMP <- tapply(X = weights.use[w == 1], INDEX = list(traj.group[traj.id %in% which(w == 1)], factor(group[w == 1],levels = 1:ngroup)), FUN = sum)
        TMP[is.na(TMP)] <- 0
        ## run unique group assignements
        for (i in which(w != 1)) {
          SUBS1 <- traj.group[traj.id == i]
          SUBS1 <- SUBS1[which.max(TMP[SUBS1, group[i]])]
          SUBS2 <- which((traj.group != SUBS1) & (traj.id == i))
          traj.id <- traj.id[-SUBS2]
          traj.group <- traj.group[-SUBS2]
        }
        w <- rep(1,length(w)) # set all
      }
    }

    ## point weights
    for (i in 1:nrow(pts)) {
      if (i==1) {
        SUBS <- traj.id[traj.group == pts$traj[i]]
      } else {
        if (pts$traj[i] != pts$traj[i-1]) {
          SUBS <- traj.id[traj.group == pts$traj[i]]
        }
      }
      pts[i, frqcols] <- unlist(tapply(X = weights.use[SUBS] * w[SUBS] * ((minx[SUBS] <= pts$x[i]) & (maxx[SUBS] >= pts$x[i])), INDEX = list(factor(group[SUBS], levels = 1:ngroup)), FUN = sum))
    }
    pts[is.na(pts)] <- 0


    ## optimisation of y axis
    if (alphabet.optim)
      {
        warning("[!] alphabet optimization is still in development")

        ## if (verbose) cat(" [>] optimize order\n")
        ## permn <- permn(1:ny)
        ## if (length(permn) > R) {
        ##   permn <- permn[sample(1:length(permn), R)]
        ## }

        ## TMP <- rep(0,length(permn))
        ## for (i in 1:length(permn))
        ##   {
        ##     pts$y.tmp <- as.integer(factor(x = pts$y, levels = permn[[i]]))
        ##     for (t in 1:ntraj)
        ##       {
        ##         SUBS <- pts$traj == t
        ##         TMP[i] <- TMP[i] + sqrt(sum(diff(pts$y.tmp[SUBS])^2 + diff(pts$x[SUBS])))
        ##       }
        ##   }
        ## pts$y <- as.integer(factor(x = pts$y, levels = permn[[which.min(TMP)]]))
        ## ylevs <- ylevs[permn[[which.min(TMP)]]]
        ## y <- factor(y, levels = ylevs)

      }


    ## Step 5: Find jitter coordinates ..................... #

    if (verbose) cat(" [>] find translation coordinates\n")

    ## needed variables
    wid.max <- apply(X = pts[,grep("n.",colnames(pts)), drop = FALSE], MARGIN = 2,FUN = tapply, list(pts$traj), max)
    propid.use <- scale(x = matrix(wid.max, ncol = ngroup), center = FALSE, scale = wid.group.use)
    propid.tot <- scale(x = matrix(wid.max, ncol = ngroup), center = FALSE, scale = wid.group.tot)
    propid.tot <- rbind(propid.tot, apply(propid.tot, 2, function(x) 1 - sum(x)))
    propid.use.max <- apply(X = propid.use, MARGIN = 1, FUN = max)
    trajord <- order(propid.use.max, decreasing = TRUE)

    seqpcplot_jitter <- function(data) {

      ## create initial grid
      ngrid <- ngrid0 <- 10
      sl <- ceiling(ngrid * sqrt(propid.use.max)) # size
      sl[sl == 0] <- 1 # minimum (relevant for weights = 0 cases)
      ngrid <- ceiling(sqrt(sum(sl ^ 2)))
      grid <- matrix(0, ngrid, ngrid) # generate initial grid
      gridsubs <- 1:(ngrid ^ 2) # identifiers
      data$xpos <- data$ypos <- NA
      potpos <- c()
      blacklist <- c()

      ## the algorithm
      for (i in 1:length(trajord)) {
        subscripts <- data$traj == trajord[i]
        count <- 0
        found <- FALSE
        while ((!found)&(count <= maxit)) {
          if ((i > 1) & (count == 0)) {
            if (sl[trajord[i]] < sl[trajord[i-1]]) {
              potpos <- c()
              blacklist <- c()
            }
          }
          count2 <- 0
          while (length(potpos) == 0) { # no potential position available
            if (count2 > 0) { # enlarge grid
              ngrid <- ngrid+1
              grid <- rbind(grid, rep(0, ngrid-1))
              grid <- cbind(grid, rep(0, ngrid))
              gridsubs <- seq(1, ngrid^2, 1)
              ## correct blacklist
              blacklist <-  blacklist + floor(blacklist/ngrid)
            }
            ## determine potential grid positions
            potpos <- (grid[gridsubs] == 0) & # occupied positions
            ((gridsubs %% ngrid) < (ngrid - sl[trajord[i]] + 2)) & # to close top
            (gridsubs < (ngrid * (ngrid - sl[trajord[i]] + 1))) # to close left
            if (sl[trajord[i]] > 1) { # filter for group with ns > 1 field
              potpos <- potpos & (gridsubs %% ngrid != 0)
            }
            potpos <- which(potpos)
            potpos <- potpos[!potpos %in% blacklist]
            count2 <- count2+1
          }
          xpotpos <- floor(potpos / ngrid) + 1
          ypotpos <- potpos %% ngrid
          ypotpos[ypotpos == 0] <- ngrid
          posind <- sample(1:length(potpos), 1) # random assignement
          pos <- potpos[posind]
          xpos <- xpotpos[posind]
          ypos <- ypotpos[posind]
          xsubs <- seq(xpos, xpos + sl[trajord[i]] - 1, 1)
          ysubs <- seq(ypos, ypos + sl[trajord[i]] - 1, 1)
          if (sum(c(grid[ysubs, xsubs]) != 0) == 0) { # assign position
            data$xpos[subscripts] <- xpos
            data$ypos[subscripts] <- ypos
            grid[ysubs,xsubs] <- trajord[i] # register assignement in matrix
            potpos <- potpos[!potpos %in% which(grid == trajord[i])]
            found <- TRUE
          } else { # nothing found
            count <- count + 1
            blacklist <- c(blacklist, pos)
            potpos <- potpos[-posind]
          }
          if (count == maxit) {
            TMP1 <- data[data$traj == i,]
            TMP2 <- createstring(y = TMP$y1, x = TMP$x1)
            stop(paste(" [!] found no grid position for trajectory: ", TMP2 ,". Please retry by setting the maxit argument", sep = ""))
          }
        }
      }
      data$xjitter <- (data$xpos - 1) / ngrid
      data$yjitter <- (data$ypos - 1) / ngrid
      ## determine central plot coordinates

      return(list(jitter = data[,c("xjitter", "yjitter")], ngrid0 = ngrid0, ngrid = ngrid))
    }

    TMP <- seqpcplot_jitter(pts)
    pts[,c("xjitter","yjitter")] <- TMP$jitter
    ngrid0 <- TMP$ngrid0
    ngrid <- TMP$ngrid

    ## Step 6: Prepare plot  ............................... #

    wdtcols <- paste("wdt", 1:ngroup, sep = ".")
    colcols <- paste("col", 1:ngroup, sep = ".")
    lwdcols <- paste("lwd", 1:ngroup, sep = ".")

    ## trajectory colors
    if (is.null(cpal)) { # colors are not specified by the user
      cpal <- brewer.pal(8, "Dark2")[rep(seq(1, 8), ceiling(ntraj / 8))][1:ntraj]
    } else {
      cpal <- rep(cpal,length.out = ntraj) # one colour for each trajectory type
    }
    cpal <- cpal[order(trajord)]
    if (is.character(cpal)) {
      cpal <- col2rgb(cpal)
      cpal <- rgb(cpal[1,], cpal[2,], cpal[3,], alpha * 255, maxColorValue = 255)
    }

    pts[, colcols] <- matrix(rep(cpal[pts$traj], ngroup), ncol = ngroup)
    col.nobs <- rep(col.nobs, ngroup)

    ## color gradients
    if (!is.null(filter)) {

      if (verbose) cat(" [>] line colouring\n")

      if (filter$type == "function") {
        TMP1 <- matrix(NA, ncol = ngroup, nrow = ntraj + 1)
        CALL <- list(filter$value)
        CALL <- append(CALL, filter$args)
        mode(CALL) <- "call"
        for (i in 1:ncol(propid.tot)) {
          CALL$x <- propid.tot[,i]
          TMP1[, i] <- eval(CALL)
        }
        TMP3 <- TMP1
        for (i in 1:nrow(TMP1)) {
          for (j in 1:ncol(TMP1)) {
            TMP3[i, j] <- colourize(c(0, TMP1[i, j], max(TMP1[, j])), hide.col, c(cpal, col.nobs[1])[i])[2]
          }
        }
        for (i in 1:ngroup) {
          pts[, colcols[i]] <- TMP3[-nrow(TMP3), i][pts$traj]
        }
        col.nobs <- TMP3[nrow(TMP3), ]
      }

      if (filter$type == "density") { # nobs has not been implemented

        ## colour by distance matrix
        TMP <- rep(NA, nrow(filter$value))
        for (i in 1:ngroup) {
          SUBS <- which(group == i)
          for (i in SUBS) TMP[i] <- 1 / (sum(weights.use[SUBS]) * filter$bw) * sum(weights.use[SUBS] * do.call(what = filter$kernel, args = list(filter$value[i, SUBS] / filter$bw)))
          TMP[SUBS] <- max(TMP[SUBS])-TMP[SUBS]
        }
        filter$value <- TMP # then coloring is done below
      }
      if (is.vector(filter)) {
        if (is.numeric(filter$value)) {

          ## colour by a vector with one numeric value for each individual
          filter$value <- tapply(filter$value[traj.id], list(traj.group, group[traj.id]), function(x) { median(x) })
          if (is.list(filter$value)) stop("[!] filter vector could not be merged")
          if (ngroup == 1) filter$value <- matrix(filter$value, ncol = 1)
          filter$value <- scale(x = filter$value, center = apply(filter$value, 2, min, na.rm = TRUE), scale = apply(filter$value, 2, function(x) { diff(range(x, na.rm = TRUE)) }))
          filter$value[is.na(filter$value)] <- 1
          for (i in 1:ntraj) {
            for (j in 1:ngroup) {
              pts[pts$traj == i, colcols[j]] <- colourize(filter$value[i,j], cpal[i], hide.col)
            }
          }
        } else if (filter$type == "sequence") {

          ## colour whole sequences
          TMP <- as.character(seqecreate(timestamp = pts$x, event = factor(pts$y, 1:ny, ylevs), id = pts$traj))
          TMP <- sapply(X = TMP, FUN = function(x) { substr(x = x, start = unlist(gregexpr("[(]", x[1]))[1], stop = nchar(x)) })
          TMP <- gsub(pattern = "-[[:digit:]]{+}-", replacement = "-",x = TMP)
          SUBS <- unique(pts$traj)[which(TMP %in% filter$value)]
          pts[!pts$traj %in% SUBS, colcols] <- hide.col
        } else if (filter$type == "subsequence") {

          ## colour subsequences
          TMP <- seqecreate(timestamp = pts$x, event = factor(pts$y, 1:ny, ylevs), id = pts$traj)
          subseq <- seqefsub(TMP, strsubseq = filter$value, minSupport = 0)
          msubcount <- seqeapplysub(subseq)
          SUBS <- unique(pts$traj)[apply(msubcount, 1, sum) > 0]
          pts[!pts$traj %in% SUBS, colcols] <- hide.col
        }
      }
    }

    ## title(s)
    if (is.null(title)) {
      if (ngroup > 1) {
        title <- paste(grouplev, ", n = ", round(wid.group.unscaled, 1), sep = "")
      } else title <- paste("n = ", round(wid.group.unscaled, 1), sep = "")
    } else {
      if ((ngroup > 1) & (length(title) < ngroup)) {
        title <- rep(title, length.out = ngroup)
      }
    }

    ## mtext: cumulative frequencies for minfreq and cumfreq filters
    if (!is.null(mtext)) {
      TMP <- apply(propid.tot, 2, function(x, level) { sum(x[filter$value(x, level) == 1]) }, level = filter$args$level)
      mtext <- paste(mtext, round(100 * TMP, 1), "%", sep = "")
    }

    ## xlab and ylab
    if (is.null(xlab))
      xlab <- if (order.align == "time") "Timestamp" else "Position"
    if (is.null(ylab)) ylab <- ""

    ## xlim and ylim parameters
    ## compute pretty xlims and ylims
    if (missing(xlim)) {
      TMP1 <- which.max(propid.tot[nrow(propid.tot),])
      TMP2 <- which.max(pts$n.1)
      TMP <- pts$n.1[TMP2] *
        sqrt(propid.tot[ntraj + 1, TMP1] /
             (pts$n.1[TMP2] / wid.group.tot[1]))
      xlim <- c(1 - sqrt(grid.scale) / 2, # 0.15 * ...  + TMP * cex
                nx + sqrt(grid.scale) / 2)
      ## note: 0.15 is a guess for xf (see plot.seqpcplot)
    }
    if (missing(ylim)) {
      TMP1 <- which.max(propid.tot[nrow(propid.tot),])
      TMP2 <- which.max(pts$n.1)
      TMP <- pts$n.1[TMP2] *
        sqrt(propid.tot[ntraj + 1, TMP1] /
             (pts$n.1[TMP2] / wid.group.tot[1]))
      ylim <- c(1 - sqrt(grid.scale) / 2,
                ny + sqrt(grid.scale) / 2)
      ## note: 0.17 is a guess for yf (see plot.seqpcplot)
    }

    ## finalize point coordinates

    ## point widths
    pts[, wdtcols] <- data.frame(sqrt(scale(x = pts[, frqcols, drop = FALSE], center = FALSE, scale = wid.group.use)) * sqrt(grid.scale) * ngrid0 / ngrid)

    ## extract coordinate pairs for plotting the lines using
    ## segments (necessary to allow different line widths)
    ntraj <- max(pts$traj)
    lns <- NULL

    for (i in 1:ntraj) {
      SUBS <- which(pts$traj == i)
      nSUBS <- length(SUBS)
      if (length(SUBS) > 1) {
        TMP <- data.frame(traj = rep(i, nSUBS-1), x0 = pts$x[SUBS[-nSUBS]], y0 = pts$y[SUBS[-nSUBS]], x1 = pts$x[SUBS[-1]], y1 = pts$y[SUBS[-1]])
        lns <- rbind(lns,TMP)
      }
    }

    if (!is.null(lns)) {
      lns[, frqcols] <- 0
      lns[, lwdcols] <- 1

      ## line widths
      for (i in 1:nrow(lns)) { # loop over all lns entries
        if (i == 1) {
          SUBS <- traj.id[traj.group == lns$traj[i]]
        } else {
          if (lns$traj[i] != lns$traj[i-1]) {
            SUBS <- traj.id[traj.group == lns$traj[i]]
          }
        }
        lns[i, frqcols] <- unlist(tapply(X = weights.use[SUBS] * w[SUBS] * ((minx[SUBS] <= lns$x0[i]) & (maxx[SUBS] >= lns$x1[i])), INDEX = list(factor(group[SUBS], levels = 1:ngroup)), FUN = sum))
      }
      lns[is.na(lns)] <- 0

      ## line widths

      lns[, lwdcols] <- data.frame(sqrt(scale(x = lns[, frqcols, drop = FALSE], center = FALSE, scale = wid.group.use)) * sqrt(grid.scale) * ngrid0 / ngrid)

      ## merge coordinates and colors
      lns <- merge(x = lns, y = pts[, c("traj", "x", "y", "xjitter", "yjitter")], by.x = c("traj", "x0", "y0"), by.y = c("traj", "x", "y"), all.x = TRUE, all.y = FALSE, sort = FALSE)
      lns <- merge(x = lns, y = pts[,c("traj","x","y","xjitter","yjitter", colcols)], by.x = c("traj","x1","y1"), by.y = c("traj","x","y"), all.x = TRUE, all.y = FALSE, sort = FALSE)
      colnames(lns)[(ncol(lns) - 3 - ngroup):ncol(lns)] <- c("x0jitter", "y0jitter", "x1jitter", "y1jitter", colcols)
    }

    if (is.na(rows) & is.na(cols)) {
      if (ngroup > 1) {
        optpr <- optimpanelraster(nx = ngroup, ny = ngroup, npanels = ngroup, c = 1)
        nxl <- optpr[1]; nyl <- optpr[2]
      } else {
        nxl <- 1; nyl <- 1
      }
    } else if ((!is.na(rows)) & (!is.na(cols))) {
      nxl <- cols; nyl <- rows;
    } else if ((!is.na(rows)) & is.na(cols)) {
      nxl <- ceiling(ngroup / rows); nyl <- rows
    } else {
      nxl <- cols; nyl <- ceiling(ngroup / cols);
    }

    ## coordinates for background rectangles
    backrect <- expand.grid(xgrid = 1:nx, ygrid = 1:ny)

    ## reset seed
    assign(".Random.seed", oldSeed, envir=globalenv())
    
    ## plot data object

    x <- list(pts = pts, # pointdata
              lns = lns, # linedata
              backrect = backrect, # background grid
              propid.use = propid.use, propid.use.max = propid.use.max,
              propid.tot = propid.tot, ntraj = ntraj,
              ngroup = ngroup, which = 1:ngroup, # which group to plot
              nid.group.tot = nid.group.tot,
              nid.group.use = nid.group.use,
              wid.group.tot = wid.group.tot,
              wid.group.use = wid.group.use,
              frqcols = frqcols, wdtcols = wdtcols, colcols = colcols,
              lwdcols = lwdcols,
              use.layout = !add, nx = nx, ny = ny,
              xlevs = xlevs, ylevs = ylevs, grouplev = grouplev,
              title = title, mtext = mtext, xlab = xlab, ylab = ylab,
              xlim = xlim, ylim = ylim,
              nxl = nxl, nyl = nyl, # layout
              ngrid = ngrid, ngrid0 = ngrid0,
              grid.col = grid.col, grid.lwd = grid.lwd,
              grid.scale = grid.scale,
              grid.fill = grid.fill, grid.shape = grid.shape,
              grid.border = grid.border, hide.col = hide.col,
              cex = cex, border = border, border.lwd = border.lwd,
              lwd = lwd, lorder = lorder, col.nobs = col.nobs,
              sf.cex = sf.cex, sf.cex.leaves = sf.cex.leaves,
              add = add, xaxis = xaxis, yaxis = yaxis, axes = axes,
              cex.plot = cex.plot)
    class(x) <- "seqpcplot" # to allow replot

  }

  ## Step 7: Execute plot .................................. #

  if (plot) {
    if (verbose) cat(" [>] plotting\n")
    plot(x, ...)
  }

  ## Step 8: Return data ................................... #

  invisible(x)
}

plot.seqpcplot <- function(x, add = NULL, which = NULL, ...) {

  if (!is.null(add)) x$add <- add; rm(add)

  if (!is.null(which)) {
    if (is.character(which)) {
      x$which <- which(x$grouplev %in% which)
    } else {
      x$which <- which
    }
  }

  if (x$ngroup != length(x$which)) x$use.layout <- FALSE

  if (x$use.layout & !((x$nxl == 1) & (x$nyl == 1))) {
    layout(matrix(1:(x$nxl * x$nyl), ncol = x$nxl, nrow = x$nyl, byrow = TRUE))
  }

  par(cex.lab = x$cex.plot)
  plist <- list(x = 1, y = 1, xlab = x$xlab, ylab = x$ylab, xlim = x$xlim, ylim = x$ylim, type = "n", axes = FALSE)
  plist <- c(plist, list(...)[!names(list(...)) %in% names(plist)])

  for (i in x$which) {

    if (!x$add) {

      do.call(plot, args = plist)
      if (x$xaxis & (((x$axes == "all") & ((par("mar")[1] > 0) | (x$ngroup - i < x$nxl))) | ((x$axes == "bottom") & (x$ngroup - i < x$nxl)))) {
        axis(1, 1:x$nx, x$xlevs)
      }
      if (x$yaxis & ((par("mar")[2] > 0) | (((i-1) %% x$nxl) == 0))) {
        axis(2, 1:x$ny, x$ylevs, las = 2)
      }
      title(main = x$title[i])
      if (length(x$mtext) != 0) mtext(x$mtext[i], side = 3)
    }

    if (i == x$which[1]) {
        x$lns[, x$lwdcols] <- x$lns[, x$lwdcols]/max((par("cxy")/par("cin")))  * 96 * x$lwd * x$cex
      }

    ## prettify plot
    ppin <- par("pin")
    pusr <- par("usr")
    xf <- abs(pusr[2L] - pusr[1L]) / ppin[1L]
    yf <- abs(pusr[4L] - pusr[3L]) / ppin[2L]
    TMP <- max(xf, yf)
    xf <- xf/TMP; yf <- yf/TMP

    ## rectangles in the background
    if (x$grid.shape == "default") { # all points within rectangle
      rect(xleft = x$backrect$xgrid - xf * sqrt(x$grid.scale) / 2, ybottom = x$backrect$ygrid - yf * sqrt(x$grid.scale) / 2, xright = x$backrect$xgrid + xf * sqrt(x$grid.scale) / 2, ytop = x$backrect$ygrid + yf * sqrt(x$grid.scale) / 2, border = x$grid.border, col = x$grid.fill, lwd = x$grid.lwd)
      ## if (x$grid.lwd != 0) abline(h = 1:x$ny, v = 1:x$nx, col = x$grid.col, lwd = x$grid.lwd)
    } else { # rectangle is proportional to point sizes
      if (x$grid.shape == "proportional") {
        rect(xleft = x$backrect$xgrid - xf * sqrt(x$grid.scale) / 2, ybottom = x$backrect$ygrid - yf * sqrt(x$grid.scale) / 2, xright = x$backrect$xgrid + xf * sqrt(x$grid.scale) / 2, ytop = x$backrect$ygrid + yf * sqrt(x$grid.scale) / 2, border = x$grid.fill, col = NULL)
        rect(xleft = x$backrect$xgrid - xf * sqrt(x$grid.scale) / 2 * x$ngrid0 / x$ngrid * x$cex, ybottom = x$backrect$ygrid - yf * sqrt(x$grid.scale) / 2 * x$ngrid0 / x$ngrid * x$cex, xright = x$backrect$xgrid + xf * sqrt(x$grid.scale) / 2 * x$ngrid0 / x$ngrid * x$cex, ytop = x$backrect$ygrid + yf * sqrt(x$grid.scale) / 2 * x$ngrid0 / x$ngrid * x$cex, border = x$grid.fill, col = x$grid.fill)
      }
    }

    ## plot the subjects without observations
    if (x$propid.tot[x$ntraj + 1, i] > 0) {
      TMP <- x$pts$wdt.1[which.max(x$pts$n.1)] *
        sqrt(x$propid.tot[x$ntraj + 1, i] /
             (x$pts$n.1[which.max(x$pts$n.1)] / x$wid.group.tot[1]))
      rect(xleft = 1 - xf * (sqrt(x$grid.scale) / 2  + TMP * x$cex),
           ybottom = 1 - yf * (sqrt(x$grid.scale) / 2 + TMP * x$cex),
           xright = 1 - sqrt(x$grid.scale) / 2 * xf,
           ytop = 1 - sqrt(x$grid.scale) / 2 * yf,
           col = x$col.nobs[i], border = 0, lwd = 0)
    }

    for (j in order(x$propid.use[,i], decreasing = (x$lorder == "background"))) {
      SUBSp <- (x$pts$traj %in% j) & (x$pts[, x$frqcols[i]] > 0)
      SUBSd <- SUBSp & (x$pts$Freq > 1)
      if (!is.null(x$lns)) {
        SUBSl <- (x$lns$traj %in% j) & (x$lns[, x$frqcols[i]] > 0)
      } else {
        SUBSl <- FALSE
      }

      ## draw the points
      if (sum(SUBSp) > 0) {
        rect(xleft = (x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2) - c(x$pts[, x$wdtcols[i]]) / 2 * xf * x$cex)[SUBSp], ybottom = (x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2) - c(x$pts[, x$wdtcols[i]]) / 2 * yf * x$cex)[SUBSp], xright = (x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2) + c(x$pts[, x$wdtcols[i]]) / 2 * xf * x$cex)[SUBSp], ytop = (x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2) + c(x$pts[, x$wdtcols[i]]) / 2 * yf * x$cex)[SUBSp], col = x$pts[SUBSp, x$colcols[i]], border = if (!is.null(x$border)) x$border else x$pts[SUBSp, x$colcols[i]], lwd = x$border.lwd)
      }

      ## draw the lines
      if (sum(SUBSl) > 0) {
        segments(x0 = (x$lns$x0 + xf * (- 0.5 + x$lns$x0jitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSl], y0 = (x$lns$y0 + yf * (- 0.5 + x$lns$y0jitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSl], x1 = (x$lns$x1 + xf * (- 0.5 + x$lns$x0jitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSl], y1 = (x$lns$y1 + yf * (- 0.5 + x$lns$y0jitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[SUBSl], lwd = x$lns[SUBSl, x$lwdcols[i]], col = x$lns[SUBSl, x$colcols[i]], lend = 0)
      }

      ## draw the flowers
      if (sum(SUBSd) > 0) {
        points(x = (x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[SUBSd], y = (x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[SUBSd], pch = 16, cex = x$sf.cex)
        i.multi <- which(SUBSd) # stolen from sunflowerplot()
        ppin <- par("pin")
        pusr <- par("usr")
        xr <- x$pts[SUBSd, x$wdtcols[i]] / 2 * x$sf.cex.leaves * xf
        yr <- x$pts[SUBSd, x$wdtcols[i]] / 2 * x$sf.cex.leaves * yf
        i.rep <- rep.int(i.multi, x$pts$Freq[SUBSd])
        z <- numeric()
        for (k in i.multi) {
          z <- c(z, 1:x$pts$Freq[k])
        }
        deg <- (2 * pi * z)/x$pts$Freq[i.rep]
        segments(x0 = (x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep], y0 = (x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep], x1 = (x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep] + xr * sin(deg), y1 = (x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid.scale) + 1/2 * (1 - sqrt(x$grid.scale)) + c(sqrt(x$propid.use.max) * sqrt(x$grid.scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep] + yr * cos(deg))
      }

    }
    if (!x$add) {
      box()
    }
  }
}

## create ordering of events
ordering <- function(x, align = "first") {
  lev <- sort(unique(x), decreasing = FALSE)
  ret <- as.integer(factor(x, levels = lev))
  ret <- ret - ifelse(align == "last", max(ret), 0)
  return(ret)
}

## create trajectory strings
createstring <- function(y, x = 1:length(y),
                         pa.left = "(", pa.right = ")",
                         pot= "^" ,con = "-") {
  subscripts <- !is.na(y);
  ret <- paste(pa.left, y, pa.right, pot, x, sep = "");
  ret <- ret[subscripts];
  ret <- paste(ret, collapse = con)
  return(ret)
}

## layout optimization function
fmin <- function(nx, ny, n, c = 0.5) {
  control <- nx * ny - n # number of empty windows
  ## optimization value: function of number of empty windows
  ## and ncol/nrow ratio
  opt <- control / n + c * abs(1 - min(nx, ny) / max(nx, ny))
  return(c(control,opt))
}

## layout optimizer
optimpanelraster <- function(nx, ny, npanels,c = 1) {
  minvalue1 <- fmin(nx - 1, ny, npanels, c)
  minvalue2 <- fmin(nx, ny - 1, npanels, c)
  while (minvalue1[1] >= 0 | minvalue2[1] >= 0) {
    if ((minvalue1[1] >= 0) & (minvalue2[1] >= 0)) {
      if (minvalue1[2] < minvalue2[2]) {
        nx <- nx - 1 } else { ny <- ny - 1 }
    } else {
      if (minvalue1[2] >= 0) {
        nx <- nx - 1
      } else {
        if (minvalue2[2] >= 0) {
          ny <- ny - 1
        }
      }}
    minvalue1 <- fmin(nx - 1, ny, npanels, c)
    minvalue2 <- fmin(nx, ny - 1, npanels, c)
  }
  return(c(nx, ny))
}

## convert input for filter argument
construct.filter <- function(x) {

  if (is.numeric(x)) {
    x <- x[1L]
    stopifnot(x >= 0 & x <= 1)
    x <- list(type = "function", value = "minfreq", level = x)
  }
  
  if (is.list(x) | inherits(x, "seqpcfilter")) {
    if (sum(names(x) %in% c("type", "value")) != 2) stop("[!] filter must contain a type and a value object")
    if (!x$type %in% c("sequence", "subsequence", "value", "density", "function")) stop("[!] unknown input for filter$type")
    if ((x$type == "function") & (length(x) > 1)) {
      if (is.character(x$value)) {
        if (x$value == "linear") {
          x$value <- linear
        } else if (x$value == "minfreq") {
          x$value <- minfreq
        } else if (x$value == "cumfreq") {
          x$value <- cumfreq
        } else {
          stop("[!] invalid x$value argument")
        }
      }
      ret <- list(type = x$type, value = x$value, args = x[!names(x) %in% c("type", "value")])
    } else {
      ret <- x
    }
  }
  if (is.function(x)) {
    ret <- list(type = "function", value = x)
  }
  if (is.vector(x)) {
    if (is.numeric(x)) {
      ret <- list(type = "value", value = x)
    }
    if (is.character(x)) { # could either be sequence, subsequence or a predefined function

      if (x == "linear") {
        ret <- list(type = "function", value = linear)
      } else if (x == "minfreq") {
        ret <- list(type = "function", value = minfreq)
      } else if (x == "cumfreq") {
        ret <- list(type = "function", value = cumfreq)
      } else {
        ret <- list(type = "sequence", value = x)
      }
    }
  }
  if (is.matrix(x)) {
    ret <- list(type = "density", value = x)
  }
  return(ret)
}

## functions for colouring the lines
## ---------------------------------

## colourramp
colourize <- function(value, col1, col2) { # define colouring function
  mp <- colorRamp(c(col1, col2))
  col <- rgb(mp(value), maxColorValue = 255)
}

## convenience function
seqpcfilter <- function(method = c("minfreq", "cumfreq", "linear"), level = 0.05) {
  value <- match.arg(method)
  if (is.null(level) && method %in% c("minfreq", "cumfreq"))
    stop("'seqpcfilter' requires an inpute for 'level'.")
  return(structure(list(type = "function", value = value,
                        level = level), class = "seqpcfilter"))
}

## linear colour gradient function
linear <- function(x, level = NULL) {
  return((x - min(x)) / diff(range(x)))
}

## minimal frequency for lines to colour
minfreq <- function(x, level = 0.05) {
  return(1*(x >= level))
}

## colour a given proportion of most frequent sequences
cumfreq <- function(x, level = 0.75) {
  TMP <- which(cumsum(sort(x, decreasing = TRUE)) >= level)[1]
  ret <- vector("logical", length(x))
  ret[order(x, decreasing = TRUE)[1:TMP]] <- TRUE
  return(1 * ret)
}

rectangular <- function(x) {dunif(x, min = -0.5, max = 0.5)}

summary.seqpcplot <- function(object, ...) {

  ## get table
  ret <- as.data.frame(object$propid.tot)

  ## rownames
  TMP3 <- object$pts[order(object$pts$traj, object$pts$x, object$pts$y), ]
  rownames(ret) <- c(as.character(seqecreate(timestamp = TMP3$x, event = factor(TMP3$y, 1:object$ny, object$ylevs), id = TMP3$traj)), "1-()") # "hidden"

  ## colnames
  if (object$ngroup > 1) {
    colnames(ret) <- paste("prop", object$grouplev, sep = ".")
  } else {
    colnames(ret) <- "prop"
  }

  ## return
  return(ret)
}
