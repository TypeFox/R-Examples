nam <- function (input, rounding = TRUE) {
    mt <- ifelse(input$sm > 0, 1, 0)
    stopifnot(is.matrix(mt))
    if (is.null(class(input)) | !(class(input) %in% c("dotinference", "gridinference"))) {
        cat("Input object is not of any class of inference \n")
        return(invisible())
        }
    stopifnot(is.logical(rounding))
    sp_num <- nrow(mt)
    stopifnot(isSymmetric(mt))
    edge <- vector("list", sp_num)
    compdn <- array(0, sp_num)
    diag(mt) <- 0
    degree <- apply(mt, 1, sum)
    for (i in 1:sp_num) edge[[i]] <- which(mt[i,] > 0)
    betscore <- NULL
    nb <- NULL
    recorded <- array(0, sp_num)
    modificado <- ignore <- array(FALSE, sp_num)
    recalcular <- 0
    component <- array(0, sp_num)
    total <- array(1, sp_num)
    betwns <- function(IDsp, edge_bet, degree) {
        if (degree[IDsp] == 0) {
            total[IDsp] <<- 2
            component[IDsp] <<- 1
            return()
        }
        pred <- matrix(0, nrow = sp_num, ncol = sp_num)
        d <- array(-1, sp_num)
        p <- array(0, sp_num)
        npred <- array(0, sp_num)
        q <- IDsp
        qpush <- 1
        qpop <- 0
        d[IDsp] <- 0
        p[IDsp] <- 1
        orden <- IDsp
        count <- 1
        while (qpush > qpop) {
            qpop <- qpop + 1
            v1 <- q[qpop]
            d2 <- d[v1] + 1
            for (v2 in edge_bet[[v1]]) {
                if (!ignore[v2]) {
                  if (d[v2] == -1) { #Not visited yet
                    d[v2] <- d2
                    p[v2] <- p[v1]
                    npred[v2] <- npred[v2] + 1
                    pred[v2, npred[v2]] <- v1
                    q <- c(q, v2)
                    qpush <- qpush + 1
                    count <- count + 1
                    orden <- c(orden, v2)
                  }
                  else if (d2 == d[v2]) {
                    p[v2] <- p[v2] + p[v1]
                    npred[v2] <- npred[v2] + 1
                    pred[v2, npred[v2]] <- v1
                  }
                }
            }
        }
        caminos <- array(0, sp_num)
        for (i in 1:count) caminos[orden[i]] <- 1
        if (count > 1) {
            for (i in count:2) {
                v1 <- orden[i]
                for (j in 1:npred[v1]) {
                  v2 <- pred[v1, j]
                  caminos[v2] <- caminos[v2] + caminos[v1] *
                    p[v2]/p[v1]
                }
            }
        }
        total[orden] <<- total[orden] + caminos[orden]
        component[IDsp] <<- count
        compdn[orden] <<- compdn[orden] + 1
    }
    msg <- "Initializing ..."
    cat(msg)
    flush.console()
    while (1) {
        component[compdn %in% recalcular] <- 0
        total[compdn %in% recalcular] <- 1
        modificado[compdn %in% recalcular] <- TRUE
        for (xp in 1:sp_num) {
            if (!ignore[xp] && modificado[xp]){
                betwns(xp, edge, degree)
            }
            backspaces <- paste(rep("\b", nchar(msg) + 1), collapse = "")
            msg <- paste(100*round(xp/sp_num, 2)," % at subnetwork ", max(recorded))
            cat(backspaces, msg)
            flush.console()
        }
        if(rounding) bvector <- zapsmall(total/2 - component, digits = 8)
        else bvector <- total/2 - component
        bvector[ignore] <- -1
        hb <- max(bvector)
        if (hb > 1e-10) {
            betscore <- c(betscore, hb)
            sp <- which(bvector == hb)
            nb <- c(nb, sp)
            recorded[-nb] <- recorded[-nb] + 1
            ignore[sp] <- TRUE
            modificado[] <- FALSE
            recalcular <- unique(compdn[sp])
            }
        else break
    }
    backspaces <- paste(rep("\b", nchar(msg) + 1), collapse = "")
    msg <- "NAM finished\n\nPlease wait, cleavogram is being created...\n\n"
    cat(backspaces, msg)
    flush.console()

    ##The next lines construct the cleavogram
    ##and calculate some indices of cohesiveness
    xsubnet <- list() # List of nodes by sub-network
    for (i in seq(0, max(recorded))) xsubnet <- c(xsubnet, list(which(recorded >= i)))
    nsub <- max(recorded) + 1
    full <- c(1:sp_num)
    numunits <- 0
    diads <- 0
    remid <- sapply(xsubnet, function(x) setdiff(full, x))
    remunq <- paste("Rem", full)
    isounq <- paste("Iso", full)
    code_sp <- array("Removed", sp_num)
    status <- array(dim = sp_num)
    analyzable <- array(FALSE, sp_num)
    groups <- c()
    for (i in xsubnet) {
      analyzable[i] <- TRUE
      while (any(analyzable)) {
        C <- Q <- match(TRUE, analyzable)
        qpush <- 1
        analyzable[Q] <- FALSE
        while (qpush > 0) {
            v1 <- Q[1]
            if (length(edge[[v1]]) == 0) {
                qpush <- 0
                next
            }
            for (i in 1:length(edge[[v1]])) {
                if (analyzable[edge[[v1]][i]]) {
                  C <- c(C, edge[[v1]][i])
                  Q <- c(Q, edge[[v1]][i])
                  qpush <- qpush + 1
                  analyzable[edge[[v1]][i]] <- FALSE
                }
            }
            Q <- Q[-1]
            qpush <- qpush - 1
        }
        if (length(C) > 2) {
           numunits <- numunits + 1
           code_sp[C] <- paste("G", numunits, sep = "")
           } else if (length(C) == 2) {
            diads <- diads + 1
            code_sp[C] <- paste("D", diads, sep = "")
           } else code_sp[C] <- "Isolated"
      }
      analyzable[] <- FALSE
      groups <- cbind(groups, code_sp)
      code_sp[] <- "Removed"
    }
    colnames(groups) <- paste("Net", seq(0, max(recorded)), sep = "")
    groups <- cbind(groups,c(1:sp_num))
    groups <- groups[do.call(order, data.frame(groups)),]
    namlast <- apply(groups, 1, function(x) match("Removed", x))
    namlast <- replace(namlast, is.na(namlast), nsub + 1)

    #Individualize each either removed or isolated element
    for(i in 1:nsub) {
      which(groups[,i] == "Removed") -> indrem
      which(groups[,i] == "Isolated") -> indiso
      groups[indrem,i] <- remunq[indrem]
      groups[indiso,i] <- isounq[indiso]
    }

    checkpost <- function(aux, x1){
      if (aux[x1] == aux[x1 + 1]) checkpost(aux, x1 + 1) else return(x1 + 1)
    }

    #Compute pairwise geodesic distances between nodes
    APD <- function(A, n = nrow(A)) {
       Z <- A%*%A
       degree <- diag(Z)
       B <- ifelse((A + Z) > 0, 1, 0)
       diag(B) <- 0
       if (all(B[lower.tri(B)]==1)) return (2*B-A)
       T <- APD(B)
       X <- T%*%A
       D <- 2*T
       for (i in 1:(n -1))
         for (j in (i + 1):n)
         if (X[i,j] < T[i,j]*degree[j]) D[i,j] <- D[j,i] <- D[i,j] - 1
       diag(D) <- 0
       return (D)
    }

    #Add some auxilliary functions
    cohesive <- function(matriz) {
      spnum <- nrow(matriz)
      auxmt <- (matriz%*%matriz)[lower.tri(matriz)] #>0 for pairs sharing neighbors
      eccen <- apply(APD(matriz, spnum), 1, max)
      ds <- mean(matriz[lower.tri(matriz)]) #density
      tr <- sum(matriz[lower.tri(matriz)]*auxmt/sum(auxmt)) #global clustering or transitivity
      di <- c(range(eccen), mean(eccen)) #eccentricity
      return(c(ds, tr, di))
    }

    #Graphical information
    verticals <- c()
    labelcoord <- c()
    components <- c()
    grname <- c()
    trunk <- array(FALSE, sp_num)
    lhor2 <- sp_num
    nlhor2 <- 1
    until <- c()
    whenadded <- array(-1, sp_num)
    l <- 0
    firstdiv <- TRUE
    prec <- array(-1, sp_num)
    for (i in 1:nsub) {
      aux <- c(groups[,i], "end")
      lhor1 <- c() # Number of nodes involved in each branch
      init <- 1
      while(aux[init] != "end") {
        finit <- checkpost(aux, init) - 1
        lhor1 <- c(lhor1, finit)
        if(prec[finit] == init) { # The group has previously been detected
           init <- finit + 1
           if (namlast[finit] < i) next # Intermediary species
           until[whenadded[finit]] <- until[whenadded[finit]] + 1
           next
        }
        #Indicates that the group has been partitioned
        if(firstdiv) {f <- 0.5*(init + finit); firstdiv <- FALSE}
        if(finit %in% lhor2) {verticals <- rbind(verticals, c(i, f, i, 0.5*(init + finit))); firstdiv <- TRUE}
        l <- l + 1
        prec[finit] <- init
        whenadded[finit] <- l
        if (init == finit) {
          components <- rbind(components, c(rep(NA, 5), i, init, finit))
          grname <- c(grname, groups[init,i])
          until <- c(until, i)
          init <- finit + 1
          next
        }
        if ((finit - init) == 1) {
          components <- rbind(components, c(c(1, 0, rep(1, 3)), i, init, finit)) #Transitivity is set to 0
          grname <- c(grname, groups[init,i])
          until <- c(until, i)
          init <- finit + 1
          next
        }
        spp <- as.integer(groups[init:finit,(nsub + 1)])
        measures <- cohesive(mt[spp,spp])
        components <- rbind(components, c(round(measures, 5), i, init, finit))
        grname <- c(grname, groups[init,i])
        until <- c(until, i)
        init <- finit + 1
      } #end while
      lhor2 <- lhor1
   } #end for
   until[grep("Rem",components[,1])] <- until[grep("Rem",components[,1])] - 0.5
   components <- cbind(components, until)
   components <- data.frame(components, row.names = grname)
   colnames(components) <- c("Density", "Transitivity", "MinEccen", "MaxEccen", "MeanEccen",
                             "FirstNet", "FromNode", "ToNode", "RefNet")
   newor <- as.integer(groups[,nsub + 1])
   cat("Done and enjoy if possible!\n")
   out <- list(mt = mt[newor, newor], LastNet = recorded[newor], namlast = namlast, Betweenness = betscore,
               leaves = names(input$occupancy)[newor], nsp = sp_num, nsub = nsub, components = components, verticals = verticals,
               kind = input$kind, occupancy = input$occupancy, coords = input$coords)
   class(out) <- "cleavogram"
   return(out)
}



