makeS <- function(pedigree, heterogametic, DosageComp = c(NULL, "ngdc", "hori", "hedo", "hoha", "hopi"), returnS = FALSE){

    if(length(unique(pedigree[,4])) > 2) stop("Error: more than 2 sexes specified")
    nPed <- numPed(pedigree[, 1:3])

    damsex <- pedigree[unique(nPed[, 2])[-1], 4]
    if(any(damsex == heterogametic)){
       pedname <- names(pedigree)
       pedigree <- pedigree[, c(1,3,2,4)]
       names(pedigree) <- pedname
       nPed <- numPed(pedigree[, 1:3])
      cat("Assuming female heterogametic (e.g., ZZ/ZW) sex chromosome system\n")
    } else cat("Assuming male heterogametic (e.g., XX/XY) sex chromosome system\n")

    sex <- rep(-998, dim(pedigree)[1])
    sex[homs <- which(pedigree[,4] != heterogametic)] <- 1
    sex[hets <- which(pedigree[,4] == heterogametic)] <- 0
    N <- dim(nPed)[1]
    N2 <- N + 1

    dc.model <- match.arg(DosageComp)
    if(is.null(dc.model)){
      warning("Assuming 'ngdc' dosage compensation model")
      dc.model <- "ngdc"
    }


    if(dc.model == "ngdc"){
          dnmiss <- which(nPed[,2] != -998)
          fsnmiss <- which(nPed[,3] != -998 & sex == 1)
          bnmiss <- which(nPed[, 2] != -998 & nPed[, 3] != -998)
          nA <- N + 2 * length(dnmiss) + 2 * length(fsnmiss)
          nA <- nA + 2 * sum(duplicated(paste(nPed[, 2], nPed[, 3])[bnmiss]) == FALSE)
          Q.col <- c(nPed[,1][dnmiss], nPed[,1][fsnmiss], 1:N) 
          Q.row <- c(nPed[,2][dnmiss], nPed[,3][fsnmiss], 1:N)
          Q.x <- c(rep(-0.5, length(dnmiss)), rep(-1, length(fsnmiss)), rep(1, N))
          ord <- order(Q.col + Q.row/(N+1), decreasing = FALSE)
          Q <- Matrix(0, N, N, sparse = TRUE)
          Q[1, 2] <- 1
          Q@i <- as.integer(Q.row[ord] -1)
          Q@p <- as.integer(c(match(1:N, Q.col[ord]), length(ord) + 1) - 1)
          Q@x <- as.double(Q.x[ord])

          nPed[nPed == -998] <- N2
          Vii <- (sex + 1)/2
          f <- c(rep(0, N), -1)

          Cout <- .C("sinv",
	    as.integer(nPed[, 2] - 1), #dam
	    as.integer(nPed[, 3] - 1),  #sire
	    as.double(f),  #f
	    as.double(Vii),  #vii
            as.integer(Q@i),  #iQP
	    as.integer(c(Q@p, length(Q@x))),  #pQP
            as.double(Q@x),  #xQP
            as.integer(N),   #nQP
	    as.integer(length(Q@x)), #nzmaxQP	
	    as.integer(rep(0, nA)), #iSP
	    as.integer(rep(0, N + 1)), #pSP
	    as.double(rep(0, nA)), #xSP
	    as.integer(nA), #nzmaxSP
            as.double(0.25), #DC
            as.double(Vii))  #sex

          Sinv <- Matrix(0, N, N, dimnames = list(as.character(pedigree[,1]), NULL))
          Sinv[1, 2] <- 1
          Sinv@i <- Cout[[10]][1:Cout[[13]]]
          Sinv@p <- Cout[[11]]
          Sinv@x <- Cout[[12]][1:Cout[[13]]]
          Vii <- Cout[[4]]
          f <- Cout[[3]][1:N]



    } else{
         if(dc.model != "hopi"){
             fdnmiss <- which(nPed[,2] != -998 & sex == 1)
             mdnmiss <- which(nPed[,2] != -998 & sex == 0)
             fsnmiss <- which(nPed[,3] != -998 & sex == 1)
             bnmiss <- which(nPed[, 2] != -998 & nPed[, 3] != -998)
             nA <- N + 2 * (length(fdnmiss) + length(mdnmiss)) + 2 * length(fsnmiss)
             nA <- nA + 2 * sum(duplicated(paste(nPed[, 2], nPed[, 3])[bnmiss]) == FALSE)
             Q.col <- c(nPed[,1][fdnmiss], nPed[,1][mdnmiss], nPed[,1][fsnmiss], 1:N) 
             Q.row <- c(nPed[,2][fdnmiss], nPed[,2][mdnmiss], nPed[,3][fsnmiss], 1:N)
             Q.x <- c(rep(-0.5, length(fdnmiss)), rep(-1, length(mdnmiss)), rep(-0.5, length(fsnmiss)), rep(1, N))
             ord <- order(Q.col + Q.row/(N+1), decreasing = FALSE)
             Q <- Matrix(0, N, N, sparse = TRUE)
             Q[1, 2] <- 1
             Q@i <- as.integer(Q.row[ord] - 1)
             Q@p <- as.integer(c(match(1:N, Q.col[ord]), length(ord) + 1) - 1)
             Q@x <- as.double(Q.x[ord])

             nPed[nPed == -998] <- N2
             Vii <- (2 - sex)
             f <- c(rep(0, N), -1)

             Cout <- .C("sinv",
	       as.integer(nPed[, 2] - 1), #dam
	       as.integer(nPed[, 3] - 1),  #sire
	       as.double(f),  #f
	       as.double(Vii),  #vii
               as.integer(Q@i),  #iQP
	       as.integer(c(Q@p, length(Q@x))),  #pQP
               as.double(Q@x),  #xQP
               as.integer(N),   #nQP
	       as.integer(length(Q@x)), #nzmaxQP	
	       as.integer(rep(0, nA)), #iSP
	       as.integer(rep(0, N + 1)), #pSP
	       as.double(rep(0, nA)), #xSP
	       as.integer(nA), #nzmaxSP
               as.double(1), #DC
               as.double(Vii))  #sex

             Sinv <- Matrix(0, N, N, dimnames = list(as.character(pedigree[,1]), NULL))
             Sinv[1, 2] <- 1
             Sinv@i <- Cout[[10]][1:Cout[[13]]]
             Sinv@p <- Cout[[11]]
             Sinv@x <- Cout[[12]][1:Cout[[13]]]
             Vii <- Cout[[4]]
             f <- Cout[[3]][1:N]


         } else{
                dnmiss <- which(nPed[,2] != -998)
                nA <- N + 2 * length(dnmiss)
                nA <- nA + 2 * sum(duplicated(paste(nPed[, 2], nPed[, 3])[dnmiss]) == FALSE)
                Q.col <- c(nPed[,1][dnmiss], 1:N) 
                Q.row <- c(nPed[,2][dnmiss], 1:N)
                Q.x <- c(rep(-0.5, length(dnmiss)), rep(1, N))
                ord <- order(Q.col + Q.row/(N+1), decreasing = FALSE)
                Q <- Matrix(0, N, N, sparse = TRUE)
                Q[1, 2] <- 1
                Q@i <- as.integer(Q.row[ord] - 1)
                Q@p <- as.integer(c(match(1:N, Q.col[ord]), length(ord) + 1) - 1)
                Q@x <- as.double(Q.x[ord])

                nPed[nPed == -998] <- N2
                Vii <- rep(1, N) 
                f <- rep(0, N)

             Cout <- .C("sinv",
	       as.integer(nPed[, 2] - 1), #dam
	       as.integer(nPed[, 3] - 1),  #sire
	       as.double(f),  #f
	       as.double(Vii),  #vii
               as.integer(Q@i),  #iQP
	       as.integer(c(Q@p, length(Q@x))),  #pQP
               as.double(Q@x),  #xQP
               as.integer(N),   #nQP
	       as.integer(length(Q@x)), #nzmaxQP	
	       as.integer(rep(0, nA)), #iSP
	       as.integer(rep(0, N + 1)), #pSP
	       as.double(rep(0, nA)), #xSP
	       as.integer(nA), #nzmaxSP
               as.double(0), #DC
               as.double(Vii))  #sex

             Sinv <- Matrix(0, N, N, dimnames = list(as.character(pedigree[,1]), NULL))
             Sinv[1, 2] <- 1
             Sinv@i <- Cout[[10]][1:Cout[[13]]]
             Sinv@p <- Cout[[11]]
             Sinv@x <- Cout[[12]][1:Cout[[13]]]
             Vii <- Cout[[4]]
             f <- Cout[[3]][1:N]

           }
      }              


    listSinv <- sm2list(Sinv, rownames = as.character(pedigree[, 1]), colnames = c("Row", "Column", "Sinverse"))
    if(returnS){
       cat("S-inverse made: Starting to make S...")
          T <- as(solve(Q), "dgCMatrix")
          S <- as(t(T) %*% Diagonal(N, Vii) %*% T, "dgCMatrix")
       cat(".done", "\n")
    } else{
         S <- NULL
      }

return(list(model = dc.model, S = S, Sinv = Sinv, listSinv = listSinv, inbreeding = f, v = Vii))
}

