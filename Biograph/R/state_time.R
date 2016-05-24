state_time <-
function (Bdata, ID){
    if (FALSE %in% (ID %in% Bdata$ID))
        stop("state_time: ID contains identification numbers that are not part of the data.")
    if (length(ID) > nrow(Bdata))
        ID = Bdata$ID
    Bdata2 <- Bdata[Bdata$ID %in% ID, ]
    nsample <- length(ID)
    if (!exists("namstates"))
        param <- Parameters(Bdata)
    namstates <- attr(Bdata, "param")$namstates
    numstates <- length(namstates)
    iagelow <- attr(Bdata, "param")$iagelow
    iagehigh <- attr(Bdata, "param")$iagehigh
    agetrans <- AgeTrans(Bdata2)
    nam <- c("-", namstates, "+")
    sjt <- array(0, c(nsample, iagehigh - iagelow + 1, numstates + 2),
      dimnames = list(ID = Bdata2$ID, Age = iagelow:iagehigh, State = nam))
    state <- matrix(data = 0, nrow = nsample, ncol = iagehigh - iagelow + 1,
      dimnames = list(ID = ID, age = iagelow:iagehigh))
    fillStateMatrixAlongID <- function(i){
        agelist.z <- unname(c(iagelow:iagehigh, agetrans$ageentry[i], agetrans$ages[i, ], agetrans$agecens[i]))
        agelist <- unique(sort(agelist.z))
        st.char <- state_age(Bdata2, age = agelist, ID = Bdata2$ID[i])$state
        state[i, ] <<- st.char[as.numeric(colnames(st.char)) %in% c(iagelow:iagehigh)]
        st <- match(st.char, nam)
        fillStateMatrixAlongAge <- function(ix){
        sjt[i, trunc(agelist[ix]) - iagelow + 1, st[ix + 1]] <<-
            sjt[i, trunc(agelist[ix]) - iagelow + 1, st[ix + 1]] + agelist[ix + 1] - agelist[ix]
          return(NULL)
        }
        nn <- sapply(1:length(agelist) - 1,fillStateMatrixAlongAge)
        ix <- length(agelist)
        sjt[i, trunc(agelist[ix]) - iagelow + 1, st[ix]] <<-
          sjt[i, trunc(agelist[ix]) - iagelow + 1, st[ix]] + trunc(agelist[ix]) + 1 - agelist[ix]
        return(NULL)
    }
    nn <- sapply(1:length(ID),fillStateMatrixAlongID)
    sjt_age_1 <- sjt
    tsjt <- matrix(nrow = attr(Bdata, "param")$nage, ncol = length(nam) + 1,
        dimnames = list(Age = attr(Bdata, "param")$namage, Case = c(nam, "Total")))
    tsjt[, 1:length(nam)] <- apply(sjt_age_1, c(2, 3), sum)
    tsjt[, (length(nam) + 1)] <- apply(tsjt[, 1:length(nam)], 1, sum)
    print("state_time: compute state.n")
    state.n <- matrix(data = 0, nrow = length(iagelow:iagehigh), ncol = length(nam) + 1,
      dimnames = list(age = iagelow:iagehigh, nam = c(nam, "Total")))
    list.state <- apply(state, 2, table)
    itt.list <- 1
    fillStateNMatrix <- function(leAge){
      if(length(na.omit(leAge))>0) {
        ri <- which(rownames(state.n)==names(list.state)[itt.list])
        ci <- match(names(leAge),colnames(state.n))
        state.n[ri,ci] <<- leAge
      }
      itt.list <<- itt.list + 1
      return(NULL)
    }
    nn <- lapply(list.state,fillStateNMatrix)
    state.n[, (length(nam) + 1)] <- apply(state.n, 1, sum)
    return(list(state = state, state.n = state.n, sjt_age_1 = sjt_age_1,
        tsjt = tsjt))
}
