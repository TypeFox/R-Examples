state_age <-
function (Bdata, age, ID) {
    if (missing(age))
        stop("No age")
    if (missing(ID))
        stop("No ID")
    if (length(ID) > nrow(Bdata))
        ID = Bdata$ID
    Bdata2 <- Bdata[Bdata$ID %in% ID, ]
    namstates <- attr(Bdata, "param")$namstates
    iagelow <- attr(Bdata, "param")$iagelow
    iagehigh <- attr(Bdata, "param")$iagehigh
    agetrans <- AgeTrans(Bdata2)
    state <- matrix(data = 0, nrow = length(ID), ncol = length(age), dimnames = list(ID = ID, age = age))
    fillStateMatrixAlongID <- function(i){
      if (length(na.omit(agetrans$ages[i, ])) == 0)
        agecens <- agetrans$agecens[i]
      else {
        if (max(agetrans$ages[i, ], na.rm = TRUE) == agetrans$agecens[i])
          agecens <- agetrans$agecens[i] + 1e-06
        else agecens <- agetrans$agecens[i]
      }
      agelist.z <- unname(c(-1, agetrans$ageentry[i], agetrans$ages[i,], agecens, 10000))
      agelist <- sort(agelist.z)
      str_char <- stringf(Bdata2$path[i])
      d <- c("-", str_char, "+", "+")
      fillStateMatrixAlongAge <- function(j){
        if (j == 1)
         state[i, j] <<- d[which.max(agelist[age[j] >= agelist])]
        else state[i, j] <<- d[which.max(agelist[age[j] > agelist])]
        return(NULL)
      }
      nn <- sapply(1:length(age),fillStateMatrixAlongAge)
      return(NULL)
    }
    nn <- sapply(1:length(ID),fillStateMatrixAlongID)
    nam <- c("-", namstates, "+")
    state.n <- matrix(data = 0, nrow = length(age), ncol = length(nam) + 1, dimnames = list(age = age, nam = c(nam, "Total")))
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
    return(list(nam = nam, state = state, state.n = state.n))
}
