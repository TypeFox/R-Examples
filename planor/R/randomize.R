# 'planor.randomize'
# A function to randomize a factorial design according to a
# specified block structure formula
# ARGUMENTS
# - blockformula: the block structure formula
# - data: a data frame
# - out.order: a list of column names of data to order the result
# - keep.initial: if TRUE, the initial row order of the design
#   is stored in column 'InitialUNITS' of the returned dataframe
# RETURN
# the input data frame after randomization
# NOTE
#   Each name in 'blockformula' must correspond to a factor
#   of the dataframe 'data'. The only exception is
#   'UNITS'. If 'UNITS' is used in 'blockformula'
#   but absent from 'data', a factor is added to 'data',
#   with one level per row. See the examples below for the usage of
#   'UNITS' in 'blockformula'.
# REFERENCES
#   Bailey, R.A., 1983. Generalized wreath products of permutation groups. Proc. London Math. Soc., 47, 69-82.
#
#   Kobilinsky A., 1989. Randomization of a cartesian block structure. Technical Report. Laboratoire de Biométrie de l'INRA Versailles.
# EXAMPLES
#   ## Block design
#   Design <- data.frame(block=rep(1:4,rep(2,4)),
#                treatment=c("A1","B1","A2","B2","A3","B3","A4","B4"))
#   planor.randomize(~block, data=Design)       ##  no within-block randomization
#   planor.randomize(~block/UNITS, data=Design) ##  blocks and units within blocks randomization
#   ## Row-Column design
#   RowColDes <- data.frame(row=rep(1:3,rep(3,3)),col=rep(1:3,3),
#            treatment=LETTERS[c(1:3,2,3,1,3,1,2)],
#            oldRow=rep(1:3,rep(3,3)),oldCol=rep(1:3,3))
#    planor.randomize(~row*col, data=RowColDes)
# --------------------------------------------------------------
planor.randomize <- function(blockformula, data, out.order, keep.initial=FALSE){
    ## PRELIMINARIES
    ## M : matrix of block factorial terms
    M <- 1*(attr(terms(blockformula),"factors") > 0)
    IdF <- rownames(M)
    NbF <- nrow(M)
    NbT <- ncol(M)
    ## calculate weights, reorder
    Weights <- apply(M>0, 2, sum)
    M <- M[,order(Weights), drop=FALSE]
    Weights <- Weights[order(Weights)]
    ## give id numbers to the factorial terms
    IdCalc <- function(x, n=length(x)){sum(x*(2^(0:(n-1))))}
    IdM <- apply(M, 2, IdCalc)

    ## 1. CLOSURE OF THE BLOCK TERMS FOR THE INTERSECTION
    S <- M[,1,drop=FALSE]
    IdS <- IdM[1]
    if(NbT > 1){
        for(i in 2:NbT){
            for(j in 1:(i-1)){
                margin.ij <- M[,i] & S[,j]
                margin.id <- IdCalc(margin.ij)
                if(!(margin.id %in% c(0,IdS))){
                    S <- cbind(S,margin.ij)
                    IdS <- c(IdS, IdCalc(margin.ij))
                }
            }
            S <- cbind(S,M[,i])
            IdS <- c(IdS, IdCalc(M[,i]))
        }
    }
    ## reorder
    Weights <- apply(S>0, 2, sum)
    S <- S[,order(Weights), drop=FALSE]
    Weights <- Weights[order(Weights)]
    IdS <- IdS[order(Weights)]

    ## 2. POSET BLOCK STRUCTURE (SEE PROPOSITION 6, Kobi 89)
    NbS <- ncol(S)
    C <- S[,1,drop=FALSE]
    IdC <- IdS[1]
    Cunion <- C
    if(NbS > 1){
        for(k in 2:NbS){
            C <- cbind(C, S[,k]&(!Cunion))
            Cunion <- Cunion | C[,k]
        }

    }
    nonvoid <- apply(C, 2, sum) > 0
    S <- S[,nonvoid,drop=FALSE]
    C <- C[,nonvoid,drop=FALSE]

    ## 3. RANDOMIZE
    N <- nrow(data)
    RowPerm <- seq(N)
    NbR <- sum(nonvoid)
    ## add temporarily a UNITS column, if needed
    if(("UNITS" %in% IdF)&!("UNITS" %in% colnames(data))){
        data <- cbind(UNITS=factor(seq(N)), data)
        add.UNITS <- TRUE
    }
    else{ add.UNITS <- FALSE }
    ## loop on the randomization strata, from the lowest (units) to the highest
    for(i in rev(seq(NbR))){
        ## nesting factors
        FsupId <- IdF[ S[,i]&(!C[,i]) ]
        if( length(FsupId)==0 ){ F.sup <- factor(rep(1,N)) }
        else{                    F.sup <- interaction(data[FsupId]) }
        ## randomized factor(s)
        FrandId <- IdF[ C[,i]==1 ]
        F.rand  <- interaction(data[FrandId])
        ## loop on the level combinations of the nesting factors
        for(j in levels(F.sup)){
            select <- seq(N)[F.sup == j]
            Lev <- unique(F.rand[select])
            ## check that there are several levels of F.rand when F.sup == j
            ## (to do: check that before entering the loop on j)
            if(length(Lev) > 1){
                LevR <- Lev[ sample(length(Lev)) ]
                oldRowPerm <- RowPerm
                ## loop on the levels of F.rand, the randomized factor(s)
                for(k in seq(Lev)){
                    Lev.Rows  <- select[F.rand[select] == Lev[k] ]
                    LevR.Rows <- select[F.rand[select] == LevR[k]]
                    RowPerm[ Lev.Rows ] <- oldRowPerm[ LevR.Rows ]
                }
            }
        }
    }
    ## application of the random permutation
    col.rand <- ! (colnames(data) %in% IdF)
    data[,col.rand] <- data[RowPerm,col.rand, drop=FALSE]

    ## FINAL STEPS : Tidying up the final result
    ## reorder the rows
    if(missing(out.order)){ out.order <- IdF }
    data <- data[ do.call(order, data[out.order]) , ]
    ## manage data columns
    if(add.UNITS){ data <- data[, -1, drop=FALSE] }
    ## keep info on the permutation
    if(keep.initial){ data <- cbind(InitialUNITS=RowPerm, data) }
    ##
    return(data)
}
