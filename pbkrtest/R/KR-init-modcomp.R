KRmodcomp_init <- function(m1, m2, matrixOK=FALSE){
    UseMethod("KRmodcomp_init")
}

KRmodcomp_init.lmerMod <-
    KRmodcomp_init.mer <-
    function(m1, m2, matrixOK=FALSE) {
        ##comparison of the mean structures of the models
        ## it is  tested for that (1) m1 is mer and (2) m2 is either mer or a matrix
        mers<- if ( .is.lmm(m1) &
                    (.is.lmm(m2) | is.matrix(m2) ) )
                   TRUE
               else
                   FALSE
        
        if (!mers) {
            cat("Error in modcomp_init\n")
            cat(paste("either model ",substitute(m1), 
                      "\n is not a linear mixed of class mer(CRAN) or lmerMod (GitHub)\n \n",sep=' '))
            cat(paste("or model ", substitute(m2),"\n is neither of that class nor a matrix",sep=''))
            stop()
        }
    
        ##checking matrixcOK is FALSE but m2 is a matrix
        if (!matrixOK & is.matrix(m2)) {
            cat ('Error in modcomp_init \n')
            cat (paste('matrixOK =FALSE but the second model: ', substitute(m2),
                       '\n is  specified via a restriction matrix \n \n',sep=''))
            stop()
        }
        
        Xlarge <- getME(m1, "X")
        rlarge <- rankMatrix(Xlarge)
        ##code <- if ('mer' %in% class(m2)) {
        ##code <- if ('lmerMod' %in% class(m2)) {
        code <- if (.is.lmm(m2)){
            Xsmall <- getME(m2,"X")
            rsmall <- rankMatrix(Xsmall)
            rboth  <- rankMatrix(cbind(Xlarge,Xsmall))
            if (rboth == pmax(rlarge,rsmall)) {
                if (rsmall< rlarge) {
                    1
                } else {
                    if (rsmall > rlarge) {
                        0
                    } else {
                        -1
                    }
                }
            } else {
                -1
            }
        } else {
            ##now model m2  is a restriction matrix
            if (rankMatrix(rbind(Xlarge,m2)) > rlarge) {
                -1
            } else {
                1
            }
        }
        code
    }


