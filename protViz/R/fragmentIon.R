#R
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/fragmentIon.R $
# $Id: fragmentIon.R 6222 2014-03-13 14:22:34Z cpanse $
# $Date: 2014-03-13 15:22:34 +0100 (Thu, 13 Mar 2014) $


defaultIon<-function(b, y){
    Hydrogen <- 1.007825
    Oxygen <- 15.994915
    Nitrogen <- 14.003074

    #yo <- fi$y - Oxygen - Hydrogen - Hydrogen
    c<- b + (Nitrogen + (3 * Hydrogen))
    z<- y - (Nitrogen + (3 * Hydrogen))
    
    return(cbind(b, y, c ,z))
}


fragmentIon <- function(sequence, FUN=defaultIon, modified=numeric(), modification=numeric()) {
    if (!is.character(sequence)) {

        R<-list()

        input.n <- length(sequence)

        out <- .C("_computeFragmentIons", n=input.n, 
            W_=as.double(sequence), 
            b_=as.double(rep(0,input.n)), 
            y_=as.double(rep(0,input.n)))

        R[[1]] <- as.data.frame(FUN(out$b, out$y))

    }else if (length(modification) > 1) {
                                                    
        FUN<-match.fun(FUN)

        R<-list()
        pim<-parentIonMass(sequence)

        C_term <- 17.002740
        N_term <- 1.007825
        Oxygen <- 15.994915
        Carbon <- 12.000000
        Hydrogen <- 1.007825
        Nitrogen <- 14.003074
        Electron <- 0.000549

        for (i in 1:length(sequence)){
            input.sequence<-sequence[i]
            input.n<-nchar(input.sequence)
            input.modified <- as.integer(strsplit(modified[i], '')[[1]])
            input.pim<-pim[i]+(sum(as.double(as.character(modification[input.modified+1]))))

            if (input.n != length(input.modified))
                stop (paste("unvalid argument",i,"- number of AA and modification config differ! stop."))

            out <- .C("computeFragmentIonsModification",
                n=as.integer(input.n),
                pepSeq=as.character(input.sequence),
                pim=as.double(input.pim),
                b=as.double(rep(0.0, input.n)),
                y=as.double(rep(0.0, input.n)),
                modified=input.modified,
                modification=as.double(as.character(modification)))

            R[[length(R)+1]] <- as.data.frame(FUN(out$b, out$y))
        }

    } else{
        FUN<-match.fun(FUN)

        R<-list()
        pim<-parentIonMass(sequence)

        C_term <- 17.002740
        N_term <- 1.007825
        Oxygen <- 15.994915
        Carbon <- 12.000000
        Hydrogen <- 1.007825
        Nitrogen <- 14.003074
        Electron <- 0.000549

        for (i in 1:length(sequence)){
            pepseq<-sequence[i]
            pepseq.pim<-pim[i]

            out <- .C("computeFragmentIons",
                n=as.integer(nchar(pepseq)),
                pepSeq=as.character(pepseq),
                pim=as.double(pepseq.pim),
                b=as.double(rep(0.0,nchar(pepseq))),
                y=as.double(rep(0.0,nchar(pepseq))))



            R[[length(R)+1]] <- as.data.frame(FUN(out$b, out$y))
        }
    }
    return(R)
}

