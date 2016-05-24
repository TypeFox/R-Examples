"compana" <- function(used, avail, test = c("randomisation", "parametric"),
                      rnv = 0.01, nrep = 500, alpha=0.1)
  {
    ### 1. Verifications
    test<-match.arg(test)
    used<-as.matrix(used)
    avail<-as.matrix(avail)
    if (any(apply(avail,2, function(x) {
        sum(x > .Machine$double.eps)
    })<2))
        stop("Availability different from zero\n for less than two animals for some habitat types. \n At least 2 animals are required for this analysis")
    if ((any(avail==0))&(test=="parametric")) {
      warning("parametric tests not suitable with 0 in avail; test has been set to \"randomisation\"")
      test<-"randomisation"
    }
    if (ncol(used)!=ncol(avail))
      stop("the two matrices should have the same dimensions")
    if (nrow(used)!=nrow(avail))
      stop("the two matrices should have the same dimensions")
    if (!all(colnames(used)==colnames(avail)))
      stop("the two matrices should have the same habitat names")
    if (is.null(colnames(used)))
      colnames(used) <- paste("Habitat", 1:ncol(used), sep = "")
    if (is.null(colnames(avail)))
        colnames(avail) <- paste("Habitat", 1:ncol(avail), sep = "")

    ## 2. Bases
    nh<-ncol(used)
    na<-nrow(used)
    proj1<-matrix(1, nrow=nrow(used), ncol=nrow(used))*(1/nrow(used))
    proj2<-matrix(0, nrow=nrow(used), ncol=nrow(used))
    if (test=="parametric")
      nrep=1
    sorties<-list()

    ## 3. First part: global test
    ##    relies on a call to the C function "aclamda"
    toto<-.C("aclambda", as.double(t(used)), as.double(t(avail)),
             as.integer(na), as.integer(nh),
             as.double(proj1), as.double(proj2), as.double(rnv),
             double(nrep), as.integer(nrep), double(nh), double(nh),
             PACKAGE="adehabitatHS")

    ## output
    vrand<-toto[[8]]
    sorties$used<-used
    sorties$avail<-avail
    sorties$type.test<-test
    if (test=="randomisation") {
        sorties$random.res<-list(sim=vrand, obs=vrand[1])
        sorties$test<-c(vrand[1], length(vrand[vrand<=vrand[1]])/nrep)
        names(sorties$test)<-c("Lambda", "P")
    }
    else {
        sorties$test<-c(vrand[1], ncol(used)-1,
                        1-pchisq(-na*log(vrand[1]), ncol(used)-1))
        names(sorties$test)<-c("Lambda", "df", "P")
    }


    ## 4. Second part: ranking matrix for habitat types
    if (test=="randomisation") {

        ## 4.1 For randomization tests:
        ##     relies on a call to the C function "rankma"
        toto<-.C("rankma", as.double(t(used)), as.double(t(avail)),
                 double(nh**2), double(nh**2), double(nh**2),
                 double(nh**2), as.integer(nh), as.integer(na),
                 as.integer(nrep), as.double(rnv), PACKAGE="adehabitatHS")

        rmp<-t(matrix(toto[[3]]/nrep, nh, nh))
        rmm<-t(matrix(toto[[4]]/nrep, nh, nh))
        rmv<-t(matrix(toto[[5]], nh, nh))
        rmnb<-t(matrix(toto[[6]], nh, nh))
    }
    else {
        ## 4.2 For parametric tests:

        ## 4.2.1 output matrices
        used[used==0]<-rnv
        rmv<-matrix(0, nh, nh)
        rmse<-matrix(0, nh, nh)
        rmm<-matrix(0, nh, nh)
        rmp<-matrix(0, nh, nh)
        rmnb<-matrix(0, nh, nh)


        for (i in 1:nh) {
            for (j in 1:nh) {

                ## The matrix of the difference of log ratios
                dlr<-log(used[,i]/used[,j])-log(avail[,i]/avail[,j])

                ## mean DLR
                rmv[i,j]<-mean(dlr)

                ## standard deviations
                rmse[i,j]<-sqrt(var(dlr)/na)

                ## t-test statistic
                if (i!=j)
                    rmv[i,j]<-rmv[i,j]/rmse[i,j]

                ## P-values
                rmp[i,j]<-pt(rmv[i,j], na-1)
                rmm[i,j]<-1-rmp[i,j]
                rmnb[i,j]<-na
            }
        }
    }

    ## preparation of the ranking matrix
    rm<-matrix("0", nh, nh)


    ## ranking matrix: keeps only the signs of the mean DLR
    for (i in 1:nh) {
        for (j in 1:nh) {
            if (rmv[i,j]<0)
                rm[i,j]<-"-"
            if (rmv[i,j]>0)
                rm[i,j]<-"+"
        }
    }

    ## adds the significance to the tests
    for (i in 1:nh) {
        for (j in 1:nh) {
            if (rmp[i,j] < (alpha/2)) {
                rm[i,j]<-"---"
            }
            if (rmm[i,j] < (alpha/2)) {
                rm[i,j]<-"+++"
            }
            if (i==j)
                rm[i,j]<-"0"
        }
    }


    ## Computes the ranks of habitat types (number of DLR >0)
    rank<-rep(0, nh)
    for (j in 1:nh) {
        for (i in 1:nh) {
            if (rmv[j,i]>0)
                rank[j]<-rank[j]+1
        }
    }

    ## row and column names for the various matrices
    names(rank)<-colnames(avail)
    rownames(rm)<-colnames(avail)
    colnames(rm)<-colnames(avail)
    rownames(rmv)<-colnames(avail)
    colnames(rmv)<-colnames(avail)
    rownames(rmp)<-colnames(avail)
    colnames(rmp)<-colnames(avail)
    rownames(rmm)<-colnames(avail)
    colnames(rmm)<-colnames(avail)
    rownames(rmnb)<-colnames(avail)
    colnames(rmnb)<-colnames(avail)

    ## output
    sorties$rmnb<-rmnb
    sorties$rank<-rank
    sorties$rm<-rm
    sorties$rmv<-rmv

    sorties$profile<-profilehab(rm, rank)
    class(sorties)<-"compana"
    return(sorties)
  }

