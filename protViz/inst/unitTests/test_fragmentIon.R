#R

test_fragmentIons <-
function(){

    # check 0
    m<-fragmentIon('HTLNQIDSVK')[[1]]

    b_assumed<-c(138.0662, 239.1139, 352.1979, 466.2409, 594.2994, 
        707.3835, 822.4104, 909.4425, 1008.5109, 1136.6058)

    lapply(1:length(m$b), function(i){
        checkEqualsNumeric(m$b[i], b_assumed[i], tolerance=1.0e-3)
    })


    checkEqualsNumeric(m$y[5], 561.3242, tolerance=1.0e-3)



    # check 1

     peptide.AA<-"KINHSFLR";
     peptide.AA.weights<-c(128.09496,113.08406,114.04293,
         137.05891,87.03203,147.06841,113.08406,156.10111);

     checkEqualsNumeric(fragmentIon(peptide.AA.weights)[[1]]$b, fragmentIon(peptide.AA)[[1]]$b, tolerance=1.0e-3)
     checkEqualsNumeric(fragmentIon(peptide.AA.weights)[[1]]$y, fragmentIon(peptide.AA)[[1]]$y, tolerance=1.0e-3)

     
# check 2
     peptide<-"XGXFNAGVGK"
     fi<-fragmentIon(c("XGXFNAGVGK"))[[1]]
     pim<-parentIonMass(peptide)

     checkEqualsNumeric(fi$b[2:8]+fi$y[8:2], rep(pim+1.008, 7), tolerance=1.0e-3)


}

test_fragmentIons()
