#R


test_aa2mass_main__verus__parentIonMass <-
function(){

    # just a test
    peptides<-c('HTLNQIDSVK', 'ALGGEDVR', 'TPIVGQPSIPGGPVR')
    C_term <- 17.002740
    N_term <- 1.007825
    H_ <- 1.008

    checkEqualsNumeric(parentIonMass(peptides),
        unlist(lapply(aa2mass(peptides), sum)) + C_term + N_term + H_, 
        tolerance=0.001)

}

test_aa2mass_main__verus__parentIonMass()
