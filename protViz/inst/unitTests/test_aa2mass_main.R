#R


test_aa2mass_main <-
function(){

    checkEqualsNumeric(aa2mass("MELIVSK")[[1]], 
    c(131.04049, 129.04259, 113.08406, 113.08406, 99.06841, 87.03203, 128.09496),
    tolerance=0.001)

}

test_aa2mass_main()
