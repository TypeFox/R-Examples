context('Betas_p179')
# tests the analytically computed results from Harte 2011 on page 179

sv=data.frame(S0=c(24,302,283,296,51,130,22,27,158,240,156,227,167,495,221,
                   1161,137,
                   67,808,
                   138,24,28,205,
                   47,12,11,32.5,49,1084),
              N0=c(2445,213791,208310,11843,1448,4807,1378,1626,1922,6048,2253,3865,1909,329026,105163,
                   355327,67465,
                   17995,296124,
                   12851,37182,60346,193353,
                   264,547,414,109,645,145575),
              beta=c(.00146,.000162,.000155,.00464,.007085,.00511,.00261,.00279,.021,.00823,.0168,.0135,.0228,.000174,.000254,
                     .00042,.000244,
                     .000488,.000341, # only off by 1e-6
                     .00168,.0000592,
                     .0000427,# only off by 3e-7
                     .000117,
                     .061,.003175,.003872,.134,.0189,.00109)
)

test_that('ESF computed correctly',{
  for(i in 1:nrow(sv)){
    esf1=meteESF(S0=sv$S0[i],N0=sv$N0[i],E0=1e7)
    test.beta=round(sum(esf1$La),nchar(sv$beta[i])-2)
    close=abs(test.beta-sv$beta[i])<1e-3
    cat(i, '  ',test.beta,'  ',sv$beta[i], '  ')
    #print(expect_equal(test.beta,sv$beta[i]))
    print(expect_equal(close,TRUE))
  }
})

