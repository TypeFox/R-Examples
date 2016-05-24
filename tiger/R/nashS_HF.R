`nashS_HF` <-
function (modelled, measured, weigth=NA)

{

    #Returns the NS with high weight on high flows and low weight on
    #low flows
        #empirical cumul. density function
        t.ecdf <- ecdf(measured)
        #Probability to be below value
        t.weigth <- t.ecdf(measured)
        t.weigth <- t.weigth*weigth
        return(nashS(measured,modelled,weigth=t.weigth))

}

