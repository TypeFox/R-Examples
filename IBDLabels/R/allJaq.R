allJaq <-
    function(){

        out <- cbind(label = 0:15,
                 jaq15 = c(1,3,NA,4,2,5,6,9,11,10,7,13,12,14,8,15),
                 jaq9 = c(1,3,NA,3,2,4,5,7,8,7,5,8,8,8,6,9))

        return( out )

    }
