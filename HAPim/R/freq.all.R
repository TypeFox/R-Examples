`freq.all` <-
function(hap){

 nb.marq 	= 	length(hap[1,])
 freq.a 	= 	list()

    for(i in 1:nb.marq)

    {       
        sy          = tabulate(na.omit(hap[,i]))     
        toty        = sum(sy)
	freqy       = sy/toty
        freq.a[[i]] = freqy

    }


 freq.a


 }

