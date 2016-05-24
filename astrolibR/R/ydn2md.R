	ydn2md = function(yr,dy) {
 
	# days before start of each month.
	ydays = c(0,31,59,90,120,151,181,212,243,273,304,334,366) + 1
 
	leap =  (((yr %% 4)==0) & ((yr %% 100)!=0)) | ((yr %% 400)==0) 
                    
        if (leap) ydays[-(1:2)] = ydays[-(1:2)] + 1

        m = findInterval(dy, ydays)

		d = dy - ydays[m] + 1

        if(any(m>12))
            stop('error in day number')

	return(list(m=m,d=d))
      }

