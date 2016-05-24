pos.charstr <-
function(string,charstr)
{  # Is charstr in string?
	z2 <- regexpr(charstr,string)
	if (z2[1]==-1) {z1 = NA} else
	{	zz <- stringf(as.character(string))
		char <- substr(charstr,1,1)
        if (char%in%zz)
           { z <- grep(char,zz)}  else
            { z <- NA }
        z1<- z[1] # returns first occurrence only
    }
  return (z1)
}
