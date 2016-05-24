string.blank.omit <-
function (string)  
#  Suppresses blanks in a string: returns string without blanks
{     k <- 1
      str <- ""
      for (j in 1:nchar(string))
      { if(substr(string,j,j)!=" ") 
          { str <- ifelse (str=="",substr(string,j,j),paste(str,substr(string,j,j),sep=""))
            k <- k + 1
          }
      }
  return (str)
}
