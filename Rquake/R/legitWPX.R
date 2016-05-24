legitWPX<-function(twpx, quiet=TRUE)
{
###  test to see if there are real picks in the WPX list
 ###  print("test legit")
 ###  print(data.frame(twpx))

  if(is.null(twpx))
    {
        return(0)
    }

  tval  = all(is.na(twpx$name))


  if(tval)
    {
      if(!quiet) cat("ERROR WPX (no picks)", sep="\n")
      return(0)

    }
  else
    {

      return(1)

    }


}
