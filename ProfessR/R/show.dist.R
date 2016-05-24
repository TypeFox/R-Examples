`show.dist` <-
function(W)
  {

  for(i in 1:length(W$LETS))
    {
      print(paste(sep=' ', i, W$LETS[i], length(W$lett[W$lett==W$LETS[i]])))
    }

print("")
  
  for(i in 1:length(W$SCRS))
    {


      print(paste(sep=' ', i, W$SCRS[i], length(W$scor[W$scor==W$SCRS[i]])))
    }
  print("")
  
  print(paste(sep=' ', "mean=",mean(W$scor)))


  }

