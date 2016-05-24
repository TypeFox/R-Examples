seebank<-function(QB)
{
  ########    show each question in question bank
  
  for(i in 1:length(QB))
    {
      cat(sep="\n", paste(sep=" ", "############", i))
      print(QB[[i]])
      readline()
    }

}

