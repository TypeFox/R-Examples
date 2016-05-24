###landscape.write.foreign
landscape.write.foreign <- function(Rland, fn = "foreign.genepop", fmt="genepop",...)
{
  if  (is.landscape(Rland))
    {
      if (fmt %in% c("GenePop","Genepop","genepop"))
        {
          landscape.write.genepop(Rland,fn,...)
        }
    }
  else
    {
      print ("rland not in landscape format")
    }
}
