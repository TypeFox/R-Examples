wrapperTraj <-
function(Data, Time, ID = FALSE)
{
  s1 = step1measures(Data, Time, ID = ID)
  s2 = step2factors(s1)
  s3 = step3clusters(s2)
  
  return(s3)
}
