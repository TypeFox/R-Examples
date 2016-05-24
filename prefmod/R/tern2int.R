# tern vector of 0/1/2
# int integer
tern2int<-function(tern)
{
  n<-length(tern);
  dec<-0
  for (k in 0:(n-1))
    dec<-dec+tern[n-k]*3^k
  dec

}
