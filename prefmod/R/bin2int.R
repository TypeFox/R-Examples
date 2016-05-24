# bin vector of 0/1
# int integer
bin2int<-function(bin)
{
  n<-length(bin);
  int<-0
  for (k in 0:(n-1))
    int<-int+bin[n-k]*2^k
  int
}
