log2int<-function(bin)
{
  n<-length(bin);
  int<-0
  for (k in 0:(n-1))
    int<-int+as.numeric(bin[n-k])*2^k
  int
}
