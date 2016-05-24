hasRDRAND<-function()
{
   .Call("Rrdrand_has_rdrand", PACKAGE="Rrdrand")
}
