# 2-4-2010 MRC-Epid JHZ

SH2B1 <- read.table("SH2B1.txt")
n9 <- matrix(0,3,3)
rs8055982 <- as.numeric(SH2B1[1,-c(1:5,10662)])
rs7498665 <- as.numeric(SH2B1[1,-c(1:5,10662)])
n <- 3552
for (i in 1:n)
{
   i3 <- 3*i
   g1 <- rs8055982[c(i3-2,i3-1,i3)]
   g2 <- rs7498665[c(i3-2,i3-1,i3)]
   n9 <- n9 + g1%o%g2
}
h12 <- klem(n9)
ld <- LD22(h12$h,2*n)
ld
