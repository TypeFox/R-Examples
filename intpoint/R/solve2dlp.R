solve2dlp <-
function(t=1, c=NULL, bA = NULL, A = NULL, bm = NULL,
	m = NULL, bM = NULL, M = NULL, z=1, ip=1, e = 1e-04,
	a1 = 1, a2 = 0.97)
{
Rest<-buildSb(m=m,bm=bm,M=M, bM=bM, A=A, bA=bA)

Rest.M<-Rest$S1
Rest.b<-Rest$b1

show2d(Rest.M, Rest.b, z=z, ip=ip, c=c, t=t)
}
