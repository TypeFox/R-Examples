#' Discrete blosum descriptor.
#'
#' \code{DiscreteBlosum} returns the sum of blosum descriptors of amino acids in
#' a protein sequence.
#'
#' @param x A string of amino acid letters
#' @return A 20 dimensional numeric vector
#'
#' @export DiscreteBlosum
#'
#' @examples
#' x = "LALHLLLLHMHMMDRSLLLH"
#' DiscreteBlosum(x)

DiscreteBlosum<-function (x)
{

AAs = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
n=nchar(x)
x=strsplit(x, split = "")
AAC=summary(factor(x[[1]],levels=AAs),maxsum=21)

    blosum = t(data.frame(
	A = c(4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0),
	R = c(-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3),
	N = c(-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3),
	D = c(-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3),
	C = c(0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1),
	E = c(-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2),
	Q = c(-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2),
	G = c(0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3),
	H = c(-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3),
	I = c(-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3),
	L = c(-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1),
	K = c(-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2),
	M = c(-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1),
	F = c(-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1),
      P = c(-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2),
	S = c(1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2),
	T = c(0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0),
	W = c(-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3),
	Y = c(-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1),
	V = c(0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4)))

AAC1=AAC[1]*blosum[1,]
AAC2=AAC[2]*blosum[2,]
AAC3=AAC[3]*blosum[3,]
AAC4=AAC[4]*blosum[4,]
AAC5=AAC[5]*blosum[5,]
AAC6=AAC[6]*blosum[6,]
AAC7=AAC[7]*blosum[7,]
AAC8=AAC[8]*blosum[8,]
AAC9=AAC[9]*blosum[9,]
AAC10=AAC[10]*blosum[10,]
AAC11=AAC[11]*blosum[11,]
AAC12=AAC[12]*blosum[12,]
AAC13=AAC[13]*blosum[13,]
AAC14=AAC[14]*blosum[14,]
AAC15=AAC[15]*blosum[15,]
AAC16=AAC[16]*blosum[16,]
AAC17=AAC[17]*blosum[17,]
AAC18=AAC[18]*blosum[18,]
AAC19=AAC[19]*blosum[19,]
AAC20=AAC[20]*blosum[20,]

result=(AAC1+AAC2+AAC3+AAC4+AAC5+AAC6+AAC7+AAC8+AAC9+AAC10+
AAC11+AAC12+AAC13+AAC14+AAC15+AAC16+AAC17+AAC18+AAC19+AAC20)/(6*n)

return(result)
}
