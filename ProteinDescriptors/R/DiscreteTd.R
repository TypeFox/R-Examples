#' Discrete 3D descriptor.
#'
#' \code{DiscreteTd} returns the sum of 3D descriptors of amino acids in
#' a protein sequence.
#'
#' @param x A string of amino acid letters
#' @return A 3 dimensional numeric vector
#'
#' @export DiscreteTd
#'
#' @examples
#' x = "LALHLLLLHMHMMDRSLLLH"
#' DiscreteTd(x)

DiscreteTd<-function (x)
{
AAs = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
n=nchar(x)
x=strsplit(x, split = "")
AAC=summary(factor(x[[1]],levels=AAs),maxsum=21)

    td = t(data.frame(
	A = c(0.189,-3.989,1.989),
	R = c(5.007,0.834,-2.709),
	N = c(7.616,0.943,0.101),
	D = c(7.781,0.030,1.821),
	C = c(-5.929,-4.837,6.206),
	E = c(7.444,1.005,-2.121),
	Q = c(5.480,1.293,-3.091),
	G = c(4.096,0.772,7.120),
	H = c(3.488,6.754,-2.703),
	I = c(-7.883,-4.900,-2.230),
	L = c(-7.582,-3.724,-2.740),
	K = c(5.665,-0.166,-2.643),
	M = c(-5.200,-2.547,-3.561),
	F = c(-8.681,4.397,-0.732),
    P = c(4.281,-2.932,2.319),
	S = c(4.201,-1.948,1.453),
	T = c(0.774,-3.192,0.666),
	W = c(-8.492,9.958,4.874),
	Y = c(-6.147,7.590,-2.065),
	V = c(-6.108,-5.341,-1.953)))


AAC1=AAC[1]*td[1,]
AAC2=AAC[2]*td[2,]
AAC3=AAC[3]*td[3,]
AAC4=AAC[4]*td[4,]
AAC5=AAC[5]*td[5,]
AAC6=AAC[6]*td[6,]
AAC7=AAC[7]*td[7,]
AAC8=AAC[8]*td[8,]
AAC9=AAC[9]*td[9,]
AAC10=AAC[10]*td[10,]
AAC11=AAC[11]*td[11,]
AAC12=AAC[12]*td[12,]
AAC13=AAC[13]*td[13,]
AAC14=AAC[14]*td[14,]
AAC15=AAC[15]*td[15,]
AAC16=AAC[16]*td[16,]
AAC17=AAC[17]*td[17,]
AAC18=AAC[18]*td[18,]
AAC19=AAC[19]*td[19,]
AAC20=AAC[20]*td[20,]


result=(AAC1+AAC2+AAC3+AAC4+AAC5+AAC6+AAC7+AAC8+AAC9+AAC10+
AAC11+AAC12+AAC13+AAC14+AAC15+AAC16+AAC17+AAC18+AAC19+AAC20)/n

return(result)
}
