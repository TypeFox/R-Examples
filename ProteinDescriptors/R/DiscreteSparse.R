#' Discrete sparse descriptor.
#'
#' \code{DiscreteSparse} returns the sum of sparse descriptors of amino acids in
#' a protein sequence.
#'
#' @param x A string of amino acid letters
#' @return A 20 dimensional numeric vector
#'
#' @export DiscreteSparse
#'
#' @examples
#' x = "LALHLLLLHMHMMDRSLLLH"
#' DiscreteSparse(x)

DiscreteSparse<-function (x)
{

AAs = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

n=nchar(x)
x=strsplit(x, split = "")
AAC=summary(factor(x[[1]],levels=AAs),maxsum=21)

    sparse = t(data.frame(
	A = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	R = c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	N = c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	D = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	C = c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	E = c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	Q = c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
    	G = c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
	H = c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
	I = c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
	L = c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
	K = c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
	M = c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
	F = c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
      P = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
	S = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
	T = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
	W = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
	Y = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),
	V = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)))

AAC1=AAC[1]*sparse[1,]
AAC2=AAC[2]*sparse[2,]
AAC3=AAC[3]*sparse[3,]
AAC4=AAC[4]*sparse[4,]
AAC5=AAC[5]*sparse[5,]
AAC6=AAC[6]*sparse[6,]
AAC7=AAC[7]*sparse[7,]
AAC8=AAC[8]*sparse[8,]
AAC9=AAC[9]*sparse[9,]
AAC10=AAC[10]*sparse[10,]
AAC11=AAC[11]*sparse[11,]
AAC12=AAC[12]*sparse[12,]
AAC13=AAC[13]*sparse[13,]
AAC14=AAC[14]*sparse[14,]
AAC15=AAC[15]*sparse[15,]
AAC16=AAC[16]*sparse[16,]
AAC17=AAC[17]*sparse[17,]
AAC18=AAC[18]*sparse[18,]
AAC19=AAC[19]*sparse[19,]
AAC20=AAC[20]*sparse[20,]


result=(AAC1+AAC2+AAC3+AAC4+AAC5+AAC6+AAC7+AAC8+AAC9+AAC10+
        AAC11+AAC12+AAC13+AAC14+AAC15+AAC16+AAC17+AAC18+AAC19+AAC20)/n

return(result)
}
