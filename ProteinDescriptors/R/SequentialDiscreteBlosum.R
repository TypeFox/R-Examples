#' Sequential discrete blosum descriptor with step size=2.
#'
#' \code{SequentialDiscreteBlosum} returns the sum of the concatenation of
#' blosum descriptors of amino acids at every step size in a protein sequence.
#'
#' @param x A string of amino acid letters
#' @return A 40 dimensional numeric vector
#'
#' @export SequentialDiscreteBlosum
#'
#' @examples
#' x = "LALHLLLLHMHMMDRSLLLH"
#' SequentialDiscreteBlosum(x)

SequentialDiscreteBlosum<-function (x)
{
AAs = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
n=nchar(x)
x1=""
x2=""
for(i in seq(1,n,2))
{
x1<-paste(x1,substr(x,i,i))
x2<-paste(x2,substr(x,i+1,i+1))
}
x1<-gsub(" ", "", x1)
x2<-gsub(" ", "", x2)


n1=nchar(x1)
x1=strsplit(x1, split = "")
AACone=summary(factor(x1[[1]],levels=AAs),maxsum=21)

n2=nchar(x2)
x2=strsplit(x2, split = "")
AACtwo=summary(factor(x2[[1]],levels=AAs),maxsum=21)

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

AACone1=AACone[1]*blosum[1,]
AACone2=AACone[2]*blosum[2,]
AACone3=AACone[3]*blosum[3,]
AACone4=AACone[4]*blosum[4,]
AACone5=AACone[5]*blosum[5,]
AACone6=AACone[6]*blosum[6,]
AACone7=AACone[7]*blosum[7,]
AACone8=AACone[8]*blosum[8,]
AACone9=AACone[9]*blosum[9,]
AACone10=AACone[10]*blosum[10,]
AACone11=AACone[11]*blosum[11,]
AACone12=AACone[12]*blosum[12,]
AACone13=AACone[13]*blosum[13,]
AACone14=AACone[14]*blosum[14,]
AACone15=AACone[15]*blosum[15,]
AACone16=AACone[16]*blosum[16,]
AACone17=AACone[17]*blosum[17,]
AACone18=AACone[18]*blosum[18,]
AACone19=AACone[19]*blosum[19,]
AACone20=AACone[20]*blosum[20,]

AACtwo1=AACtwo[1]*blosum[1,]
AACtwo2=AACtwo[2]*blosum[2,]
AACtwo3=AACtwo[3]*blosum[3,]
AACtwo4=AACtwo[4]*blosum[4,]
AACtwo5=AACtwo[5]*blosum[5,]
AACtwo6=AACtwo[6]*blosum[6,]
AACtwo7=AACtwo[7]*blosum[7,]
AACtwo8=AACtwo[8]*blosum[8,]
AACtwo9=AACtwo[9]*blosum[9,]
AACtwo10=AACtwo[10]*blosum[10,]
AACtwo11=AACtwo[11]*blosum[11,]
AACtwo12=AACtwo[12]*blosum[12,]
AACtwo13=AACtwo[13]*blosum[13,]
AACtwo14=AACtwo[14]*blosum[14,]
AACtwo15=AACtwo[15]*blosum[15,]
AACtwo16=AACtwo[16]*blosum[16,]
AACtwo17=AACtwo[17]*blosum[17,]
AACtwo18=AACtwo[18]*blosum[18,]
AACtwo19=AACtwo[19]*blosum[19,]
AACtwo20=AACtwo[20]*blosum[20,]


result1<-(AACone1+AACone2+AACone3+AACone4+AACone5+AACone6+AACone7+AACone8+AACone9+AACone10+
AACone11+AACone12+AACone13+AACone14+AACone15+AACone16+AACone17+AACone18+AACone19+AACone20)

result2<-(AACtwo1+AACtwo2+AACtwo3+AACtwo4+AACtwo5+AACtwo6+AACtwo7+AACtwo8+AACtwo9+AACtwo10+
AACtwo11+AACtwo12+AACtwo13+AACtwo14+AACtwo15+AACtwo16+AACtwo17+AACtwo18+AACtwo19+AACtwo20)

result<-c(result1,result2)/(3*n)

return(result)

}
