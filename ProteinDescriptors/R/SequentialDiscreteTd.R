#' Sequential discrete 3D descriptor with step size=2.
#'
#' \code{SequentialDiscreteTd} returns the sum of the concatenation of 3D
#' descriptors of amino acids at every step size in a protein sequence.
#'
#' @param x A string of amino acid letters
#' @return A 6 dimensional numeric vector
#'
#' @export SequentialDiscreteTd
#'
#' @examples
#' x = "LALHLLLLHMHMMDRSLLLH"
#' SequentialDiscreteTd(x)

SequentialDiscreteTd<-function (x)
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


AACone1=AACone[1]*td[1,]
AACone2=AACone[2]*td[2,]
AACone3=AACone[3]*td[3,]
AACone4=AACone[4]*td[4,]
AACone5=AACone[5]*td[5,]
AACone6=AACone[6]*td[6,]
AACone7=AACone[7]*td[7,]
AACone8=AACone[8]*td[8,]
AACone9=AACone[9]*td[9,]
AACone10=AACone[10]*td[10,]
AACone11=AACone[11]*td[11,]
AACone12=AACone[12]*td[12,]
AACone13=AACone[13]*td[13,]
AACone14=AACone[14]*td[14,]
AACone15=AACone[15]*td[15,]
AACone16=AACone[16]*td[16,]
AACone17=AACone[17]*td[17,]
AACone18=AACone[18]*td[18,]
AACone19=AACone[19]*td[19,]
AACone20=AACone[20]*td[20,]

AACtwo1=AACtwo[1]*td[1,]
AACtwo2=AACtwo[2]*td[2,]
AACtwo3=AACtwo[3]*td[3,]
AACtwo4=AACtwo[4]*td[4,]
AACtwo5=AACtwo[5]*td[5,]
AACtwo6=AACtwo[6]*td[6,]
AACtwo7=AACtwo[7]*td[7,]
AACtwo8=AACtwo[8]*td[8,]
AACtwo9=AACtwo[9]*td[9,]
AACtwo10=AACtwo[10]*td[10,]
AACtwo11=AACtwo[11]*td[11,]
AACtwo12=AACtwo[12]*td[12,]
AACtwo13=AACtwo[13]*td[13,]
AACtwo14=AACtwo[14]*td[14,]
AACtwo15=AACtwo[15]*td[15,]
AACtwo16=AACtwo[16]*td[16,]
AACtwo17=AACtwo[17]*td[17,]
AACtwo18=AACtwo[18]*td[18,]
AACtwo19=AACtwo[19]*td[19,]
AACtwo20=AACtwo[20]*td[20,]


result1=(AACone1+AACone2+AACone3+AACone4+AACone5+AACone6+AACone7+AACone8+AACone9+AACone10+
AACone11+AACone12+AACone13+AACone14+AACone15+AACone16+AACone17+AACone18+AACone19+AACone20)/n

result2=(AACtwo1+AACtwo2+AACtwo3+AACtwo4+AACtwo5+AACtwo6+AACtwo7+AACtwo8+AACtwo9+AACtwo10+
AACtwo11+AACtwo12+AACtwo13+AACtwo14+AACtwo15+AACtwo16+AACtwo17+AACtwo18+AACtwo19+AACtwo20)/n

result<-c(result1,result2)
return(result)

}
