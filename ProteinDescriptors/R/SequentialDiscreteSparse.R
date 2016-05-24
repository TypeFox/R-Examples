#' Sequential discrete sparse descriptor with step size=2.
#'
#' \code{SequentialDiscreteSparse} returns the sum of the concatenation of
#' sparse descriptors of amino acids at every step size in a protein sequence.
#'
#' @param x A string of amino acid letters
#' @return A 40 dimensional numeric vector
#'
#' @export SequentialDiscreteSparse
#'
#' @examples
#' x = "LALHLLLLHMHMMDRSLLLH"
#' SequentialDiscreteSparse(x)

SequentialDiscreteSparse<-function (x)
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

AACone1=AACone[1]*sparse[1,]
AACone2=AACone[2]*sparse[2,]
AACone3=AACone[3]*sparse[3,]
AACone4=AACone[4]*sparse[4,]
AACone5=AACone[5]*sparse[5,]
AACone6=AACone[6]*sparse[6,]
AACone7=AACone[7]*sparse[7,]
AACone8=AACone[8]*sparse[8,]
AACone9=AACone[9]*sparse[9,]
AACone10=AACone[10]*sparse[10,]
AACone11=AACone[11]*sparse[11,]
AACone12=AACone[12]*sparse[12,]
AACone13=AACone[13]*sparse[13,]
AACone14=AACone[14]*sparse[14,]
AACone15=AACone[15]*sparse[15,]
AACone16=AACone[16]*sparse[16,]
AACone17=AACone[17]*sparse[17,]
AACone18=AACone[18]*sparse[18,]
AACone19=AACone[19]*sparse[19,]
AACone20=AACone[20]*sparse[20,]

AACtwo1=AACtwo[1]*sparse[1,]
AACtwo2=AACtwo[2]*sparse[2,]
AACtwo3=AACtwo[3]*sparse[3,]
AACtwo4=AACtwo[4]*sparse[4,]
AACtwo5=AACtwo[5]*sparse[5,]
AACtwo6=AACtwo[6]*sparse[6,]
AACtwo7=AACtwo[7]*sparse[7,]
AACtwo8=AACtwo[8]*sparse[8,]
AACtwo9=AACtwo[9]*sparse[9,]
AACtwo10=AACtwo[10]*sparse[10,]
AACtwo11=AACtwo[11]*sparse[11,]
AACtwo12=AACtwo[12]*sparse[12,]
AACtwo13=AACtwo[13]*sparse[13,]
AACtwo14=AACtwo[14]*sparse[14,]
AACtwo15=AACtwo[15]*sparse[15,]
AACtwo16=AACtwo[16]*sparse[16,]
AACtwo17=AACtwo[17]*sparse[17,]
AACtwo18=AACtwo[18]*sparse[18,]
AACtwo19=AACtwo[19]*sparse[19,]
AACtwo20=AACtwo[20]*sparse[20,]


result1=(AACone1+AACone2+AACone3+AACone4+AACone5+AACone6+AACone7+AACone8+AACone9+AACone10+
AACone11+AACone12+AACone13+AACone14+AACone15+AACone16+AACone17+AACone18+AACone19+AACone20)

result2=(AACtwo1+AACtwo2+AACtwo3+AACtwo4+AACtwo5+AACtwo6+AACtwo7+AACtwo8+AACtwo9+AACtwo10+
AACtwo11+AACtwo12+AACtwo13+AACtwo14+AACtwo15+AACtwo16+AACtwo17+AACtwo18+AACtwo19+AACtwo20)

result=c(result1,result2)/(n/2)
return(result)

}


