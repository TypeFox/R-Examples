#' Discrete sequential sparse descriptor with split number=5.
#'
#' \code{DiscreteSequentialSparseFiveParts} returns the concatenation of the sum
#' of sparse descriptors of amino acids in each split of a protein sequence.
#'
#' @param x A string of amino acid letters
#' @return A 100 dimensional numeric vector
#'
#' @export DiscreteSequentialSparseFiveParts
#'
#' @examples
#' x = "LALHLLLLHMHMMDRSLLLH"
#' DiscreteSequentialSparseFiveParts(x)

DiscreteSequentialSparseFiveParts<-function (x)
{
AAs = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
n=nchar(x)
of=round(n/5)
x1=substring(x,1,of)
x2=substring(x,of+1,2*of)
x3=substring(x,2*of+1,3*of)
x4=substring(x,3*of+1,4*of)
x5=substring(x,4*of+1,n)

x1=strsplit(x1, split = "")
x2=strsplit(x2, split = "")
x3=strsplit(x3, split = "")
x4=strsplit(x4, split = "")
x5=strsplit(x5, split = "")

AACone=summary(factor(x1[[1]],levels=AAs),maxsum=21)
AACtwo=summary(factor(x2[[1]],levels=AAs),maxsum=21)
AACthree=summary(factor(x3[[1]],levels=AAs),maxsum=21)
AACfour=summary(factor(x4[[1]],levels=AAs),maxsum=21)
AACfive=summary(factor(x5[[1]],levels=AAs),maxsum=21)

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

AACthree1=AACthree[1]*sparse[1,]
AACthree2=AACthree[2]*sparse[2,]
AACthree3=AACthree[3]*sparse[3,]
AACthree4=AACthree[4]*sparse[4,]
AACthree5=AACthree[5]*sparse[5,]
AACthree6=AACthree[6]*sparse[6,]
AACthree7=AACthree[7]*sparse[7,]
AACthree8=AACthree[8]*sparse[8,]
AACthree9=AACthree[9]*sparse[9,]
AACthree10=AACthree[10]*sparse[10,]
AACthree11=AACthree[11]*sparse[11,]
AACthree12=AACthree[12]*sparse[12,]
AACthree13=AACthree[13]*sparse[13,]
AACthree14=AACthree[14]*sparse[14,]
AACthree15=AACthree[15]*sparse[15,]
AACthree16=AACthree[16]*sparse[16,]
AACthree17=AACthree[17]*sparse[17,]
AACthree18=AACthree[18]*sparse[18,]
AACthree19=AACthree[19]*sparse[19,]
AACthree20=AACthree[20]*sparse[20,]


AACfour1=AACfour[1]*sparse[1,]
AACfour2=AACfour[2]*sparse[2,]
AACfour3=AACfour[3]*sparse[3,]
AACfour4=AACfour[4]*sparse[4,]
AACfour5=AACfour[5]*sparse[5,]
AACfour6=AACfour[6]*sparse[6,]
AACfour7=AACfour[7]*sparse[7,]
AACfour8=AACfour[8]*sparse[8,]
AACfour9=AACfour[9]*sparse[9,]
AACfour10=AACfour[10]*sparse[10,]
AACfour11=AACfour[11]*sparse[11,]
AACfour12=AACfour[12]*sparse[12,]
AACfour13=AACfour[13]*sparse[13,]
AACfour14=AACfour[14]*sparse[14,]
AACfour15=AACfour[15]*sparse[15,]
AACfour16=AACfour[16]*sparse[16,]
AACfour17=AACfour[17]*sparse[17,]
AACfour18=AACfour[18]*sparse[18,]
AACfour19=AACfour[19]*sparse[19,]
AACfour20=AACfour[20]*sparse[20,]

AACfive1=AACfive[1]*sparse[1,]
AACfive2=AACfive[2]*sparse[2,]
AACfive3=AACfive[3]*sparse[3,]
AACfive4=AACfive[4]*sparse[4,]
AACfive5=AACfive[5]*sparse[5,]
AACfive6=AACfive[6]*sparse[6,]
AACfive7=AACfive[7]*sparse[7,]
AACfive8=AACfive[8]*sparse[8,]
AACfive9=AACfive[9]*sparse[9,]
AACfive10=AACfive[10]*sparse[10,]
AACfive11=AACfive[11]*sparse[11,]
AACfive12=AACfive[12]*sparse[12,]
AACfive13=AACfive[13]*sparse[13,]
AACfive14=AACfive[14]*sparse[14,]
AACfive15=AACfive[15]*sparse[15,]
AACfive16=AACfive[16]*sparse[16,]
AACfive17=AACfive[17]*sparse[17,]
AACfive18=AACfive[18]*sparse[18,]
AACfive19=AACfive[19]*sparse[19,]
AACfive20=AACfive[20]*sparse[20,]

result1=(AACone1+AACone2+AACone3+AACone4+AACone5+AACone6+AACone7+AACone8+AACone9+AACone10+
AACone11+AACone12+AACone13+AACone14+AACone15+AACone16+AACone17+AACone18+AACone19+AACone20)

result2=(AACtwo1+AACtwo2+AACtwo3+AACtwo4+AACtwo5+AACtwo6+AACtwo7+AACtwo8+AACtwo9+AACtwo10+
AACtwo11+AACtwo12+AACtwo13+AACtwo14+AACtwo15+AACtwo16+AACtwo17+AACtwo18+AACtwo19+AACtwo20)

result3=(AACthree1+AACthree2+AACthree3+AACthree4+AACthree5+AACthree6+AACthree7+AACthree8+AACthree9+AACthree10+
AACthree11+AACthree12+AACthree13+AACthree14+AACthree15+AACthree16+AACthree17+AACthree18+AACthree19+AACthree20)

result4=(AACfour1+AACfour2+AACfour3+AACfour4+AACfour5+AACfour6+AACfour7+AACfour8+AACfour9+AACfour10+
AACfour11+AACfour12+AACfour13+AACfour14+AACfour15+AACfour16+AACfour17+AACfour18+AACfour19+AACfour20)

result5=(AACfive1+AACfive2+AACfive3+AACfive4+AACfive5+AACfive6+AACfive7+AACfive8+AACfive9+AACfive10+
AACfive11+AACfive12+AACfive13+AACfive14+AACfive15+AACfive16+AACfive17+AACfive18+AACfive19+AACfive20)

result<-c(result1,result2,result3,result4,result5)/n
return(result)
}

