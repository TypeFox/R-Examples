#' Sequential 3D descriptor.
#'
#' \code{SequentialTd} returns the concatenation of 3D descriptors of amino
#' acids in a protein sequence.
#'
#' @param x A string of amino acid letters
#' @return A 3*n dimensional numeric vector where n is the protein length
#'
#' @export SequentialTd
#'
#' @examples
#' x = "LALHLLLLHMHMMDRSLLLH"
#' SequentialTd(x)


SequentialTd<-function(x)
{

n=nchar(x)
x=unlist(strsplit(x, split = ""))

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


n=length(x)
result<-c()
for(i in 1:n)
{
result<-c(result,(td[rownames(td)==x[i]]))
}

return(result)
}

