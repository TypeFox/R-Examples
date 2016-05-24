## ----OptionsAndLibraries, include=FALSE, message=FALSE------------------------------------------------------------------------------------
# knitr has to be loaded for 'set_parent' and CRAN checks and also for opt_chunk during build process.
library(knitr)
if (exists("opts_chunk")) {
  opts_chunk$set(concordance=TRUE)
  opts_chunk$set(tidy.opts=list(keep.blank.line=FALSE, width.cutoff=95))
  opts_chunk$set(size="footnotesize")
  opts_chunk$set(cache=TRUE)
  opts_chunk$set(autodep=TRUE)
}
library(gMCP, quietly=TRUE)
options(width=140)
options(digits=4)
gMCPReport <- function(...) {invisible(NULL)}
graphGUI <- function(...) {invisible(NULL)}
#options(prompt="> ", continue="+ ")


## ----echo=FALSE,results='asis'------------------------------------------------------------------------------------------------------------
graph <- BonferroniHolm(3)
cat(graph2latex(graph, scale=0.7, labelTikZ="near start,fill=blue!20", scaleText=FALSE))

## ----echo=FALSE,results='asis'------------------------------------------------------------------------------------------------------------
graph <- BonferroniHolm(3)
cat(graph2latex(graph, scale=0.7, nodeTikZ="minimum size=1.2cm", scaleText=FALSE))
cat("\\\\$\\downarrow$ reject $H_1$\\\\")
graph <- rejectNode(graph, "H1")
cat(graph2latex(graph, scale=0.7, nodeTikZ="minimum size=1.2cm", scaleText=FALSE))
cat("\\\\$\\downarrow$ reject $H_3$\\\\\\ \\\\")
graph <- rejectNode(graph, "H3")
cat(graph2latex(graph, scale=0.7, nodeTikZ="minimum size=1.2cm", scaleText=FALSE))

## ----echo=TRUE, eval=FALSE----------------------------------------------------------------------------------------------------------------
#  library(gMCP)
#  graphGUI()

## ----echo=TRUE----------------------------------------------------------------------------------------------------------------------------
library(gMCP)
graph <- BonferroniHolm(3)
gMCP(graph, pvalues=c(0.01,0.07,0.02), alpha=0.05)

## ----echo=FALSE,results='asis'------------------------------------------------------------------------------------------------------------
graph <- BretzEtAl2011()
cat(graph2latex(graph, scale=0.7, scaleText=FALSE))

## ----echo=TRUE, tidy=FALSE----------------------------------------------------------------------------------------------------------------
m <- rbind(H11=c(0,   0.5, 0,   0.5, 0,   0  ),
           H21=c(1/3, 0,   1/3, 0,   1/3, 0  ),
           H31=c(0,   0.5, 0,   0,   0,   0.5),
           H12=c(0,   1,   0,   0,   0,   0  ),
           H22=c(0.5, 0,   0.5, 0,   0,   0  ),
           H32=c(0,   1,   0,   0,   0,   0  ))

graph <- matrix2graph(m)
graph <- setWeights(graph, c(1/3, 1/3, 1/3, 0, 0, 0))

## ----echo=TRUE----------------------------------------------------------------------------------------------------------------------------
print(graph)

## ----echo=TRUE----------------------------------------------------------------------------------------------------------------------------
graph@nodeAttr$X <- c(H11=100, H21=300, H31=500, H12=100, H22=300, H32=500)
graph@nodeAttr$Y <- c(H11=100, H21=100, H31=100, H12=300, H22=300, H32=300)	

## ----echo=TRUE----------------------------------------------------------------------------------------------------------------------------
graph <- placeNodes(graph, nrow=2)

## ----echo=TRUE----------------------------------------------------------------------------------------------------------------------------
cat(graph2latex(graph))

## ----echo=TRUE----------------------------------------------------------------------------------------------------------------------------
edgeAttr(graph, "H11", "H21", "labelX") <- 200
edgeAttr(graph, "H11", "H21", "labelY") <- 80

## ----echo=TRUE----------------------------------------------------------------------------------------------------------------------------
graphGUI("graph")

## ----echo=TRUE,tidy=FALSE-----------------------------------------------------------------------------------------------------------------
graph <- BretzEtAl2011()

# We can reject a single node:
print(rejectNode(graph, "H11"))

# Or given a vector of pvalues let the function gMCP do all the work:  
pvalues <- c(0.1, 0.008, 0.005, 0.15, 0.04, 0.006)
result <- gMCP(graph, pvalues)
print(result)

## ----echo=FALSE,results='asis'------------------------------------------------------------------------------------------------------------
cat(graph2latex(result@graphs[[4]], scale=0.7, scaleText=FALSE))

## ----echo=TRUE, eval=FALSE----------------------------------------------------------------------------------------------------------------
#  gMCPReport(result, "Report.tex")

## ----include=FALSE------------------------------------------------------------------------------------------------------------------------
d1 <- c(1.68005156523844, 1.95566697423700, 0.00137860945822299, 0.660052238622464, 
		1.06731835721526, 0.39479303427265, -0.312462050794408, 0.323637755662837, 
		0.490976552328251, 2.34240774442652)

d2 <- c(0.507878380203451, 1.60461475524144, 2.66959621483759, 0.0358289240280020, 
		-1.13014087491324, 0.792461583741794, 0.0701657425268248, 3.15360436883856, 
		0.217669661552567, 1.23979492014026)

d3 <- c(-1.31499425534849, 1.62201370145649, 0.89391826766116, 0.845473572033649, 
		2.17912435223573, 1.07521368050267, 0.791598289847664, 1.58537210294519, 
		-0.079778759456515, 0.97295072606043)

est <- c(0.860382, 0.9161474, 0.9732953)
s <- c(0.8759528, 1.291310, 0.8570892)
pval <- c(0.01260, 0.05154, 0.02124)/2

df <- 9
# Statistics:
st <- qt(pval/2, df=df, lower.tail=FALSE)
# Estimates:
est <- st*s/sqrt(10)

## ----confint, echo=TRUE, tidy=FALSE-------------------------------------------------------------------------------------------------------
# Estimates:
est <- c("H1"=0.860382, "H2"=0.9161474, "H3"=0.9732953)
# Sample standard deviations:
ssd <- c("H1"=0.8759528, "H2"=1.291310, "H3"=0.8570892)

pval <- c(0.01260, 0.05154, 0.02124)/2

simConfint(BonferroniHolm(3), pvalues=pval, 
		confint=function(node, alpha) {
			c(est[node]-qt(1-alpha,df=9)*ssd[node]/sqrt(10), Inf)
		}, estimates=est, alpha=0.025, mu=0, alternative="greater")

# Note that the sample standard deviations in the following call
# will be calculated from the pvalues and estimates.
simConfint(BonferroniHolm(3), pvalues=pval, 
		confint="t", df=9, estimates=est, alpha=0.025, alternative="greater")

## ----echo=FALSE,results='asis'------------------------------------------------------------------------------------------------------------

cat(graph2latex(parallelGatekeeping(), nodeTikZ="minimum size=1.2cm", tikzEnv=FALSE))
cat(graph2latex(improvedParallelGatekeeping(), nodeTikZ="minimum size=1.2cm", tikzEnv=FALSE, offset=c(300, 0), nodeR=27))


## ----echo=TRUE, tidy=FALSE----------------------------------------------------------------------------------------------------------------
m <- rbind(H1=c(0,           0,           0.5,           0.5          ),
           H2=c(0,           0,           0.5,           0.5          ),
           H3=c("\\epsilon", 0,           0,             "1-\\epsilon"),
           H4=c(0,           "\\epsilon", "1-\\epsilon", 0            ))

graph <- matrix2graph(m)
#graph <- improvedParallelGatekeeping()
graph
substituteEps(graph, eps=0.001)

gMCP(graph, pvalues=c(0.02, 0.04, 0.01, 0.02), eps=0.001)

## ----tidy=FALSE---------------------------------------------------------------------------------------------------------------------------
m <- rbind(H1=c(0, 0, 1, 0, 0),
           H2=c(0, 0, 1, 0, 0),
           H3=c(0, 0, 0, 0.9999, 1e-04),
           H4=c(0, 1, 0, 0, 0),
           H5=c(0, 0, 0, 0, 0))
weights <- c(1, 0, 0, 0, 0)
subgraph1 <- new("graphMCP", m=m, weights=weights)

m <- rbind(H1=c(0, 0, 1, 0, 0),
           H2=c(0, 0, 1, 0, 0),
           H3=c(0, 0, 0, 1e-04, 0.9999),
           H4=c(0, 0, 0, 0, 0),
           H5=c(1, 0, 0, 0, 0))
weights <- c(0, 1, 0, 0, 0)
subgraph2 <- new("graphMCP", m=m, weights=weights)

weights <- c(0.5, 0.5)
graph <- new("entangledMCP", subgraphs=list(subgraph1, subgraph2), weights=weights)

## ----eval=TRUE, tidy=FALSE----------------------------------------------------------------------------------------------------------------

cr <- rbind(H1=c(1   , 0.5 , 0.3 , 0.15),
            H2=c(0.5 , 1   , 0.15, 0.3 ),
            H3=c(0.3 , 0.15, 1   , 0.5 ),
            H4=c(0.15, 0.3 , 0.5 , 1   ))


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------------
#  
#  library(bindata)
#  
#  n1 <- 20
#  n2 <- 1e4
#  
#  pvals <- t(replicate(n2,
#                       sapply(colSums(rmvbin(n1, margprob = c(0.35, 0.4, 0.25, 0.3), bincorr = cr)),
#                              function(x, ...) {binom.test(x, ...)$p.value}, n=n1, alternative="less")
#                       ))
#  

## ----include=FALSE------------------------------------------------------------------------------------------------------------------------
load("pvals.RData")

## ----eval=TRUE----------------------------------------------------------------------------------------------------------------------------

graph <- generalSuccessive(gamma=0, delta=0)
out <- graphTest(pvalues=pvals, graph = graph)
extractPower(out)


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------------
#  # Other example:
#  
#  library(copula)
#  cop <- mvdc(normalCopula(0.75, dim=4), c("norm", "norm", "exp", "exp"), list(list(mean=0.5), list(mean=1), list(rate = 2), list(rate = 3)))
#  x <- rMvdc(250, cop)
#  # ...
#  
#  # Boot example?
#  # Permutation test example?
#  #

## ----echo=TRUE, tidy=FALSE----------------------------------------------------------------------------------------------------------------

# Bonferroni adjustment
G <- diag(2)
weights <- c(0.5,0.5)
corMat <- diag(2)+matrix(1,2,2)
theta <- c(1,2)
calcPower(weights, alpha=0.025, G, theta, corMat)
calcPower(weights, alpha=0.025, G, 2*theta, 2*corMat)


## ----echo=FALSE,results='asis'------------------------------------------------------------------------------------------------------------

cat(graph2latex(generalSuccessive(), scale=0.7, scaleText=FALSE))


## ----echo=TRUE, tidy=FALSE----------------------------------------------------------------------------------------------------------------
graph <- generalSuccessive()
graph

## ----echo=FALSE,results='asis'------------------------------------------------------------------------------------------------------------

cat(graph2latex(result@graphs[[3]], pvalues=pvalues, scale=0.7, scaleText=FALSE))


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
data(hydroquinone)

pvalues <- c() 
x <- hydroquinone$micronuclei[hydroquinone$group=="C-"]
for (dose in c("30 mg/kg", "50 mg/kg", "75 mg/kg", "100 mg/kg", "C+")) {
	y <- hydroquinone$micronuclei[hydroquinone$group==dose]
	result <- wilcox.test(x, y, alternative="less", correct=TRUE)
	pvalues <- c(result$p.value,  pvalues)
}
pvalues

library(coin, quietly=TRUE)
pvalues <- c() 
for (dose in c("30 mg/kg", "50 mg/kg", "75 mg/kg", "100 mg/kg", "C+")) {
	subdata <- droplevels(hydroquinone[hydroquinone$group %in% c("C-", dose),])
	result <- wilcox_test(micronuclei ~ group, data=subdata, distribution="exact")
	pvalues <- c(pvalue(result), pvalues)
}
pvalues

## ----echo=FALSE,fig.width=8,fig.height=6--------------------------------------------------------------------------------------------------
data(hydroquinone)
boxplot(micronuclei~group, data=hydroquinone)

## ----echo=FALSE, results='asis'-----------------------------------------------------------------------------------------------------------

graphs <- list(BonferroniHolm(4),
               parallelGatekeeping(),
               improvedParallelGatekeeping(),
               BretzEtAl2011(),
               HungEtWang2010(),
               HuqueAloshEtBhore2011(),
               HommelEtAl2007(),
               HommelEtAl2007Simple(),
               MaurerEtAl1995(),
               improvedFallbackI(weights=rep(1/3, 3)),
               improvedFallbackII(weights=rep(1/3, 3)),
               #cycleGraph(nodes=paste("H",1:4,sep=""), weights=rep(1/4, 4)),
               fixedSequence(5),
               fallback(weights=rep(1/4, 4)),
               generalSuccessive(weights = c(1/2, 1/2)),
               simpleSuccessiveI(),
               simpleSuccessiveII(),
               truncatedHolm(),
               BauerEtAl2001(),
               BretzEtAl2009a(),
               BretzEtAl2009b(),
               BretzEtAl2009c(),
               FerberTimeDose2011(times=5, doses=3, w=1/2),
               Ferber2011(),
               WangTing2014()
               #Entangled1Maurer2012(),
               #Entangled2Maurer2012()
)

make.url <- function (x) {
    return(gsub("(https?://[^ ]+)", "\\\\url{\\1}", x))
}

texify <- function(x, linebreak="\n\n") {
  x <- make.url(x)
  x <- gsub(linebreak, "\\\\\\\\\n", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("\\\\epsilon", "$\\\\epsilon$", x)
  x <- gsub("\\\\tau", "$\\\\tau$", x)
  x <- gsub("\\\\nu", "$\\\\nu$", x)
  return(x)
}

# first.line("a\nb\nc")
first.line <- function(x) {
  return(strsplit(x, split="\n")[[1]][1])
}

count <-0
for (g in graphs) {
  descr <- attr(g, "descr")
  cat(graph2latex(g, fig=TRUE, scaleText=TRUE, scale=0.7, fig.caption=texify(descr, linebreak = "\n"), fig.caption.short=texify(first.line(descr))))
  count <- count + 1
  if(count%%6==0) cat("\\cleardoublepage\n")
}



