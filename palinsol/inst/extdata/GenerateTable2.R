# Copyright (c) 2012 Michel Crucifix <michel.crucifix@uclouvain.be>

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject

# the following conditions:

# The above copyright notice and this permission notice shall be
# incluudedin all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND INFRINGEMENT
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR

# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# When using this package for actual applications, always
# cite the authors of the original insolation solutions 
# Berger, Loutre and/or Laskar, see details in man pages

# ------------------------------------------------------------------
# R Code developed for R version 2.15.2 (2012-10-26) -- "Trick or Treat"
# ------------------------------------------------------------------ 

GenerateTable2 <- function(Table1, Table4, Table5, sol='BER78')
{
#
# Re-Generates the precession amplitudes and frequencies based on 
# files  provided by Berger 1978 
# (labelled Table 2 in the original publication)
# The present routine is an implementation of the 
# original article by Berger and Loutre : 
# A. Berger and M. F. Loutre, Origine des fr\'equences des \'el\'ements 
# astronomiques intervenant dans l'insolation, 
# Bull. Classe des Sciences, 1-3, 45-106  1990
#

if (sol == 'BER78')
{
  P<- 50.439273  
  zeta <- 3.392506
} else if (sol == 'BER90')
{
  P <- 50.41726176   # cf. eq. (31) in BER90
  zeta <- 1.60075265 # cf. eq. (30) in BER90
}
  ## PsiBar

  g <- Table4$V3
  M <- Table4$V2
  beta <- Table4$V4
  F <- Table5$V2/60./60.*pi/180
  f <- Table5$V3
  delta <- Table5$V4

  ## division in 3 groups, as in Table 13 and Table 14 of Berger & Loutre, Ac. Roy. 1990

  Fre <- c(g+P,outer(g,f,"+")+P,outer(g,f,"-")+P)
  Amp <- c(M,outer(M,F,"*")/2.,outer(M,F,"*")/2. )
  Pha <- c(beta+zeta,outer(beta,delta,"+")+zeta,outer(beta,delta,"-")+zeta)


  ## regroup similar frequencies

  tol <- 0.0001
  Ntrun <- 200.

  Order <- order(abs(Amp),decreasing=TRUE)

  Amp <- Amp[Order]
  Fre <- Fre[Order]
  ## truncates the first 200 terms
  Fre <- Fre[1:Ntrun]
  Amp <- Amp[1:Ntrun]


  N <- length(Fre)
  for (i in 1:(N-1)) {
   for (j in (i+1):N) {
   if (abs(Fre[j] - Fre[i]) < tol) {
    Amp[i] <- Amp[i]+Amp[j]
    Amp[j] <- 0
    } } }

  Order <- order(abs(Amp),decreasing=TRUE)

  Amp <- Amp[Order]
  Fre <- Fre[Order]
  Pha <- (Pha[Order]+180)%%360.
  Per <- 360*60*60/Fre/1000.
  Table2 <<- data.frame(Index=seq(1,length(Amp)),Amp=Amp,Fre=Fre,Pha=Pha,Per=abs(Per[Order])) 

}
