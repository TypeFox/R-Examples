`getlet2` <-
function(ggrades, divs)
{
M = length(divs)
  N = M-4
  E = ggrades>=divs[N-1]&ggrades<divs[N]

  letts = rep("E", length(ggrades))
  scores  = rep(0, length(ggrades))

  SCRS = seq(from=100, by=(-4), length=13)
  LETS = c("A+", "A", "A-", "B+", "B", "B-", "C+", "C", "C-", "D+", "D", "D-", "F")

adiff = diff(divs)
adiff/3
J = rep(0, 13)
J[1] = 0

for(i in 2:(M-1))
  {
    j = 3*(i-1)-1
    J[j] = divs[i]+0*adiff[i]/3
    J[j+1] = divs[i]+1*adiff[i]/3
    J[j+2] = divs[i]+2*adiff[i]/3
  }

bdiff = diff(J)

finval = findInterval(ggrades, divs)
Jinval = findInterval(ggrades, J, all.inside = TRUE)

rLETS = rev(LETS)
rSCRS = rev(SCRS)

JLET = rLETS[Jinval]

scores = rSCRS[Jinval]+ 4*(ggrades-J[Jinval])/(J[Jinval+1]-J[Jinval])

letts = JLET

scores[scores>=100] = 100
letts[scores>=100] =LETS[1]

##  cbind(D1$grades[o], D1$lett[o], D1$scor[o], G2$lett[o], format(G2$scor[o], digits=5) )        

 return(list(grades=ggrades, lett=letts, scor=scores, divs=divs, LETS=LETS, SCRS=SCRS))
}

