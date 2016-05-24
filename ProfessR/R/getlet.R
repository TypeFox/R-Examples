`getlet` <-
function(ggrades, divs)
{
############################
  ############################   OLD WAY
M = length(divs)
  N = M
  AP =  ggrades>=divs[N]-(divs[N]-divs[N-1])/3
  A =   ggrades>= (divs[N]-2*(divs[N]-divs[N-1])/3)& ggrades<divs[N]-(divs[N]-divs[N-1])/3
  AM =   ggrades>= (divs[N-1])& ggrades<divs[N]-2*(divs[N]-divs[N-1])/3

  N = M-1
  BP =  ggrades>=divs[N]-(divs[N]-divs[N-1])/3 & ggrades<divs[N]
  B =   ggrades>= (divs[N]-2*(divs[N]-divs[N-1])/3) & ggrades<divs[N]-(divs[N]-divs[N-1])/3
  BM =   ggrades>= (divs[N-1])& ggrades<divs[N]-2*(divs[N]-divs[N-1])/3

  N = M-2
  CP =  ggrades>=divs[N]-(divs[N]-divs[N-1])/3 & ggrades<divs[N]
  C =   ggrades>= (divs[N]-2*(divs[N]-divs[N-1])/3)& ggrades<divs[N]-(divs[N]-divs[N-1])/3
  CM =   ggrades>= (divs[N-1])& ggrades<divs[N]-2*(divs[N]-divs[N-1])/3

  N = M-3
  DP =  ggrades>=divs[N]-(divs[N]-divs[N-1])/3 & ggrades<divs[N]
  D =   ggrades>= (divs[N]-2*(divs[N]-divs[N-1])/3)& ggrades<divs[N]-(divs[N]-divs[N-1])/3
  DM =   ggrades>= (divs[N-1])& ggrades<divs[N]-2*(divs[N]-divs[N-1])/3

  N = M-4
  E = ggrades>=divs[N-1]&ggrades<divs[N]

  oletts = rep("E", length(ggrades))
  oscores  = rep(0, length(ggrades))

  SCRS = seq(from=100, by=(-4), length=13)
  LETS = c("A+", "A", "A-", "B+", "B", "B-", "C+", "C", "C-", "D+", "D", "D-", "E")

  oscores[AP] = SCRS[1]
  oscores[A] = SCRS[2]
  oscores[AM] = SCRS[3]

  oscores[BP] = SCRS[4]
   oscores[B] = SCRS[5]
   oscores[BM] = SCRS[6]

  oscores[CP] = SCRS[7]
  oscores[C] = SCRS[8]
  oscores[CM] = SCRS[9]

  oscores[DP] = SCRS[10]
  oscores[D] = SCRS[11]
  oscores[DM] = SCRS[12]

 oscores[E] = SCRS[13]-SCRS[13]*(divs[2]-ggrades[E])/divs[2]
  
  oletts[AP] = LETS[1]
  oletts[A] = LETS[2]
  oletts[AM] = LETS[3]

  oletts[BP] = LETS[4]
  oletts[B] = LETS[5]
  oletts[BM] = LETS[6]

  oletts[CP] = LETS[7]
  oletts[C] = LETS[8]
  oletts[CM] = LETS[9]

  oletts[DP] = LETS[10]
  oletts[D] = LETS[11]
  oletts[DM] = LETS[12]

  oletts[E] = LETS[13]
############################
#####################################################  Donna's Way
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
############################
######################################

 return(list(grades=ggrades, lett=letts, scor=scores,  olett=oletts, oscor=oscores, divs=divs, LETS=LETS, SCRS=SCRS))
}

