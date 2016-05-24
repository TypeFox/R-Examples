DAISIE_ExpEIN = function(t,pars,M)
{
   pars1 = pars
   lac = pars1[1]
   mu = pars1[2]
   ga = pars1[4]
   laa = pars1[5]
   if(!is.na(pars1[11]))
   {
       M2 = M - roundn(pars1[11] * M)
   } else {
       M2 = M
   }
   A = mu - lac
   B = lac + mu + ga + laa
   C = laa + 2 * lac + ga
   if(t == Inf)
   {
      Imm = ga * M2 / B
      End = (C - ga)/A * Imm
   } else {
      End = M2 * ga * (laa + 2 * lac) * (1/(A * B) - exp(-A*t) / (A * C) + exp(-B*t)/(B * C))
      Imm = M2 * ga / B * (1 - exp(-B * t))
   }
   All = End + Imm
   expEIN = list(End,Imm,All)
   names(expEIN) = c("ExpE","ExpI","ExpN")
   return(expEIN)
}