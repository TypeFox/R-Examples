LRFN.plot <-
function(M, Left.fun = NULL, Right.fun = NULL, xlim = NULL, ylim = NULL, lwd = NULL, lty = NULL, col = NULL, add = NULL)
{
if ( length(col) == 0) {col=1}

if ( messages(M) != 1 )  { return( messages(M) ) }

m = M[1]
m_l = M[2]
m_r = M[3]

x = seq(xlim[1], xlim[2], len = 1000)

if ( M[4] == 0 ) { y = Left.fun((m-x)/m_l) * (x<=m) + Right.fun((x-m)/m_r) * (m<x) }
  else if ( M[4] == 1 ) { y = Right.fun((m-x)/m_l) * (x<=m) + Left.fun((x-m)/m_r) * (m<x) }
  else if ( M[4] == 0.5 ) { y = Left.fun((m-x)/m_l) * (x<=m) + Left.fun((x-m)/m_r) * (m<x) }
  else { return( noquote( paste0("The fourth element of each LR fuzzy number must be 0 or 0.5 or 1!" ) ) ) }

y = y * (0<=y & y<=1)
par(new = add)
return( plot.default(x, y, type = 'l', xlim = xlim, ylim = ylim, lwd = lwd, lty = lty, col = col, new = add) )
# if (Left.fun == Right.fun + 100 ) print(2) # Yek jomleye alaki choon CRAN majburet karde bud ke 
# argument-haye Left.fun va Right.fun ra baraye tabe alaki tarif koni va az an-ha estefade koni
}
