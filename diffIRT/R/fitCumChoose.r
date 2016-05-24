fitCumChoose=function(nit,order){   # function s(r) from Eq 3 Maydeu - Olivares & Joe 2005 Let op! Deze functie telt de mogelijkheid 'alle items fout' niet mee (i.e., r in Eq 3 loopt vanaf 1 en niet vanaf 0)
  ou=0
if(order==0) return(0)
for(i in 1:order){
  ou=ou+choose(nit,i)
}
return(ou)
}
