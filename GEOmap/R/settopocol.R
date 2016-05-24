`settopocol` <-
function()
  {
##########  set up a color map good for topographic displays


cmat = matrix( c(
-20000,  70,      130 ,    180,  -3000,   141,     182,     205  ,  
-3000,   162,     181 ,    205,  -2000   ,  188  ,     210   ,    238     , 
-2000,   202,     225 ,    255 , -0.9 ,   176 ,    196 ,    222  ,   
-0.9,    107,     142,     35   ,  0.1 ,    85 ,     107  ,   47 ,     
0.1,     85 ,     107 ,    47 ,   300 ,    143 ,    188  ,   143   ,  
300 ,    105  ,   139 ,    105 ,  600 ,    180,     238 ,    180 ,    
600 ,    193,     255 ,    193 ,  1000 ,   255 ,    211 ,    155     ,
1000 ,   238,     197,     145,   2000 ,   255 ,    255 ,    255   ,  
2000 ,   255,     255,     255 ,  3500 ,   255 ,    255 ,    255  ),
  ncol = 8, byrow=TRUE)

notes = c(
'#SteelBlue,LightSkyBlue3',
 '#LightSteelBlue3,LightSteelBlue2',
'#LightSteelBlue1,LightSteelBlue',
 '#OliveDrab,DarkOliveGreen',
'#DarkOliveGreen,DarkSeaGreen',
'#DarkSeaGreen4,DarkSeaGreen2',
'#DarkSeaGreen1,burlywood1',
'#burlywood2,White',
'#White,White')


calcol=list(z1=cmat[,1], r1=cmat[,2],g1=cmat[,3],b1=cmat[,4],
     z2=cmat[,5],  r2=cmat[,6],g2=cmat[,7],b2=cmat[,8], note=notes)




    
 coltab = cbind(calcol$r1, calcol$g1, calcol$b1,calcol$r2, calcol$g2, calcol$b2)

    coltab = rbind(coltab, coltab[length(calcol$r1),])


    return(list(calcol=calcol , coltab=coltab))


}

