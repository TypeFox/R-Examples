`foc.color` <-
function(i, pal=0)
  {
    if(missing(pal)) {  pal=0  }
####   setXMCOL() ; XMCOL[c(18, 23, 22, 21, 19, 20, 2)+1]
###  strikeslip_col = 18;
###  rev_oblSS_col = 23;
###  obl_rev_col = 22;
###  reverse_col = 21;
###  norm_oblSS_col = 19;
###  oblq_norm_col = 20;
###  normal_col = 2;
###    fcolors =  XMCOL[c(18, 23, 22, 21, 19, 20, 2)+1 ]
##########  default colorscheme from geotouch
    
    ## fcolors =   c("#C1EEC1","#32FFFF","#B9FFFF","#739BFF","#DFFF61","#FFD732","#FF7878")
    
    
    fcolors =c('#CCFFC3','#C2FFEA','#CBECFF','#D4D9FF','#FFECCF','#E9FFBE','#FFD1DA')
    
    if(identical(pal,1) )
      {
        fcolors=c("DarkSeaGreen", "cyan1","SkyBlue1" , "RoyalBlue" ,"GreenYellow","orange","red")
      }
####  fcolors =   c("DarkSeaGreen", "cyan1",        "SkyBlue1" ,    "RoyalBlue" ,   "GreenYellow",
####   "orange",       "red")
    
    
########  col2rgb(fcolors)
    
#### /* strikeslip */
####  /* rev-obl strk-slp */
####  /* oblique reverse */
####  /* reverse */
####  /* norm-oblq strkslp */
####  /* oblq norm */
####  /* normal */
    return( fcolors[i] )
  }

