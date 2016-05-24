# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


ipainfo = function(wanted){
  out = list()
  out$symbol =  -as.hexmode(c('0069','0079','0268','0289','026F','0075',
  '026A','028F','028A','0065','00F8','0258','0275','0264','006f','0259',
  '025B','0153','025C','025E','028C','0254','00E6','0250','0061','0276',
  '0251','0252','025A','025D'))
  out$description = c("close front unrounded","close front rounded",
  "close central unrounded","close central rounded","close back unrounded",
  "close back rounded","lax close front unrounded","lax close front rounded", 
  "lax close back rounded","front close-mid unrounded","front close-mid rounded",
  "close-mid schwa","rounded schwa","close-mid back unrounded",
  "close-mid back rounded","schwa" ,"open-mid front unrounded",
  "front open-mid rounded","open-mid central","open-mid central rounded",
  "open-mid back unrounded","open-mid back rounded", "raised open front unrounded",
  "open-mid schwa","front open unrounded","front open rounded","open back unrounded",        
  "open back rounded","rhotacized schwa","rhotacized open-mid central")
  out$num = 1:30
  tmp1 = c(4.0,4.0,4.0,4.0,4.0,4.0,3.5,3.5,3.5,3.0,3.0,3.0,3.0,
  3.0,3.0,2.5,2.0,2.0,2.0,2.0,2.0,2.0,1.5,1.5,1.0,1.0,1.0,1.0,2.5,2.5)
  tmp2 = c(1.00,1.00,2.00,2.00,3.00,3.00,1.50,1.50,2.50,1.00,1.00,2.00,2.00,3.00,
  3.00,2.25,1.00,1.00,2.00,2.00,3.00,3.00,1.00,2.00,1.00,1.00,3.00,3.00,2.75,1.75)
  tmp3 = c(0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,0,0,1,0,1,0,0)
  out$chooseplot = data.frame (height = tmp1, frontness = tmp2, rounded = tmp3)
  out$sampa = c('i','y','1','}','M','u','I','Y','U','e','2','@\\','8','7','o',
  '@','E','9','3','3\\','V','O','{','6','a','&','A','Q','@\'','3\'')
  invisible (out) 
}
