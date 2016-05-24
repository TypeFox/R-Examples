rowBUTTONS=function (labs, col = 6, pch = 4, cex = 1, boxsize = -1) 
{
  # set some parameters
  usr = par("usr")
  fin = par("fin")
  pin = par("pin")
  ymar.in = (fin[2] - pin[2])/2
  nlabs=length(labs)

  # process boxsize
  if(min(boxsize) < 0 | length(boxsize)==0 | is.na(min(boxsize)) ){
    boxsize = strwidth(labs, units='inches', cex=cex) + 0.15
    boxsize[boxsize<0.4] = 0.4
  }
  if(length(boxsize) == 1){ boxsize = rep(boxsize,nlabs) }

  # define transfer functions between inches and user coordinates
  if(!par('xlog')){
    in2usrx = function(xin,usr,pin){usr[1] + diff(usr[1:2])/pin[1] * xin} # linear x-scale
  } else{
    in2usrx = function(xin,usr,pin){10^usr[1] * (10^usr[2]/10^usr[1])^(xin/pin[1])}# log x-scale
  }
  if(!par('ylog')){
    in2usry = function(yin,usr,pin){usr[3] + diff(usr[3:4])/pin[2] * yin} # linear y-scale
  } else{
    in2usry = function(yin,usr,pin){10^usr[3] * (10^usr[4]/10^usr[3])^(yin/pin[2])}# log y-scale
  }

  # create vector of x locations for boxes (inches from right edge of plot (NOT right edge of figure))
  xin = pin[1] - cumsum(c(0,boxsize[1:(nlabs-1)]) + boxsize)/2 # spacing between box n and box n+1 = boxsize[n]/2 + boxsize[n+1]/2
  bottoms = xin < 0
  tops = xin >= 0
  xin[bottoms] = pin[1] - cumsum( c(0,boxsize[bottoms][min(1,sum(bottoms)-1):(sum(bottoms)-1)]) + boxsize[bottoms]) / 2
  bottoms = bottoms & !(xin < 0)
  xin[xin < 0] = NaN # not enough room for these buttons, so omit them

  # create vector of y locations for boxes (inches from bottom edge of plot (NOT bottom edge of figure))
  yin=0*xin
  yin[tops] = pin[2] + 0.5 * ymar.in
  yin[bottoms] = -0.5 * ymar.in

  # make box location vectors in user coordinates
  px = in2usrx(xin,usr,pin) # box centers
  x1 = in2usrx(xin - 0.9 * boxsize/2, usr, pin) # left edges
  x2 = in2usrx(xin + 0.9 * boxsize/2, usr, pin) # right edges

  py = in2usry(yin,usr,pin) # box centers
  y1 = in2usry(yin - ymar.in * 0.3, usr, pin) # top edges
  y2 = in2usry(yin + ymar.in * 0.3, usr, pin) # bottom edges
  # note: unlike old function, boxes are symmetrical about centers

  # plot buttons
  points(px, py, pch = pch, col = col, xpd = TRUE)
  buttons = list(N = length(px), labs = labs, x1 = x1, x2 = x2, y1 = y1, y2 = y2)
  rect(buttons$x1, buttons$y1, buttons$x2, buttons$y2, border = col, xpd = TRUE)
  text((buttons$x1 + buttons$x2)/2, buttons$y1, labels = labs, pos = 3, col = col, xpd = TRUE, cex = cex)

  return(buttons)
}

    
  
