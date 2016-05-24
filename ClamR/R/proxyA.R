proxyA <- function(ax, ay, xin)
  {

    n = length(ax)
    xout = seq(from=min(ax), to=max(ax), length=4*n)
    hxy = smooth.spline(ax, ay)

    predy = predict(hxy, xout)
    deltax = xout[2]-xout[1]
    fzed = fft(predy$y)
    N = length(predy$x)
    m = floor(N/2)+1   
    f = seq(from=0, to=0.5, length=m)*(1/(deltax))

    Phs = Arg(fzed[2])
    Bphs = sqrt( Mod(fzed[2]) )
    ef = f[2]

    Pos = mean(ay)
    Amp = diff(range(ay))/2
    Prd = 1/ef
    xin = c(Phs,Pos,Amp,Prd)
    
    theY=ay
    myEx=ax
	v1=-2*(max(ax)-min(ax))
	v2=-2*(max(ay)-min(ay))
    
    fr<-function(x)
      {

        Phs = x[1]
        Wmid  = myEx
        Pos    = x[2]
        Amp = x[3]
        Prd = x[4]

        wi = (Amp/2) * sin( (Wmid -Phs)*2*pi/Prd ) + Pos
        ssum = sum( (theY-wi)^2 )
        return(ssum)
      }

### the change is made here###
	  
	#FOUT = optim(xin , fr)
	ui = rbind (c(0,0,1,0),c(0,0,0,1),c(0,0,-1,0),c(0,0,0,-1))
	ci = c(0,0,v2,v1)
	FOUT = constrOptim(xin,fr,grad=NULL,ui,ci) 

    return(FOUT)

  }
