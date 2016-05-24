limit <-
function(x,nx,y,ny,conflev,lim,t){
      
      z = qchisq(conflev,1)
      px = x/nx
      score= 0
      while ( score < z){
        a = ny*(lim-1)
        b = nx*lim+ny-(x+y)*(lim-1)
        c = -(x+y)
        p2d = (-b+sqrt(b^2-4*a*c))/(2*a)
        p1d = p2d*lim/(1+p2d*(lim-1))
        score = ((nx*(px-p1d))^2)*(1/(nx*p1d*(1-p1d))+1/(ny*p2d*(1-p2d)))*(nx+ny-1)/(nx+ny) ## added *(nx+ny-1)/(nx+ny)
        ci = lim
        if(t==0) { lim = ci/1.001 }
        else{ lim = ci*1.001 }
        } 
 return(ci)
}

