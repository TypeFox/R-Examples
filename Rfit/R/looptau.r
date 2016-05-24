looptau = function(delta,abdord,wtord,const,n) {

      m= length(abdord)
      icmax = 0
      inddelta = round(delta*m)
#
# dbug
#     inddelta = 135
# dbug
      guess = abdord[inddelta]
      itmax = m
      y = hstar(abdord,wtord,const,n,guess)
      ikeep = inddelta
      ic = 0
      ierror = 0
      gohome = 0
      istage = 1
#
#  good guess 
#
      if(y == delta){
          ic = 1
          looptau = guess
          gohome = 1
      }
      if(gohome==0){
         if(y > delta){istage = -1}  
      }
#
#  Set the stage
#
while(gohome == 0){
#
#   Iterate
#
      while(ic < .5){
         icmax = icmax + 1
         if(istage == -1){
              if(y > delta){
                   ikeep = ikeep - 1
                   if(ikeep < 1){
                      ic = 1
                      looptau = guess
                      ierror = 1
                   } else {
                      guess = abdord[ikeep]
                      y = hstar(abdord,wtord,const,n,guess)
                   }
              } else {
                 ic = 1 
                 looptau = guess
                 gohome = 1
              }
         } 
         if(istage == 1){
              if(y < delta){
                   ikeep = ikeep + 1
                   if(ikeep > m){
                      ic = 1
                      looptau = guess
                      ierror = 1
                   } else {
                      guess = abdord[ikeep]
                      y = hstar(abdord,wtord,const,n,guess)
                   }
              } else {
                 ic = 1
                 looptau = guess
                 gohome = 1
              }
         }
         icmax = icmax + 1
         if(icmax > itmax){
              ic = 1
              looptau = guess
              ierror = 2
              gohome = 1
         }
      }
}
      list(quan=guess,prob=y,ierror=ierror,icmax=icmax,gohome=gohome)
}
