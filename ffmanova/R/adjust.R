### $Id: adjust.R 45 2006-08-15 13:11:29Z bhm $
# %=============== adjust.m ====================
# % File: adjust.m
# %
# % Call from Matlab:
# %     Xadj = adjust(X,Y)
# %
# % Purpose:
# %     X is "adjusted for Y"
# %     The output matrix Xadj is an orthogonal
# %     basis orthogonal to Y
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright                  %
# % Oyvind Langsrud, MATFORSK  %
# % 2001                       %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function Xadj = adjust(X,Y)
# orthY = myorth(Y);
# orthX = myorth(X);
# Xadj = X(:,[]);
# rankXadj = size(myorth([orthX,orthY]),2) - size(orthY,2);
# if(rankXadj==0)
#    return;
# end;
# Xadj = myorth(orthX - orthY*(orthY'*orthX));
##########################################################
adjust = function(X,Y){
orthY = myorth(Y)
orthX = myorth(X)
rankXadj = myrank(cbind(orthX,orthY)) - dim(orthY)[2]
if(rankXadj==0)
   return(X[,numeric(0),drop = FALSE])
Xadj = myorth(orthX - (orthY%*%(t(orthY)%*%orthX)))
}# end adjust
