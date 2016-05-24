XZ_BSPLINE.f <-
function(x, par.lme, ...){
#   ************************************************************************************************

		knt_x 	= par.lme$knt_x
		ord_x 	= par.lme$ord_x

		knt_z 	= par.lme$knt_z
		ord_z 	= par.lme$ord_z

   		BS_x  	= BSplines(knots=knt_x,ord = ord_x, der = 0, x = x)   ;	# fix(BS.x)
   		BS_z  	= BSplines(knots=knt_z,ord = ord_z, der = 0, x = x)   ;

   		X	= BS_x[1:nrow(BS_x),1:ncol(BS_x)-1,drop=F]     ;	# --> Ey(x=1)=0
    	Z	= BS_z[1:nrow(BS_z),1:ncol(BS_z)-1,drop=F]     ;

		return(list(x = x, X = X, Z = Z))
	}
