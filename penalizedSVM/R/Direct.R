
Direct<-function(Problem,
                 bounds,
                 # options
                 maxits =     500,         #% maximum of iterations
                 maxevals =   1000,       #% maximum # of function evaluations
                 maxdeep =    100,        #% maximum number of side divisions
                 testflag =   0,          #%  the optimal value is unknown
                 globalmin =  0,          #% minimum value of function
                 ep =         1e-4,       #% global/local weight parameter.
                 tol =        0.01,       #% allowable relative error if f_reach is set
                 showits =    c("none", "final", "all"), #%  plot iteration stats: none, final iteration, all iterations 
                 verbose = 		TRUE, 			#% print  iteration stats: none, final itertation, all iterations 
                 impcons =    0,          #% flag for using implicit constraint handling
                 pert =       1e-6,       #% pertubation for implicit constraint handling
                 # maxflag =   0,          #% set to 1 for max problems, 0 for min problems
                 # sizeconst = 0.5,        #% constant on rectangle size function
                 # distance =  1,          #% 1/0 for distance/volume measure of size
                 # minlength = 1e-4,       #% stop if best rectangle has all sides 1ess than this
                 # minevals =  0          #% but must evaluate at least this many points
                 
                 # plot parameter
                 pdf.name=NULL, 
                 pdf.width=12, pdf.height=12,
                 my.mfrow=c(1,1),  
                 ...){

#% Implementation taken from:
#% D.R. Jones, C.D. Perttunen, and B.E. Stuckman. "Lipschitzian
#% Optimization Without the Lipschitz Constant". Journal of
#% Optimization Theory and Application, 79(1):157-181, October 1993
#%
#%------------------------------------------------------------------%
#
#%-- Initialize the variables --------------------------------------%
lengths<-c <- fc <- vector()
con <- szes <- feas_flags <- vector()
om_lower     <- bounds[,"lower", drop=FALSE]
om_upper     <- bounds[,"upper", drop=FALSE]
fcncounter   <- 0
perror       <- 0
itctr        <- 1
done         <- 0
#g_nargout    <- nargout;  ?????
n            <-nrow(bounds)

#% Determine option values
theglobalmin = globalmin


# %-- New 06/08/2004 Pre-allocate memory for storage vectors
if (testflag == 0){
    lengths    = matrix(0,n,c(maxevals + floor(.10*maxevals)))
    c          = lengths;
    fc         = matrix(0,1,c(maxevals + floor(.10*maxevals))) 
    szes       = fc
    con        = fc
    feas_flags = fc
}

#%-- Call DIRini ---------------------------------------------------%
# define the multy dim hypercube with center point c, fc=f(c) 
#[thirds , lengths, c , fc, con, feas_flags minval,point.xatmin,perror,...
#        History,szes,fcncounter,calltype] =...
#        DIRini(Problem,n,bounds(:,1),bounds(:,2),...
#        lengths,c,fc,con, feas_flags, szes,...
#        theglobalmin,maxdeep,testflag,g_nargout, impcons, varargin{:});
DIRini.list<-.DIRini(Problem,n, a=bounds[,"lower"], b=bounds[, "upper"],
        # p_lengths=lengths, 
        param.names= rownames(bounds),
        c=c,fc=fc, con=con, 
        feas_flags=feas_flags,  szes=szes,
        theglobalmin, maxdeep, testflag, impcons, ...)
      #  theglobalmin, maxdeep, testflag, impcons, fmin, fit.gp, muX, muY)
 
        
thirds <- DIRini.list$thirds
#%  lengths =  length array! will store number of slices in each dimension for
#%-- each rectangle. dimension will be rows; 
#%-- each rectangle will be a column
lengths <-DIRini.list$lengths
c <-DIRini.list$c 
fc <-DIRini.list$fc
con <-DIRini.list$con
feas_flags <-DIRini.list$feas_flags
minval <-DIRini.list$minval
point.xatmin <-DIRini.list$point.xatmin
perror <-DIRini.list$perror
History <-DIRini.list$History
szes <-DIRini.list$szes                 # number of regions
fcncounter <-DIRini.list$fcncounter
calltype <-DIRini.list$calltype
 
ret_minval = minval
ret_point.xatmin = point.xatmin


if(showits !="none" & !is.null(pdf.name)) { 
	pdf(pdf.name, pdf.width, pdf.height)
}
par(mfrow=my.mfrow)

#%-- MAIN LOOP -----------------------------------------------------%
minval = fc[1] + con[1]
while (perror > tol){
   #%-- Create list S of potentially optimal hyper-rectangles
   S <- .find_po(fc=fc+con,
       				 lengths= lengths,
       				 minval=minval, ep=ep, szes=szes)
	# if we don't find potentially hyper-rectanges --> break!
	if (ncol(S)==0) break
	
   #%-- Loop through the potentially optimal hyper-rectangles ------%
   #%-- and divide -------------------------------------------------%
   for (i in 1:ncol(S)){
     # [lengths,fc,c,con,feas_flags,szes,fcncounter,success] = ...
     #     .DIRdivide(bounds(:,1),bounds(:,2),Problem,S(1,i),thirds,lengths,...
     #     fc,c,con,feas_flags,fcncounter,szes,impcons,calltype,varargin{:});
     
     # plot options:
     # don't plot if not requested
     if ( (showits =="none") ) { 
     	showits.flag <-  FALSE
     } else {
	     # plot  last iteration's step  if requested
	     if ((showits =="final")) showits.flag<- ifelse ( (i == ncol(S)), TRUE, FALSE )
	     # plot all iterations' steps  if requested
	     if ( (showits =="all") )showits.flag <-  TRUE
     }
    
     tmp.list.divide<-  .DIRdivide(a=bounds[,1], b=bounds[,2],Problem=Problem,
														index=S[1,i], thirds=thirds, lengths=lengths,
   													fc=fc,c=c,con=con, feas_flags=feas_flags,
          									p_fcncounter=fcncounter, szes=szes,
          									impcons=impcons, calltype=calltype, showits.flag=showits.flag, itctr=itctr, i=i, ...)
                          # impcons=impcons, calltype=calltype, showits.flag=showits.flag, itctr=itctr, i=i, fmin, fit.gp, muX, muY )
	  lengths <- tmp.list.divide$lengths
	  fc  <- tmp.list.divide$fc
	  c  <- tmp.list.divide$c
	  con <- tmp.list.divide$con
	  feas_flags <- tmp.list.divide$feas_flags
	  szes <-  tmp.list.divide$szes
	  fcncounter <-  tmp.list.divide$fcncounter
	  success <- tmp.list.divide$pass     
  }
 
 
  
   #%-- update minval, point.xatmin --------------------------------------%
   # [minval,fminindex] =  min(fc(1:fcncounter)+con(1:fcncounter)); sicher ????
 	 minval<- min(fc + con )
   fminindex<- which.min(fc + con)
   penminval = minval + con[fminindex]
   point.xatmin = (om_upper - om_lower)*c[,fminindex] + om_lower
   if ( (con[fminindex] > 0)|(feas_flags[fminindex] != 0) ){
       #%--- new minval is infeasible, don't do anything
   }else {
       #%--- update return values
       ret_minval <- minval;
       ret_point.xatmin <- point.xatmin;
   } 
 
   #%--see if we are done ------------------------------------------%
   if (testflag == 1){
      #%-- Calculate error if globalmin known
      perror<- ifelse ((theglobalmin != 0), 
      									100*(minval - theglobalmin)/abs(theglobalmin),
      									100*minval )
   }else{
      #%-- Have we exceeded the maxits?
      if (itctr >= maxits){
         if (verbose) print("Exceeded max iterations. Increase maxits")
         done <- 1
      }
      #%-- Have we exceeded the maxevals?
      if (fcncounter > maxevals){
         if (verbose) print("Exceeded max fcn evals. Increase maxevals")
         done <- 1
      }
      if (done == 1)
         perror = -1
   } # end of if else
   
   if (max(as.vector(lengths)) >= maxdeep ){
      #%-- We've exceeded the max depth
      if (verbose) print("Exceeded Max depth. Increse maxdeep")
      perror = -1
   }
   
   #%-- Store History
   History<-rbind(History, 
                 c(itctr, fcncounter, minval))
  
  #%-- New, 06/09/2004
  #%-- Call replaceinf if impcons flag is set to 1
  if (impcons == 1) {
    fc <- .replaceinf(lengths=lengths,c=c,fc=fc,con=con,
                    flags=feas_flags, pert=pert)
  }

  #%-- show iteration stats
  if (verbose)  print(paste("Iter:", itctr, "f_min:", minval, "fn evals:", fcncounter, sep="   "))
    
  itctr  = itctr + 1

} # end  of while (perror > tol)

if(showits !="none" & !is.null(pdf.name)) dev.off()

#%-- Return values      #################
#%-- return x*
final_point.xatmin <- ret_point.xatmin;

#%-- chop off (abschneiden) 1st row of History
History<-History[-1,]

return (list(final_point.xatmin=final_point.xatmin,
							minval =minval, 
							c=c, fc=fc, 
							History=History))
}

###########################################################################################################
        
.DIRini<-function(Problem,n,a,b,
								 param.names=c(1:length(a)),
                 #p_lengths,
                 c, fc, con, feas_flags, szes,
                 theglobalmin,
                 maxdeep,testflag,impcons,...){

#%------------------------------------------------------------------%
#% Function:   DIRini                                               %
#% Written by: Dan Finkel                                           %
#% Created on: 10/19/2002                                           %
#% Purpose   : Initialization of Direct                             %
#%             to eliminate storing floating points                 %
#%------------------------------------------------------------------%
#function [l_thirds,l_lengths,l_c,l_fc,l_con, l_feas_flags, minval,point.xatmin,perror,...
#        History,szes,fcncounter,calltype] = 


# DIRECT begins the optimization by transforming the domain of the problem into the unit
# hyper-cube. The algorithm works in this normalized space, referring to the original space only when
# making function calls. The center of this space is c1, and we begin by fnding f(c1).

	#%-- start by calculating the thirds array
	#%-- here we precalculate (1/3)^i which we will use frequently
	l_thirds<-rep(NA, maxdeep)
	l_thirds[1] <- 1/3
	for (i in 2:maxdeep){
	   l_thirds[i]= (1/3)*l_thirds[i-1];
	}
	#%--lengths=  length array will store # of slices in each dimension for
	#%-- each rectangle. dimension will be rows; 
	#%-- each rectangle will be a column
	
	#%-- first rectangle is the whole unit hyperrectangle
	l_lengths <- matrix(0,n,1);
	rownames(l_lengths)<-  param.names
	
	#%01/21/04 HACK
	#%-- store size of hyperrectangle in vector szes
	szes = 1
	names(szes)<- "start"
	
	#%-- first element of c is the center of the unit hyperrectangle
	#l_c(:,1) = matrix(1/2,n,1)
	# erster Spalte
	l_c <- matrix(1/2,n,1)
	colnames(l_c)<-"start"
	rownames(l_c)<-  param.names
	
	
	#%-- Determine if there are constraints
	calltype = .DetermineFcnType(Problem,impcons);
		
	#%-- first element of f is going to be the function evaluated
	#%-- at the center of the unit hyper-rectangle.
	#%om_point   = abs(b - a)*l_c(:,1)+ a;
	#%l_fc(1)    = feval(f,om_point,varargin{:});
	func.List<- .CallObjFcn(Problem, point.x=l_c[,1, drop=FALSE],a, b, impcons, calltype, ...)
	
	## debugging  
	# func.List<- .CallObjFcn(Problem, point.x=l_c[,1, drop=FALSE],a, b, impcons, calltype, fmin, fit.gp, muX, muY)
	## end of debugging
	
	l_fc<- func.List$fcn_value
	l_con<- func.List$con_value
	l_feas_flags<- func.List$feas_flag
	fcncounter = 1
		
	#%-- initialize minval and point.xatmin to be center of hyper-rectangle !  (NOT in the original intervals!!!!)
	point.xatmin = l_c[,1, drop=FALSE]      
	minval   = l_fc[1]
	if (testflag == 1){
	    if (theglobalmin != 0){
	        perror = 100*(minval - theglobalmin)/abs(theglobalmin);
	    }else{
	        perror = 100*minval;
	    }
	}else {
	   perror = 2
	}
	#%-- initialize History
	History<-t(matrix(c(0,0,0)))  
  colnames(History)<- c("Iteration Nr", "Function Count", "f_min"  )
  
return(list( thirds=l_thirds, lengths=l_lengths,
						 c=l_c, fc=l_fc, con=l_con, 
					 	 feas_flags=l_feas_flags,
		    		 minval=minval,point.xatmin=point.xatmin,perror=perror,
		         History=History, szes=szes,
		         fcncounter=fcncounter,calltype=calltype ))
}

####################################################################################################
 
.find_po<-function(fc,lengths,minval,ep,szes){

#%--------------------------------------------------------------------%
#% Function   :  find_po                                              %
#% Written by :  Dan Finkel                                           %
#% Created on :  10/19/2002                                           %
#% Purpose    :  Return list of potentially optimal hyper-rectangles  %
#%--------------------------------------------------------------------%
#function rects = find_po(fc,lengths,minval,ep,szes)

#%-- 1. Find all rects on hub
# diff_szes = sum(lengths,1);
diff_szes = colSums(lengths) # col sum or row sums? nicht sicher ??????
tmp_max = max(diff_szes)
j=1
hull<-vector()   ##??? nicht sicher !!!!
sum_lengths = colSums(lengths)

for (i in 1:(tmp_max+1)){
    tmp_idx <- which(sum_lengths == (i-1))
    if(length(tmp_idx)>0){
	    tmp_n <- min(fc[tmp_idx])
	    hullidx <- which.min(fc[tmp_idx]) 
	    if (length(hullidx) > 0 ){
	        hull[j] <- tmp_idx[hullidx]
	        j=j+1;
	        #%-- 1.5 Check for ties
	        ties <- which(abs(fc[tmp_idx]- tmp_n) <= 1e-13)
	        if (length(ties) > 1){
	            mod_ties <- which(tmp_idx[ties] != hull[j-1])
	            hull <- c(hull, tmp_idx[ties[mod_ties]])
	            j <- length(hull)+1;
	        } # end of the if length(ties) > 1
	    } # end of the  if length(hullidx) > 0 
    } # end of if length(tmp_idx)>0
} # end of for 


#%-- 2. Compute lb and ub for rects on hub
lbound <- .calc_lbound(lengths,fc,hull,szes)
ubound = .calc_ubound(lengths,fc,hull,szes)

#%-- 3. Find indeces of hull who satisfy
#%--    1st condition
maybe_po <- which(lbound-ubound <= 0)

#%-- 4. Find indeces of hull who satisfy
#%--    2nd condition
t_len  <- length(hull[maybe_po])
if (minval != 0){
    po = which( ( (minval-fc[hull[maybe_po]])/abs(minval) +
                   szes[hull[maybe_po]] * ubound[maybe_po]/abs(minval)  ) >= ep)
}else {
    po = which (( fc[hull[maybe_po]] - 
                  szes[hull[maybe_po]] * ubound[maybe_po] ) <= 0)
} 

final_pos  <- hull[maybe_po[po]]


rects <- rbind(final_pos, "szes"=szes[final_pos])  
return(rects)
}

####################################################################################################
.calc_ubound<-function(lengths,fc,hull,szes){
#%------------------------------------------------------------------%
#% Function   :  calc_ubound                                        %
#% Written by :  Dan Finkel                                         %
#% Created on :  10/19/2002                                         %
#% Purpose    :  calculate the ubound used in determing potentially %
#%               optimal hrectangles                                %
#%------------------------------------------------------------------%
#function ub = calc_ubound(lengths,fc,hull,szes)
	ub<- vector()
	hull_length  <- length(hull)
	hull_lengths <- lengths[,hull, drop=FALSE]
	for (i in 1:hull_length){
	    tmp_rects = which(colSums(hull_lengths)< sum(lengths[,hull[i] ] ))
	    if (length(tmp_rects) > 0){
	        tmp_f     <- fc[hull[tmp_rects]]
	        tmp_szes  <- szes[hull[tmp_rects]]
	        tmp_ubs   <- (tmp_f - fc[hull[i]])/(tmp_szes - szes[hull[i]])
	        ub[i]     <- max(tmp_ubs);
	    }else{
	        ub[i]     <- 1.976e14;
	    }
	}
return(ub)
}
######################################################################################################
.calc_lbound<-function(lengths,fc,hull,szes){
#%------------------------------------------------------------------%
#% Function   :  calc_lbound                                        %
#% Written by :  Dan Finkel                                         %
#% Created on :  10/19/2002                                         %
#% Purpose    :  calculate the lbound used in determing potentially %
#%               optimal hrectangles                                %
#%------------------------------------------------------------------%
##function lb = calc_lbound(lengths,fc,hull,szes)
	lb<- vector()
	hull_length  <- length(hull)
	hull_lengths <- lengths[,hull, drop=FALSE]
	for (i in 1:hull_length){
	    tmp_rects = which(colSums(hull_lengths)> sum(lengths[,hull[i] ] ))
	    if (length(tmp_rects) > 0){
	        tmp_f     <- fc[hull[tmp_rects]]
	        tmp_szes  <- szes[hull[tmp_rects]]
	        tmp_lbs   <- (fc[hull[i]]-tmp_f)/(szes[hull[i]]-tmp_szes)
	        lb[i]     <- max(tmp_lbs);
	    }else{
	        lb[i]     <- -1.976e14;
	    }
	}
	return(lb)
}

#lbound <- calc_lbound(lengths,fc,hull,szes)

######################################################################################################




.DIRdivide<-function(a,b,Problem,index,thirds,
										lengths, fc, c, con, 
								    feas_flags, p_fcncounter, szes,
								    impcons, calltype,showits.flag=TRUE, itctr="", i="", ...){
								    
#%------------------------------------------------------------------%
#% Function   :  DIRdivide                                          %
#% Written by :  Dan Finkel                                         %
#% Created on :  10/19/2002                                         %
#% Purpose    :  Divides rectangle i that is passed in              %
#%------------------------------------------------------------------%
#function [lengths,fc,c,con,feas_flags,szes,fcncounter,pass] = ...
#    DIRdivide(a,b,Problem,index,thirds,p_lengths,p_fc,p_c,p_con,...
#    p_feas_flags,p_fcncounter,p_szes,impcons,calltype,varargin)

fcncounter <- p_fcncounter

#%-- 1. Determine which sides are the largest #########################
li     <- lengths[,index]  
biggy  <- min(li)
ls     <- which(li==biggy)
lssize <- length(ls)

#%-- 2. Evaluate function in directions of biggest size #########################
#%--    to determine which direction to make divisions
oldc       <- c[,index,drop=FALSE]
delta      <- thirds[biggy+1]   

# add or create new centers? c_old +/- delta*e_i --> 4 new centers(left, right, up, down)!   
#newc_left  <- oldc[,matrix(1, nrow=1,ncol=lssize)];  ### ??? 
#newc_right <- oldc(:,ones(1,lssize));                ### ???  

newc_left  <- newc_right <- matrix(rep(oldc, lssize), ncol=lssize  ) 

# initialize
f_left <- con_left <- fflag_left   <- rep(NA, lssize)
f_right <- con_right <- fflag_right    <- rep(NA, lssize) 

# for each dimention (parameter) in ls create new centers : left and right from the old center
for (i in 1:lssize){
    lsi  <- ls[i]
    
    # c_i +/- delta*e_i  , e_i has at the ith position 1, rest=0
    newc_left[lsi,i]  = newc_left[lsi,i] - delta;
    newc_right[lsi,i] = newc_right[lsi,i] + delta;
    
    # f(new_centers_left)
    func.left.list<- .CallObjFcn(Problem, point.x=newc_left[,i,drop=FALSE],a, b, impcons, calltype, ...)
		f_left[i]<-  func.left.list$fcn_value
		con_left[i]<-  func.left.list$con_value
		fflag_left[i]<-  func.left.list$feas_flag
    
    # f(new_centers_right)
    func.right.list<- .CallObjFcn(Problem, point.x=newc_right[,i,drop=FALSE],a, b, impcons, calltype, ...)
		f_right[i]<-  func.right.list$fcn_value
		con_right[i]<-  func.right.list$con_value
		fflag_right[i]<-  func.right.list$feas_flag
    
		# counter := add 2 	
    fcncounter <- fcncounter + 2
}
  
#%---- 2.1 Calculate w - min of function in new centers #############
#w = [min(f_left, f_right)' ls]; 
# like in DIRECTUserGuide:
# w_i := min( f_right, f_left ) for i in 1:N
# best function valueS ! in the largest space
# it means for each dimention find separat min ! 
w = apply(cbind(f_left, f_right), 1, min)


#%-- 3. Sort w for division order #########################
tmp.sort<-sort(w,index.return=TRUE)
V<-tmp.sort$x; order<-tmp.sort$ix

#%-- 4. Make divisions in order specified by order #########################
for (i in 1:length(order) ){

   newleftindex  = p_fcncounter+2*(i-1)+1
   newrightindex = p_fcncounter+2*(i-1)+2
   #%-- 4.1 create new rectangles identical to the old one ########
   oldrect <- lengths[,index, drop= FALSE]
   lengths<- cbind(lengths, oldrect, oldrect)

   #%-- old, and new rectangles have been sliced in order(i) direction
   lengths[ls[order[i]],newleftindex ] <- lengths[ls[order[i]],index] + 1
   lengths[ls[order[i]],newrightindex] <-  lengths[ls[order[i]],index]  + 1;
   lengths[ls[order[i]],index]         <-  lengths[ls[order[i]],index]  + 1;

   #%-- add new columns to c
   c<-cbind(c, newc_left[,order[i]], newc_right[,order[i]] )
   colnames(c)[(ncol(c)-1):ncol(c)]<- c("left", "right")
   
   #%-- add new values to fc
   fc<- c(fc, f_left[order[i]], f_right[order[i]] )
   
   #%-- add new values to con
   con<- c(con, con_left[order[i]], con_right[order[i]] )
  
   #%-- add new flag values to feas_flags
   feas_flags<- c(feas_flags, fflag_left[order[i]], fflag_right[order[i]] )
     
   #%-- 01/21/04 Dan Hack
   #%-- store sizes of each rectangle    ### sicher ????  A -vector or matrix
   #n = norm(A), a matrix  returns the largest singular value (!) of A, max(svd(A)$d).
   # if A vector:  norm(A)=sum(abs(A).^p)^(1/p), p=2
   # szes(1,newleftindex)  = 1/2*norm((1/3*ones(size(lengths,1),1)).^(lengths(:,newleftindex)));
   
		#	   # if matrix
		#		 tmp<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,newleftindex, drop=FALSE])
		#		 max(svd(tmp)$d)
		#		 # if vector
		#		 tmp<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,newleftindex])
		#	   sum(abs(tmp)^2)^(1/2)
		#		 ## the same result :-) funny!!!! use matrix!
		
		tmp.szes.l<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,newleftindex, drop=FALSE])
		tmp.szes.l<- 1/2*max(svd(tmp.szes.l)$d)	
		
		tmp.szes.r<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,newrightindex, drop=FALSE])
		tmp.szes.r<- 1/2*max(svd(tmp.szes.r)$d)	
		   
    szes<-c(szes, tmp.szes.l, tmp.szes.r )   ## sicher ???
   names(szes)[(length(szes)-1):length(szes)]<- c("left", "right")
} #end of for
  

## plot old and new centers ####
if (showits.flag ){
	my.col<- rep("black", ncol(c))
	my.col[grep("left",colnames(c))]<- "blue"
	my.col[grep("right",colnames(c))]<- "red"
	
	if (nrow(c)==1){
		toPlot<- rbind(c, fc) 
		plot(as.data.frame(t(toPlot)), pch=20, type="p", xlim=c(0,1), col=my.col, main=paste("Iteration", itctr, "step", i  ) )
		text(x=toPlot[1,],y=toPlot[2,], round(fc,1), pos=1, cex=0.3)		
		text(x=toPlot[1,],y=toPlot[2,], c(1:ncol(c)), pos=3, cex=0.3, col="#008080")	
		legend("topright", c("old center", "new left", "new right", "center's number"), fill=c("black", "blue", "red", "#008080"), cex=0.5)
	} # end for 1 tuning parameter
	
	# for 2 tuning parameters
	if (nrow(c)==2){
	toPlot<- c
	plot(as.data.frame(t(toPlot)), pch=20, type="p", xlim=c(0,1), ylim=c(0,1), col=my.col, main=paste("Iteration", itctr, "step", i  )  )
	text(x=toPlot[1,],y=toPlot[2,], round(fc,1), pos=1, cex=0.3)		
	text(x=toPlot[1,],y=toPlot[2,], c(1:ncol(c)), pos=3, cex=0.3, col="#008080")	
	legend("topright", c("old center", "new left", "new right", "center's number"), fill=c("black", "blue", "red", "#008080"), cex=0.5)
	} # end for 2 tuning parameters
}
## end of plot old and new centers #### 		


tmp.szes.ind<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,index, drop=FALSE])
tmp.szes.ind<- 1/2*max(svd(tmp.szes.ind)$d)	
szes[index] <- tmp.szes.ind
pass = 1


return(list(lengths=lengths,fc=fc,c=c,con=con, feas_flags=feas_flags,
			 szes=szes,fcncounter=fcncounter,pass=pass))

}

####################################################################################################

 .CallConstraints<- function(Problem,x,a,b,...){  
#%------------------------------------------------------------------%
#% Function   :  CallConstraints                                    %
#% Written by :  Dan Finkel                                         %    OK
#% Created on :  06/07/2004                                         %
#% Purpose    :  Evaluate Constraints at pointed specified          %
#%------------------------------------------------------------------%
#function ret_value = CallConstraints(Problem,x,a,b,varargin)

	#%-- Scale variable back to original space
	point = abs(b - a)*x+ a;
	
	ret_value = 0;
	if ( ("constraint" %in% names(Problem))  ){
	    if (length(Problem$constraint)>0){ # slot constraints exists, look at their slots
	        for (i in 1:Problem$numconstraints){
	            if (length(Problem$constraint[[i]]$func) == length(Problem$f)){
	                if ( Problem$constraint[[i]]$func == Problem$f ){                    
	                    #%-- Dont call constraint; value was returned in obj fcn
	                    con_value = 0;
	                }else{
	                    con_value = eval(parse(text=Problem$constraint[[i]]$func))(point, ...)
	                } 
	            }else{
	               con_value = eval(parse(text=Problem$constraint[[i]]$func))(point, ...)
	            }
	            if (con_value > 0) {
	                #%-- Infeasible, punish with associated pen. param
	                ret_value = ret_value + con_value * Problem$constraint[[i]]$penalty
	            } 
	        } # end of for
	   } # end of (length(Problem$constraint)>0)
	}
	return (ret_value)
}

###########################################################################################################

.CallObjFcn<-function(Problem,point.x,a,b,impcon,calltype,...){
#%------------------------------------------------------------------%
#% Function   :  CallObjFcn                                         %
#% Written by :  Dan Finkel                                         %   OK
#% Created on :  06/07/2004                                         %
#% Purpose    :  Evaluate ObjFcn at pointed specified               %
#%------------------------------------------------------------------%
#function [fcn_value, con_value, feas_flag] = ...
#    CallObjFcn(Problem,point.x,a,b,impcon,calltype,varargin)

## point.x = vector of values for tuning parametr(s)  
## in arguments of Problem function: point - vector of tuning parameters at the first place

	con_value = 0;
	feas_flag = 0;
	
	#%-- Scale variable back to original space
	point = abs(b - a)*point.x+ a
	
	if (calltype == 1){
	    #%-- No constraints at all
	    # find the functions value at 'point'
	    fcn_value = eval(parse(text=Problem$f))(point, ...)
	    
	    ## debug
	    # fcn_value = eval(parse(text=Problem$f))(point, fmin, fit.gp, x.svm,y.svm)
	    # end of debug
	}
	if (calltype == 2){
	    #%-- f  and   all constraints
	    tmp.list2<- eval(parse(text=Problem$f))(point, ...)
	    fcn_value<-tmp.list2$fcn_value
	    cons<-tmp.list2$cons
	   # [fcn_value, cons] = feval(Problem.f,point,varargin{:});
	    for (i in 1:length(cons)){
	        if (cons > 0){
	         con_value <- con_value + Problem$constraint[[i]]$penalty * cons(i);
	        }
	    }
	}
	if (calltype == 3){    
	    #%-- f returns no constraint values
	    fcn_value <- eval(parse(text=Problem$f))(point, ...)
	    con_value <- .CallConstraints(Problem,point.x,a,b,...);
	}
	if (calltype == 4){  
	    #%-- f returns feas flag
	    tmp.list4<- eval(parse(text=Problem$f))(point, ...)
	    fcn_value<-tmp.list4$fcn_value
	    feas_flag<-tmp.list4$feas_flag
	}
	if (calltype == 5){
	    #%-- f returns feas flags, and there are constraints
	    tmp.list5<- eval(parse(text=Problem$f))(point, ...)
	    fcn_value<-tmp.list5$fcn_value
	    feas_flag<-tmp.list5$feas_flag
	     con_value <- .CallConstraints(Problem,point.x,a,b,...);
	}
	
	if (feas_flag == 1){
		fcn_value = 10^9
	  con_value = 0
	}
	return(data.frame("fcn_value"=fcn_value, "con_value"=con_value, "feas_flag"=feas_flag))
}

##################################################################################################################

 .replaceinf<- function(lengths,c,fc,con,flags,pert){

#%------------------------------------------------------------------%
#% Function   :  replaceinf                                         %
#% Written by :  Dan Finkel                                         %
#% Created on :  06/09/2004                                         %
#% Purpose    :  Assign R. Carter value to given point              %
#%------------------------------------------------------------------%
# 

#%-- Initialize fcn_values to original values
fcn_values <- fc

#%-- Find the infeasible points
infeas_points <- which(flags == 1)

#%-- Find the feasible points
feas_points   = which(flags == 0)

#%-- Calculate the max. value found so far
maxfc<- ifelse( length(feas_points)>0, 
								max(fc[feas_points] + con[feas_points]),
								max(fc + con) )

if (length(infeas_points)>0){
	for (i in 1:length(infeas_points) ){
	
	    if (length(feas_points)==0){
	        #%-- no feasible points found yet
	        found_points <-found_pointsf <- vector()
	        index <- infeas_points[i];
	    } else {
	        index = infeas_points[i]
	
	        #%-- Initialize found points to be entire set
	        found_points  <- c[,feas_points, drop=FALSE ]
	        found_pointsf <- fc[feas_points] + con[feas_points]
	
	        #%-- Loop through each dimension, and find points who are close enough
	        for (j in 1:nrow(lengths) ){
	            neighbors <- which(abs(found_points[j,] - c[j,index]) <=  3^(-lengths[j,index]))
	            if (length(neighbors)>0 ){
	                found_points  <- found_points[,neighbors]
	                found_pointsf <- found_pointsf[neighbors]
	            } else{
	                found_points <-found_pointsf <- vector()
	                break
	            } 
	        } 
	    } # end of if else
	
	    #%-- Assign Carter value to the point
	    if (length(found_pointsf)>0) {
	        #%-- assign to index the min. value found + a little bit more
	        fstar <- min(found_pointsf);
	        fcn_values[index] <- ifelse ((fstar != 0), 
	        															fstar + pert*abs(fstar),
	        															fstar + pert*1 )
	    }else {
	        fcn_values(index) = maxfc+1
	        maxfc             = maxfc+1
	    }
	} # end of for
} # end of if  (length(infeas_points)>0)
return (fcn_values)
} 

####################################################################################################

.DetermineFcnType<-function(Problem,impcons){
#%------------------------------------------------------------------%
#% Function   :  DetermineFcnType                                   %
#% Written by :  Dan Finkel                                         %
#% Created on :  06/25/2004                                         %
#% Purpose    :  Determine how constraints are handled              %
#%------------------------------------------------------------------%

	retval = 0;
	if ( !("constraint" %in% names(Problem)) & (!impcons) ){
	    # 1. %-- No constraints at all  and  (no implicit constraint)
	    retval = 1
	}
	if ("constraint" %in% names(Problem)){
		#%-- There are explicit constraints. Next determine where
    #%-- they are called
    if (length(Problem$constraint)>0){ # slot constraints exists, look at their slots
		        if (length(Problem$constraint[[1]]$func) == length(Problem$f))
		            #%-- Constraint values may be returned from objective
		            #%-- function. Investigate further
		            if (Problem$constraint[[1]]$func == Problem$f ){  
		                #%-- f returns constraint values
		                retval = 2
		            }else {
		               # %-- f does not return constraint values
		               retval = 3;
		            } # end of There are explicit constraints
		} else{ 
			# thre is a EMPTY slot named constraints     
      if (impcons){
		  	retval = 0;
		  }else{
		  	retval = 1;
		  } 
		}# end of no constraints
	} #no slot named constraints    
	
	if (impcons){
		    if (retval== 0 ){
		        #%-- only implicit constraints
		        retval = 4;
		    }else{
		        #%-- both types of constraints
		        retval = 5;
		    }
	}
	
	return(retval)
}


#%------------------------------------------------------------------%
#% Versions  : 1.0 - 1st successful implemenation of DIRect
#%           : 2.0 - Removed floating point arithmetic
#%                   duplicated Table 5 of Jones et al.
#%           : 2.1 - increased speed by storing size calcs.
#%           : 2.2 - utitilized linked lists to increase speed
#%           : 2.3 - rewrote ubound to increase speed
#%           : 2.4 - rewrote lbound to increase speed
#%           : 2.5 - removed call to calcsize
#%           : 2.6 - added check_for_ties
#%           : 2.7 - rewrote check_for_ties to compare fp correctly
#%           : 2.8 - changed output arguments, rewrote help
#%           : 3.0 - simplified input/output. Put on web.
#%           : 3.1 - Performanced Tuned! Tremendous speed increase
#%           : 3.2 - Removed llists; performance tuned
#%                   Many thanks to Ray Muzic and Paul Fackler
#%                   for their suggestions to improve this code
#%           : 4.0 - Sped up code, and added 2 constraint handling
#%                   mechanisms.
#%------------------------------------------------------------------%
