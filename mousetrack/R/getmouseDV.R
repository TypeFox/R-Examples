## mother function  to compute mouse movement dependent measures
## written by Moreno I. Coco, 2013, (moreno.cocoi@gmail.com)
## GPL > 3; i.e., only non-profit research.

## based on Matlab code by Nick Duran (nduran2@ucmerced.edu)
## with contributions of Michael Spivey and Rick Dale

## Arguments:

## x,y:      coordinate points of the mouse trajectory
## t:        a vector with time indexes
## unit:     the ms sampling x time-point
## counterb: the counterbalancing position of yes-no button to
##           to flip responses accordingly
## dwellfin: region around final click where dwell measures are computed, in pixels
## velajbin: number of timesteps in which to compute/average velocity
## escape:   the amount of pixel to escape when trimming
## escapeinit: region around origin in which initial angle is measured, in pixels

## NOTE: 1) remember to adjust sampen defaults if desired
##       2) remember to twick the counterbalancing (button) variable
##          according to the dataset used

.packageName <- 'mousetrack'

getmouseDV <- function(x, y, t, unit, counterb, refcounterb,
                       dwellfin, velajbin, escape, escapeinit){

    tryCatch({  ## set up an exception handler so that algo does not stop
      
        y = -1 * (y - y[1])  # normalize so all trajs begin at 0
        x = x - x[1]      
        t = t * unit         # indeces of timecourse * the timepoint unit
    
    if (counterb == refcounterb){ # for one of the counterbalanced trial
                          # flip trajectories to put everything on the same side
        x = -1*x;
    }
    
            ## distance (euclidean) ##
            distx = (x[2:length(x)]-x[1:length(x)-1])^2 # (x2-x1)^2
            disty = (y[2:length(y)]-y[1:length(y)-1])^2 # (y2-y1)^2
            DVdist = sum(sqrt(distx+disty));           # distance traveled by one trajectory

        if (DVdist > escape){ ## cases where there is no distance traveled
                              ## or it is smaller than escaping factor
            
            DVtotaltime = tail(t,1)
            
            ## latency time, motion time ##
            
            res = trim0(x,y,escape); # keep only x,y coordinates after "escaping"
                                        # with an given thshd from the move origin
            
            latindx = res$latindx; xtrim = res$xtrim; ytrim = res$ytrim
            DVlatency = t[latindx-1]            # latency time (ms) to begin moving
            DVinmot = t[length(t)] - DVlatency; # motion time (ms) - time spent moving
            # the truncated x trajectory only in latency region (e.g., 100 pixel radius)
            xlatency = x[1:(latindx-1)]
            ylatency = y[1:(latindx-1)]
            
            # only x,y coordinates within N pixel range around the final click "response dwell region"
            res = trim1(x, y, dwellfin);
            latindx3 = res$latindx
            
            DVdwell = t[length(t)] - t[latindx3 + 1]; # response dwell time (ms) - to commit to a final response


            ## dist w/ no latency ##
    
            distlatx = (xtrim[2:length(xtrim)] - xtrim[1:length(xtrim)-1])^2 # (x2-x1)^2
            distlaty = (ytrim[2:length(ytrim)] - ytrim[1:length(ytrim)-1])^2 # (y2-y1)^2
            DVdistinmot = sum(sqrt(distlatx+distlaty)); #distance traveled by one trajectory... outside the 100 pixel region!  


############# velocity, acceleration, and jerk ##################################

            res = velaj(t, x, y, velajbin)         # within 6-step time bins
            veloc = res$veloc; accel = res$accel; jerk = res$jerk 
            
            DVvelmax = max(veloc)
            DVvelmaxstart = t[which(veloc == DVvelmax)[1]]
           # when max velocity occured

            DVaccmax = max(accel) # maximum acceleration
            DVaccmaxstart = t[which(accel == DVaccmax)[1] + velajbin]
            # at what point does maximum acceleration occur
            # DVaccmin = min(accel) % minimum acceleration
            # DVaccminstart = t[which( accel==DVaccmin)[1] + velajbin];
            # at what point is there the greatest change in decreasing velocity

            # DVjerk = 'NA'; # not sure how to interpret jerk

############# properties over trajs and angle ###################################
            
            ##need to make sure not to consider 0 values otherwise R assign pi
            zero = intersect(which(x == 0), which(y == 0) )
            tanpoint = atan2(x,y)

            if (length(tanpoint) > 0){
                tanpoint[zero] = 0
            }
            
            trajAng = tanpoint*(180/pi)
            
            # the trajectories converted to angle plane, with implicit origin
            # at 0,0; compared to y-axis (at 0 degrees), fans out between 0 to
            # 90/-90, with 45/-45 degree angle being direct route to target,
            # see Scherbaum et al, 2010 (Cognition)

            ## arclength after motion was initiated ##

            arclength = cumsum(sqrt(diff(xtrim)^2 + diff(ytrim)^2))
        
            arclengthtotal = arclength[length(arclength)]

            ## calculate max path offset
            maxpathoff = pathoffset(xtrim, ytrim)


######################## xflips ##################################################
        
            DVxflplat = xflip(xlatency)
            DVxflpmot = xflip(xtrim)   # flipping after the first 100 pixel radius
            #  change in direction (flips) for x (latency, motion, dwell) %%%%
            # change in direction for angle trajectories (latency, motion, dwell) %%%% 
            DVaflplat = xflip(trajAng[1:(latindx-1)]) 
            DVaflpmot = xflip(trajAng[latindx:length(trajAng)] )
            
            xlatency3 = x[(latindx3 + 1):length(x)]
            
            # angle flipping in response dwell region
            # the truncated x trajectory only in reponse dwell time
            # (last N pixel radius - see trim1 above)            
            DVxflpdwl = xflip(xlatency3); # x flipping in response dwell region
            DVaflpdwl = xflip(trajAng[(latindx3 + 1):length(trajAng)])

######################### sample entropy ########################################

            # entropy along x-axis 
                
            res = sampen(x, 5, .2, 1, 1)
            DVxe = res$e; DVxse = res$se

            # entropy along y-axis
            res = sampen(y, 5, .2, 1, 1);
            DVye = res$e; DVyse = res$se

            # entropy along trajectory angle

            res = sampen(trajAng, 5, .2, 1, 1);
            DVae = res$e; DVase = res$se            
            
################################################################################# 
 #following measures depend on whether target is on right or left side, properties over trajs and angle 

            ## area under the curve
            DVAUC = areaunder(x,y)
#            DVAUC = auc(x,y)

            ## Below alternative measurement for the area
            ##  to be better tested
#            if ( sign(x[length(x)]) == 1 ){
 # answer on the right side (x > 0), so get deviation toward the left side (x < 0, from point of origin)
                ##take the difference between consecutive points
#                ydiff = diff(y); xdiff = diff(x)
#                ## pad with zero the deviation towards the target
#                xdiff[which(xdiff < 0)] = 0
                
#                DVarea = trapz(y[y >= 0], x[x <= 0]*-1)     
               # grabs portion of trajectory that is x < 0 and y > 0
               # get integral using trapezoidal method
               # DV_areax = trapz(x[x<=0]*-1) # alternative method
               ## trapz in R does not allow integration over a single vector 
#            } else if (sign(x[length(x)]) == -1){
                
               # answer on the left side (x < 0), so get deviation
               # toward the right side (x > 0, from point of origin)

                ##take the difference between consecutive points
 #               ydiff = diff(y); xdiff = diff(x)
                ## pad with zero the deviation towards the target
 #               xdiff[which(xdiff < 0)] = 0
                
 #               DVareadiff = trapz(ydiff[ydiff >= 0], xdiff[xdiff >= 0] )
                
 #               DVarea = trapz(y[y >= 0], x[x >= 0] )
               # grabs portion of trajectory that is x > 0 and y > 0
               # get integral using trapezoidal method
#                DVareax = trapz(max(x,0))
#           }

#            if (is.na(DVarea)){ DVarea = 0 } ## there was a correct in trajectory

### how much PULL to the distractor response - absolute values ###
### pull max timing ###
          
      if ( sign(x[length(x)]) == 1 ){ # rightside
          DVmaxpull = abs(min(x));
          # toward left, either 0 (no pull), or larger
          DVmaxpullstart = t[which( x == min(x))[1]]
          # either 1 (no pull), or larger
      } else if (sign(x[length(x)]) == -1){ # leftside
          DVmaxpull = abs(max(x))
          # toward right, either 0 (no pull), or larger
          DVmaxpullstart = t[which(x == max(x))[1]]
      }

### severity of ANGLE to the distractor response while in motion (no latency region) - absolute values
###  angle max timing ##
            
            motTrajAng = trajAng[latindx:length(trajAng)]
# just the region of trajectory outside of latency region

            if (sign(x[length(x)]) == 1){ # rightside 

                DVmaxang = min(motTrajAng)
            # toward left distractor, should be between 0 and -90
            
                if (DVmaxang >= 0){
            # if DV_maxang is above 0, then in rightside region and no deviation
                    DVmaxang = 0
                    DVmaxangstart = 0
                }  else {
                    DVmaxangstart = t[which(motTrajAng == min(motTrajAng))[1]]
                }
            
            } else if (sign(x[length(x)]) == -1){ # leftside
                DVmaxang = max(motTrajAng)
                # toward right distractory, should be between 0 and 90

                if (DVmaxang <= 0){
                # if DV_maxang is below 0, then in leftside region and no deviation
                    DVmaxang = 0
                    DVmaxangstart = 0
                } else {
                    DVmaxangstart = t[which(motTrajAng == max(motTrajAng))[1]]
                }
            }

            DVmaxang = abs(DVmaxang)
            # get absolute value for angle deviation

            
### severity of INITIAL ANGLE to the distractor response, after committing
## to N pixels (some initial region) - absolute values

            res = trim0(x, y, escapeinit)

            latindx2 = res$latindx; xtrimInit = res$xtrim; ytrimInit = res$ytrim 
            # only x,y coordinates after "escaping" an initial ESCAPE_INIT
            # pixel range around the origin of movement
            DVinitang = trajAng[latindx2-1]
            # what is the angle when trajectory leaves the ESCAPE_INIT
            # region, should be between -90 to 90

            if ((sign(x[length(x)]) == 1) & DVinitang >= 0){
            #if Target is on right, and DV_initang is above 0
            # then no initial deviation

                DVinitang = 0
                
            } else if ( (sign(x[length(x)]) == -1) & DVinitang <= 0){
           # if Target is on left, and DV_initang is below 0,
           # then no initial deviation
                
                DVinitang = 0
            }
            
DVinitang = abs(DVinitang); # get absolute value for initial angle deviation

#################################################################################
            
## Get a few useful potential outlier measures and/or quick analysis
## of quality of trajectories

## max x, y; min x, y; %%% detect wild trajectories

## most positive coordinate (how close to "yes" did they get) ##
#       DVhighx = max(abs(x)) is absolute value needed?
    
DVmaxx = max(x); DVminx = min(x); DVmaxy = max(y); DVminy = min(y);

### percentage of trajectory that is not moving toward or away
## from attractor/distractor, negative movements (y < 0)
            
trajangpos = trajAng;

less1 = which(trajangpos< -90)
trajangpos[less1] = -90

more1 = which(trajangpos > 90)
trajangpos[more1] = 90
            
OLnegmove = (length(less1) + length(more1))/length(trajAng)

## Is motion time longer than latency time? Gets at whether cognitive
## processing is mostly occuring in latency region.   

if (DVinmot < DVlatency){ 
    OL1 = 1
} else {  OL1 = 0 }
      
## Is max velocity inside of latency region? Gets at whether strongest
## commitment to a response occurs in latency region. 

if (which(veloc == DVvelmax)[1] < latindx) {
    OL2 = 1
} else { OL2 = 0 }

##  Is max acceleration inside of latency region? Gets at whether strongest
## commitment to a response occurs in latency region.

if (which(accel == DVaccmax)[1] < latindx){
    OL3 = 1
} else { OL3 = 0 }

## Does trajectory dive below y axis after escapting latency region?
## Gets at whether the trajectory is particulaly wild. 

if (DVmaxang >= 90 | DVinitang >= 90){
    OL4 = 1
} else {  OL4 = 0 }
            
            finres =  list(
                
                DVtotaltime = DVtotaltime, DVlatency = DVlatency,
                DVinmot = DVinmot, DVdwell = DVdwell, 
                DVdist = DVdist, DVdistinmot = DVdistinmot,
                DVvelmax = DVvelmax, DVvelmaxstart = DVvelmaxstart,
                DVaccmax = DVaccmax, DVaccmaxstart = DVaccmaxstart,
                arclengthtotal =  arclengthtotal, maxpathoff = maxpathoff,
                DVxflplat = DVxflplat, DVxflpmot = DVxflpmot,
                DVafllat = DVaflplat, DVaflpmot = DVaflpmot,
                DVxflpdwl = DVxflpdwl, DVaflpdwl = DVaflpdwl,
                DVxe = DVxe, DVxse = DVxse,
                DVye = DVye, DVyse = DVyse,
                DVae = DVae, DVase = DVase,
                trajang = trajAng,
                DVAUC = DVAUC,
                DVmaxpull = DVmaxpull, DVmaxpullstart = DVmaxpullstart,
                DVmaxang = DVmaxang, DVmaxangstart = DVmaxangstart,
                DVinitang = DVinitang,
                DVmaxx = DVmaxx, DVminx = DVminx,
                DVmaxy = DVmaxy, DVminy = DVminy,
                OLnegmove = OLnegmove,
                OL1 = OL1, OL2 = OL2, OL3 = OL3, OL4= OL4
                )
    
        return( finres )
    
        } else { print (
            paste("Distance travelled =", round( DVdist, digits = 4),
                  "< than escape =", escape)
            )
                 
                 return ( vector() )
                 
             }
#    }, warning = function(war) { ## here exception can be handled
    
       # warning handler picks up where error was generated
#       print(paste("WARNING:", war))

#        return ( "check" )
#        return( finres )


    }, error = function(err) {
 
     # warning handler picks up where error was generated
        print(paste("ERROR:",err))
        
        return( vector() )
    })
}
    

##############################################################################
## TODO filter the data as for Jason Friedmann tutorial
    
## create a butterworth filter

## Arguments: set of x,y coordinates
##            butparam = a vector with (n, W)
#             where n = is the polynomial order of the filter
#                   W = is the critical frequency (with a value ranging 0,1       

    
#    res = filtervel(coordinate, butparam, sampling)
#

