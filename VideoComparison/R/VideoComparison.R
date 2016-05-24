VideoComparison<-function(mm,m2,stp=10,nprocesses=0){
  #
  # Match curvature function
  part2part<-function(x,y,lmin=100,nprocesses=0) {
   
    fastpart2partcc<-function(x,y,lmin=100,score=0,nprocesses=0){
      ## Returns the max score among the scores of the crosscorrelation between y and the different subvectors of x
      ## score is x's subvector's length times the normalized xcorr between x's subvector and y value 
      ## WARNING: fastpart2partcc(X,Y) != fastpart2partcc(Y,X)
      ## ARGS:
      #   x,y: vectors
      #   lmin: minimal subvector length (the absolute minimal length is 2, as subvectors of length 1 have variance == 0)
      #   score: minimal score
      #   nprocesses: Number of processes that should be spawned 
      #       <1: Leave it to the function's choice (requires parallel package)
      #        1: Only one process (defaults to this option if no parallel package is found)
      #       >1: This amount of processes (requires parallel package)


      #############################################################################################################################
      # Author (Optimization: FFT): Jose M Perez Ramos (josem.perez.ramos@gmail.com) (Working for ASASEC)
      # Author (Optimization: N^3, parallelization): Jose M Perez Ramos (josem.perez.ramos@gmail.com) (Working for ASASEC)
      # Date: 2013.11.15
      # Complexity: O(N^3)

      # This function calculates the relation between two vectors
      # The relation is measured by a score where 
      #   score = max(normalized_xcorr(subvector(x),y)*length(subvector(x))) for every subvector of x with at least min_length values

      # To find this score, every x subvector is xcorr'd with y. The complexity of this approach is O(N^4) with a straigthforward method and O(N^3*logN)
      # if the xcorr is calculated using FFT (Notice that the straightforward method fares better with short subvectors because of FFT requiremente of
      # having both vectors the same length)

      # So where's the trick to achieve a N^3 algorithm? 
      # -The delay var in the xcorr and the requirement of substract the mean of the subvector x to the values of that subvector

      # X & Y don't have to be the same length
      # X: [----------------------------------------------------------------------------------------]
      # Y: [----------------------------------------------------------------------------------------]

      # Straigthforward method
      #         X: [----------------------------------------------------------------------------------------]
      # Subvector1(X): [................................--------------------------------........................]
      # Subvector2(X): [................................--------------------------------N.......................]

      # XCorr with subvector1:
      #   Xcorr[i]=   
      #     Subvector1(X): [................................OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO........................]
      #     Y :            [--------------------------------OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO------------------------]
      #   Xcorr[i+1]=   
      #     Subvector1(X): [................................11111111111111111111111111111111........................]
      #     Y :           [---------------------------------11111111111111111111111111111111-----------------------]
      #   Xcorr[i+2]=   
      #     Subvector1(X): [................................22222222222222222222222222222222........................]
      #     Y :          [----------------------------------22222222222222222222222222222222----------------------]
      #   ...
      # XCorr with subvector2:
      #   Xcorr[i]=   
      #     Subvector2(X): [................................OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOON.......................]
      #     Y :            [--------------------------------OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOON-----------------------]
      #   Xcorr[i+1]=   
      #     Subvector2(X): [................................11111111111111111111111111111111N.......................]
      #     Y :           [---------------------------------11111111111111111111111111111111N----------------------]
      #   Xcorr[i+2]=   
      #     Subvector2(X): [................................22222222222222222222222222222222N.......................]
      #     Y :          [----------------------------------22222222222222222222222222222222N---------------------]
      #   ...   
        
      # Notice how for subvector2 xcorr, the values marked had already been calculated (disregarding the fact that it's needed to substract the mean for every subvector and they are not the same),
      # so it would be very nice to use the previous calculated values instead of calculating each one again.

      # Parenthesis (to explain why the mean of x is not a problem):
      # (y has already been normalized substracting its mean, xcorr is seminormalized because it needs to be divided by a factor related to the standar deviations of both subvector(x) and Y)
      # sx = subvector(X)
      # sy = &(Y[delay]) 
      # seminormalized_xcorr[delay] = sumatory( (sx[i] - mean(sx))*sy[i])
      #                               sumatory( sx[i]*y[i] - mean(sx) * sy[i]) 
      #                               sumatory( sx[i]*y[i]  ) - mean(sx) * sumatory(sy[i])
      #                               dotproduct(x,y) - mean(sx) * sumatory(sy[i])
      # \Parenthesis

      # Sumatories of x and y values only need to be calculated once (independent from the delay), then they are updated with the new values, this includes
      #   mean (sumatory divided by length)
      #   sum of cuadratic differences (standard deviation * number of elements)
      #     Explaination:
      #       mean(x) = sum[i=1,i=n](x[i]) / n
      #       sum_cuadratic_differences = sum[i=1,n]((x[i]-mean(x))**2)
      #       sum_cuadratic_differences = sum[i=1,n](x[i]**2 + mean(x)**2 - 2*x[i]*mean(x))
      #       sum_cuadratic_differences = n * mean(x)**2  +  sum[i=1,n](x[i]**2) - 2*mean(x)*(n*mean(x))
      #       sum_cuadratic_differences = sum[i=1,n](x[i]**2) - n * mean(x)**2
      # Also, the dot product of vectors when increasing their lengths are included in this group (the delays mess this up)

      # At this point, a choice needs to be made: use memory to store the values or move the delay loop as the outermost loop, my choice is the second one.

      # My method:
      #             Y: [----------------------------------------------------------------------------------------]
      #             X: [----------------------------------------------------------------------------------------]
      # Subvector1(X): [................................--------------------------------........................]
      # Subvector2(X): [................................--------------------------------N.......................]

      # For each delay, we delete that amount of values from the front of Y (in c, this can be done in O(1) with Y=Y+delay, store the original Y to free later if needed)
      # The it turns out in: for each X subvector, dotmultiply it with the front of the Y vector
      # XCorr with subvector1:
      #   Xcorr[i]=
      #                     Subvector1(X): [................................OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO........................]
      #                     Y :            [deleted_deleted_deleted_deleted_OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO------------------------]
      #   Xcorr[i+1]=   
      #                     Subvector1(X): [................................11111111111111111111111111111111........................]
      #                     Y :           [deleted_deleted_deleted_deleted_d11111111111111111111111111111111-----------------------]
      #   Xcorr[i+2]=   
      #                     Subvector1(X): [................................22222222222222222222222222222222........................]
      #                     Y :          [deleted_deleted_deleted_deleted_de22222222222222222222222222222222----------------------]
      #   ...
      # XCorr with subvector2:
      #   Xcorr[i]=   
      #                     Subvector2(X): [................................OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOON.......................]
      #                     Y :            [deleted_deleted_deleted_deleted_OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOON-----------------------]
      #   Xcorr[i+1]=   
      #                     Subvector2(X): [................................11111111111111111111111111111111N.......................]
      #                     Y :           [deleted_deleted_deleted_deleted_d11111111111111111111111111111111N----------------------]
      #   Xcorr[i+2]=   
      #                     Subvector2(X): [................................22222222222222222222222222222222N.......................]
      #                     Y :          [deleted_deleted_deleted_deleted_de22222222222222222222222222222222N---------------------]
      #   ...   
      
      # Notice how simple is to add another values to change from Xcorr[n] with subvector1 to the xcorr[n] with subvector2 (N values)  

      # Obviously, it's arranged in a different way: having the subvector_length as the innermost loop:

      # delay i
      #   start j
      #     length l
      #                       Subvector1(X): [................................OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO........................]
      #                       Y :            [deleted_deleted_deleted_deleted_OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO------------------------]

      #     length l+1 
      #                       Subvector2(X): [................................OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOON.......................]
      #                       Y :            [deleted_deleted_deleted_deleted_OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOON-----------------------]
      #     ...
      #   ...
      # delay i+1
      #   start j
      #     length l
      #                       Subvector1(X): [................................11111111111111111111111111111111........................]
      #                       Y :           [deleted_deleted_deleted_deleted_d11111111111111111111111111111111-----------------------]
      #     length l+1
      #                       Subvector2(X): [................................11111111111111111111111111111111N.......................]
      #                       Y :           [deleted_deleted_deleted_deleted_d11111111111111111111111111111111N----------------------]
      #     ...
      #   ...
      # ...
      # delay and subvector_start loops can be switched, they mean the same for the different vectors.

      # Another method could be to 'slide' one vector over the other, and then evaluate all the subvectors of the intersection between both vectors
      # In c this other method spends roughly the same time (the complexity is still O(N^3)), but in R it would need a third subvector (xy product elementwise), with the following overhead
      #############################################################################################################################

      ## PART2PARTCC ##############################################################################################################
      xsize = length(x);
      ysize = length(y);
      lmin = as.integer(max(min(lmin,xsize,ysize)),2); # lmin must be > 1 (subvectors of length 1 have variance == 0)

      y = y - mean(y);
      y_sumcuadraticdiff = sum(y^2);

      maxscore<-function(delays,thisscore){
        res = list(cmax=thisscore,pos1=0,pos2=0,lngth=0)

        for (delay in delays){
          for(subvector_start in 0:(xsize-lmin)){
            # All the subvectors lengths are checked at the same time (thanks to R's vectorized functions)
            xy_maxsize = min(xsize-subvector_start,ysize-delay); 
            
            ## x and y subvectors #########################
            # Read only vectors (in c: O(1), working with pointers)
            xx = x[(subvector_start+1) : (subvector_start+xy_maxsize)];
            yy = y[(          delay+1) : (          delay+xy_maxsize)];
            ###############################################
            
            ## Cumulative values ##########################
            # Related between iterations (in c it's not needed to calc them every time, in r it's better to recalc them due to the overhead)
            xydot_subvector = cumsum(xx*yy);
            x_sum_subvector = cumsum(xx);
            x_sumcuadratic_subvector = cumsum(xx^2);
            y_sum_subvector = cumsum(yy);
            ###############################################
            
            subvector_length = 1:xy_maxsize; #lengths
            x_mean = x_sum_subvector / subvector_length; #for each length, get the average of x's subvector
            denom = sqrt((x_sumcuadratic_subvector - x_mean * x_sum_subvector)*y_sumcuadraticdiff); #for each length, get the denom
            scorecand = (xydot_subvector - y_sum_subvector * x_mean)/denom;
            scorecalc = scorecand*subvector_length;
            newscore = max(thisscore,scorecalc[lmin:xy_maxsize],na.rm = TRUE); #remove the first lmin-1 elements from max (accounts for ~5.5% of the time), and remove the NaN (denom == 0,required)
            ## thisscore is included to avoid empty vector warnings (from scorecand after na are removed)

            if (newscore>thisscore){
              thisscore = newscore
              lngth = tail(which(scorecalc==newscore),1) # Longer subvector lengths are better
              res$lngth= lngth
              res$cmax = scorecand[lngth]
              res$pos1 = subvector_start +1; # position of the first element of subvector inside vector
              res$pos2 = delay +1; # position of the first element of subvector inside vector
            }
          } 
        }
        return (res);
      }
      #############################################################################################################################

      ## GO PARALLEL! #############################################################################################################
      if (nprocesses!=1){
      
        if(nprocesses < 1){ ## I have the power to schedule!
          cores = parallel::detectCores();
          if(is.na(cores))
            nprocesses = 1
          else{
            nprocesses = max(min(cores-1,floor(cores*0.8)),1); # Use 80% of the machine, leaving at least 1 core free and using at least 1 core
            # min is redundant because with floor() and any factor below 1.0, at least a processor will be always free.
            # max is redundant because if min(X)<1, the flow will go through the sequential route. 
          
            nprocesses = min(nprocesses,ysize-lmin+1) # avoid spawning more processes than delays
          }
        }
      }

      if (nprocesses>1){
        # SCHEDULE JOBS
        # The delays are scheduled in a round-robin fashion among the processes.
        # RR was chosen because its very simple to implement and the number of operations required for the evaluation of each delay
        # stays the same or decreases when the delay increases
        # A function to calculate the cost (number of iterations required for each evaluation) is provided in case a custom scheduling is required in the future
        ###############################################################
        evaleffortbydelay<-function(delays,x_size,y_size,lmin){
          # iterations per delay and start
          # delays\starts
          #       0 1 2 3 4 5 6 <- start
          #   0   6 5 4 3 2 1 0
          #   1   6 5 4 3 2 1 0
          #   2   5 5 4 3 2 1 0 <- x&y subvector_length (number of elements to compute)
          #   3   4 4 4 3 2 1 0
          #   4   3 3 3 3 2 1 0
          #   5   2 2 2 2 2 1 0
          #   6   1 1 1 1 1 1 0
          #   7   0 0 0 0 0 0 0
          # lmin is the number of columns and rows to delete

          # estimates the effort per delay (assuming lmin=0)
          # a delay performs (n-k+1)k+(k-1)k/2 operations if k<=n
          # where n=x_size and k=y_size-delay
          # if k>n a delay performs (n-k+1)k+(k-1)k/2
          # where n=x_size and k=n
          ###############
          # so: a delay performs (n-k+1)k+(k-1)k/2 operations
          # where n=x_size and k=min(y_size-delay,n)
          #:: effort = k*(n+(1-k)/2)
          k = apply(as.matrix(delays),1,function(delay) return(min(y_size-delay,x_size)))
          tot = k*(x_size+(1-k)/2)
          
          if (lmin!=0)
            return (tot-evaleffortbydelay(delays,lmin-1,y_size,0)) ##remove last columns (lmin columns)

          return(tot)
        }
        evalallefforts<-function(x_size,y_size,lmin){
          # Returns a vector with the effort per delay (delays = 0:(ysize-lmin), so effort(delay) = vector[delay+1]  )
          if (x_size==0)
            return(seq(from=0,to=0,length.out=y_size+1))
        
          sums = cumsum(seq(from=x_size,to=1))
          lensums = length(sums)

          if(lensums>y_size){
            lensums = y_size
            length(sums) = y_size
            vec = sums
          }else{
            last = sums[lensums]
            vec = seq(from=last,to=last,length.out=y_size)
            vec[1:lensums] = sums
          }
          res = c(rev(vec),0)

          if (lmin!=0){
            res=res-evalallefforts(lmin-1,y_size,0) ##remove last columns (lmin columns)
            length(res) = y_size -lmin+1 ##remove last rows (lmin rows)
          }
          
          return(res)
        }
        ###############################################################

        paralleleval<-function(id){
          if(ysize-lmin<id) #not enough work for this process
            return(NULL) 

          return(maxscore(seq(from=id,to=ysize-lmin,by=nprocesses),score))
        } # The seq is created in the thread
        jobs = lapply(1:(nprocesses-1), function(children) { # Schedule work on other processes 
          parallel::mcparallel(paralleleval(children))
        });
        thisres = paralleleval(0) # Schedule work on this process
        
        lapply(parallel::mccollect(jobs),function(res){ ##Get the best of the results (could use reduction... but is just the number of processes, shouldn't be very big)
          if ((!is.null(res))&&(res$cmax * res$lngth > thisres$cmax * thisres$lngth))
            thisres <<- res
        })
     
        return (thisres);
      }
      #############################################################################################################################

      return (maxscore(0:(ysize-lmin),score)); # Check all the delays
    }

    #
    # Start the analysis by determining who is aguja and pajar
    aguja<-x
    pajar<-y
    nominal<-TRUE
    if (length(x) > length(y)) {
      nominal<-FALSE
      aguja<-y
      pajar<-x
    }

    res = fastpart2partcc(pajar,aguja,lmin,nprocesses=nprocesses)
    #res$pos1 is for pajar, and res$pos2 is for aguja

    return(list(orderok=nominal,pos1=res$pos2, pos2=res$pos1,cmax=res$cmax, lngth=res$lngth))
    #return$pos1 is for aguja and return$pos2 is for pajar
  }
  #
  intcurv<-function(x) {
    gg<-function(l,y){
      return(sum(rollmean(y[1:l],2)))
    }
    d1<-diff(x,1)
    d3<-c(d1[1],d1)
    d31<-diff(d3,1)
    crv<-c(d31,d31[length(d31)])
    abcrv<-abs(crv)
    cvint<-c(abcrv[1],apply(as.matrix(2:length(x)),
                            1,gg,abcrv)+abcrv[1])
    return(list(crv=crv,crv_abs=abcrv,crv_int=cvint))
  }
  #
  interpolate<-function(x,paso=10,eps=0.0001) {
    gg<-function(x,eps){
      vx<-x
      siz<-length(x)
      for (i in 2:siz){
       if (( vx[i] - vx[i-1] ) < eps ) {
          vx[i] <- vx[i-1] + eps
        }
      }
      return(vx)
    }
    siz<-length(x)
    res<-intcurv(x)
    rsvint<-res[["crv_int"]]
    rsv   <-res[["crv"]]
    vx<-gg(rsvint,eps)
    vy<-1:siz
    xinterp<-linspace(vx[1],vx[siz],floor((vx[siz]-vx[1])/paso))
    intk<-xinterp
    ints<-approx(x=vx,y=vy,xout=xinterp,method="linear")
    crvk<-approx(x=vy,y=rsv,xout=ints$y,method="linear")
    return(list(k=crvk$y,s=ints$y,intk=intk))
  }
  #
  #
  if (length(mm) == 0 ) {
    rs<-list(sc=0, pos1=0, pos2=0, lngth=0)
    return(rs)
  }
	imm<-interpolate(mm,paso=stp)
	im2<-interpolate(m2,paso=stp)
	res<-part2part(imm$k,im2$k,floor(0.75*min(length(imm$k),length(im2$k))),nprocesses=nprocesses)
	if ( res$pos1 == 0 & res$pos2 == 0 ) {
		return(list(sc=0,pos1=0,pos2=0,lngth=0))
	}
  if ( res$orderok ) {
		lg <- floor(imm$s[res$lngth+res$pos1-1]-imm$s[res$pos1])
   	rs <- list(sc=res$cmax,pos1=floor(imm$s[res$pos1]),pos2=floor(im2$s[res$pos2]),lngth=lg)
	} else {
		lg <- floor(im2$s[res$lngth+res$pos1-1]-im2$s[res$pos1])
   	rs <- list(sc=res$cmax,pos1=floor(im2$s[res$pos1]),pos2=floor(imm$s[res$pos2]),lngth=lg)
	}
	return(rs)
}
