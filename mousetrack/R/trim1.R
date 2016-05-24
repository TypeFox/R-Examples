## another trimming function

.packageName <- 'mousetrack'

trim1 <-function(vect1, vect2, thresh){


    o1 = tail(vect1, 1) # get final value of trajectory
    o2 = tail(vect2, 1) 

    flipo1 = vect1[length(vect1):1] # take original vector, and reverse
    flipo2 = vect2[length(vect2):1]

# how far (in Euclidean distance, pixels) is x,y (reversed) from the final click at each time step
    dists = sqrt((flipo1-o1)^2 + (flipo2-o2)^2)

    latindx = which(dists > thresh)[1] # where is the x,y (reversed) trajectory THRESH pixels (aprox) from the final click
    latindx2 = length(vect1)-latindx # where is the x,y (non-reversed) trajectory TRHESH pixels from the final click

    trimmed1 = vect1[latindx2:length(vect1)] # based on latindex2, grab the x,y (non-reversed) trajectory for the last THRESH pixels. 
    trimmed2 = vect2[latindx2:length(vect2)]

    return( list(latindx = latindx2, trimmed1 = trimmed1, trimmed2 = trimmed2) )

}
