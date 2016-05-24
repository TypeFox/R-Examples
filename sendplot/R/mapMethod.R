
#
# finds bounding region based on given information
#   current methods for retrieve are mode or median 
#

mapMethod <- function(automap.method, temp){


  # determin column locations of where tifs differ
  row.count = rowSums(temp, na.rm=TRUE)
  col.loc = which(row.count>0)

    
  # store largest break between locations to seperate lower and upper bound
  len = length((col.loc[1]):(col.loc[2]))
  brk = 1
  for(i in 2:(length(col.loc)-1)){
    l = length((col.loc[i]):(col.loc[i+1]))
    if(l > len){
      brk = i + 1
      len = l
    }
  }
  # start regions (columns) of potential lower and upper bound
  s.reg.1 = col.loc[1]
  s.reg.2 = col.loc[brk]

  # full region of upper and lower bound (columns)
  col.reg.1 = (col.loc[1]):(col.loc[(brk-1)])
  col.reg.2 = (col.loc[brk]):(col.loc[(length(col.loc))])


 
  
  # if find by median
  if(automap.method == "median"){
      #
      # find upper bound
      #
     
      # column loc is median 
      r1.col = floor(median(col.reg.1))
      # row location is length of different in determine column
      r1.row = which(temp[r1.col,]>0)[1]
      add.r = floor((row.count[r1.col])/2) 
      up.left.row = r1.row + add.r
      # add additional to compensate for row
      add.c = floor((length(which(row.count[col.reg.1] == row.count[r1.col])))/2)
      up.left.col = r1.col + add.c 
  
      #
      # find lower bound
      #
        
      # column loc is median
      r2.col = floor(median(col.reg.2))
      # row location is length of different in determine column
      r2.row = which(temp[r2.col,]>0)[1]
      add2.r = floor((row.count[r2.col])/2) 
      low.right.row = r2.row + add2.r
      # add additional to compensate for row
      add2.c = floor((length(which(row.count[col.reg.2] == row.count[r2.col])))/2)
      low.right.col = r2.col + add2.c 

    } # end if median

    
   # if find by mode (default)
    if(automap.method == "mode"){
  
      # 
      # find upper bound
      #

      # finds which length repeats most often
      vec = row.count[col.reg.1] 
      u = unique(vec)
      fr = u[1]
      len = length(which(vec == fr))
      for(i in u[-1]){
        if(i != 0){
          c = length(which(vec == i))
          if(c > len){
            fr = i
            len = c
          }
        }        
      }

      # column is which repeats most often 
      r1.col = which(row.count == fr)[1]
      # row location is length of different in determine column
      r1.row = which(temp[r1.col,]>0)[1]
      add.r = floor(fr/2) 
      up.left.row = r1.row + add.r
      # add additional to compensate for row
      add.c = floor((len/2))
      up.left.col = r1.col + add.c
      
      # if no row repeats, uses median method      
      if(len == 1){
        #dx1 = floor(length(col.reg.1)/2)
        #r1.col = col.reg.1[dx1]

        # column loc is median
        r1.col = floor(median(col.reg.1))

        # make sure median is not one of the places that was 0 
        if(!is.na(match(r1.col, col.reg.1[which(vec == 0)]))){
          sDx = which(col.reg.1 == r1.col)
          upDx = sDx
          unF = TRUE
          iM = 0
          while(unF & (upDx < length(col.reg.1))) {
            iM = iM + 1
            unF = !is.na(match(col.reg.1[sDx+iM], col.reg.1[which(vec == 0)]))
            upDx = upDx+1
          }
          lwDx = sDx
          lnF = TRUE
          iM = 0
          while(lnF & (lwDx > 0)) {
            iM = iM + 1
            lnF = !is.na(match(col.reg.1[sDx-iM], col.reg.1[which(vec == 0)]))
            lwDx = lwDx - 1
          }
          if( (upDx-sDx) <= (sDx-lwDx) ) {
            r1.col= col.reg.1[upDx]
          }else{
            r1.col= col.reg.1[lwDx]
          }  
        }        
        # row location is length of different in determine column
        r1.row = which(temp[r1.col,]>0)[1]
        add.r = floor((row.count[r1.col])/2) 
        up.left.row = r1.row + add.r
        # add additional to compensate for row
        add.c = floor((length(which(row.count[col.reg.1] == row.count[r1.col])))/2)
        up.left.col = r1.col + add.c 
      }

         
      #
      # lower bound
      #
      
      # finds which length repeats most often
      vec2 = row.count[col.reg.2] 
      u2 = unique(vec2)
      fr2 = u2[1]
      len2 = length(which(vec2 == fr2))
      for(i2 in u2[-1]){
        if(i2 != 0){
          c2 = length(which(vec2 == i2))
          if(c2 > len2){
            fr2 = i2
            len2 = c2
          }
        }        
      }
      
      # column is which repeats most often 
      r2.col =  which(row.count == fr2)
      r2.col = r2.col[which(r2.col > (s.reg.1 + length(col.reg.1)))][1]
      # row location is length of different in determine column
      r2.row = which(temp[r2.col,]>0)[1]
      add2.r = floor(fr2/2) 
      low.right.row = r2.row + add2.r
      # add additional to compensate for row
      add2.c = floor((len2/2))
      low.right.col = r2.col + add2.c 

      # if no row repeats, uses median method  
      if(len2 == 1){
        
        # column loc is median
        r2.col = floor(median(col.reg.2))
        # make sure median is not one of the places that was 0 
        if(!is.na(match(r2.col, col.reg.2[which(vec2 == 0)]))){
          sDx = which(col.reg.2 == r2.col)
          upDx = sDx
          unF = TRUE
          iM = 0
          while(unF & (upDx < length(col.reg.2))) {
            iM = iM + 1
            unF = !is.na(match(col.reg.2[sDx+iM], col.reg.2[which(vec2 == 0)]))
            upDx = upDx+1
          }
          lwDx = sDx
          lnF = TRUE
          iM = 0
          while(lnF & (lwDx > 0)) {
            iM = iM + 1
            lnF = !is.na(match(col.reg.2[sDx-iM], col.reg.2[which(vec2 == 0)]))
            lwDx = lwDx - 1
          }
          if( (upDx-sDx) <= (sDx-lwDx) ) {
            r2.col= col.reg.2[upDx]
          }else{
            r2.col= col.reg.2[lwDx]
          }  
        }
        
        # row location is length of different in determine column
        r2.row = which(temp[r2.col,]>0)[1]
        add2.r = floor((row.count[r2.col])/2) 
        low.right.row = r2.row + add2.r
        # add additional to compensate for row
        add2.c = floor((length(which(row.count[col.reg.2] == row.count[r2.col])))/2)
        low.right.col = r2.col + add2.c 
      }
    } # end if mode

    # make bounding object
    bounds = list()
    bounds$up.left = c(up.left.row,up.left.col)
    bounds$low.right = c(low.right.row,low.right.col)

    # return bounding region 
    return(bounds)
    
  }
