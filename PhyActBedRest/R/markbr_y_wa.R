# v added argument for length of initial averages
# 3 adds line to make sure sleep column doesn't get longer than initial matrix
# a adds code to check sleep period is checked against nmin & nmax if end of vmwin is reached 

markbr_y_wa =  function(dset, TS, col, rstdr, rstfl, TH, wd, nmin, nmax, nvm)
{    
  if (missing(dset)){
    return("You must specify the dset")}
  
  else{
    
    if(missing(TS))  {TS ="TS"}
    if(missing(col)) {col = "Axis1"}
    
    if(missing(rstdr))  {rstdr = getwd()}
    if(missing(rstfl)) {rstfl = "subj_slp_sum"}
    
    if(missing(TH))        {TH = 20}
    if(missing(wd)) {wd = 500}
    
    if(missing(nmin)) {nmin = 60}
    if(missing(nmax)) {nmax = 60}
    if(missing(nvm))  {nvm = 50}
    #if(missing(per))  {per = 60}
    per = 60
    
    # set length of window to be searched for waking (minutes) #
    lagwin = 2*per
    # set length of window to be searched for sleep (minutes) #
    lagslp = 2*per
        
    
    ts = as.vector(dset[,TS])
    vm = as.numeric(dset[,col])
    l = length(ts)
    rd_l =trunc(l/per, 0) 
    rm_l = l - (per * rd_l)
    lstchng = 1
    sleep = rep("",l)
    slp = 0 #awake
    sleep.end = NULL
    sleep.start = NULL
    H = rd_l
    
    #Calculate the [hour] average counts
    h_avg = rep (NA, rd_l +1)
    for (g in 1:rd_l)
    {h_avg[g] = mean(vm[(1+((g-1)*per)):(g*per)])}
    
    if(rm_l>0){
      h_avg[rd_l +1] = mean(vm[(1+(rd_l*per)):l])
      H = rd_l +1}
    
    if (h_avg[1] < TH) 
    {slp = 1
     sleep.start = ts[1]} 
    
    for (m in 2:H)
    { 
      if(slp == 0) {
        if(h_avg[m-1] >= TH & h_avg[m] < TH )
        {i = (m * per) 
         begwin = i - lagslp +1
         if(i > l) {i = l}
         if(begwin < 1){begwin =1}
         winvm = vm[begwin:i]
         j=length(winvm)
         while (j>1){
           if(winvm[j] > nvm & winvm[j-1] > nvm)
           {falloff = begwin + j
            if(falloff > l){falloff = l} #makes sure sleep vector doesn't too long
            if(falloff <= lstchng)
            {falloff = lstchng
             sleep.end = sleep.end[1:length(sleep.end)-1]
             slp=1
             j=0
            }
            else {sleep = replace(sleep, lstchng:falloff, "a")
                  slp = 1
                  lstchng = falloff
                  if(lstchng > l){lstchng =l} #makes sure sleep vector doesn't too long
                  lc =  substr(ts[lstchng],1,19)
                  sleep.start = c (sleep.start, lc)
                  j = 0
            }
           }
           j = j - 1
         }#end of "while"
         
         if (j==1) 
         {falloff = i - lagslp + j
          if(falloff < 1){falloff =1}
          sleep = replace(sleep, lstchng:falloff, "a")
          slp = 1
          lstchng = falloff
          if(lstchng > l){lstchng =l} #makes sure sleep vector doesn't too long
          lc =  substr(ts[lstchng],1,19)
          sleep.start = c (sleep.start, lc)
         }
        }
      }
      
      if(slp == 1) {
        if(h_avg[m-1] < TH & h_avg[m] > TH )
        {  i = (m * per)
           begwin = i - lagwin +1
           if(begwin<1) {begwin=1}
           if(i > l) {i = l}
           winvm = vm[begwin:i]
           lw=length(winvm)
           
           k=2	
           while(k < lw)
             
           {if(winvm[k-1] > wd & winvm[k] > wd)
           {wake = begwin + k - 1
            if(wake-lstchng < nmin)
            {
              sleep = replace(sleep, lstchng:wake, "a")
              lstchng=wake
              if(lstchng > l){lstchng =l} #makes sure sleep vector doesn't too long
              sleep.start = sleep.start[1:length(sleep.start)-1]
            }
            else
            {
              if(wake-lstchng >= nmax) 
              {
                sleep = replace(sleep, lstchng:wake, "br")
                lstchng=wake
                lc =  substr(ts[lstchng],1,19)
                sleep.end = c (sleep.end, lc)
              }
              else	 
              {
                sleep = replace(sleep, lstchng:wake, "n")
                sleep.start = sleep.start[1:length(sleep.start)-1]
                lc =  substr(ts[lstchng],1,19)
                lstchng=wake
                if(lstchng > l){lstchng =l} #makes sure sleep vector doesn't too long
                lc =  substr(ts[lstchng],1,19)
              }
            }
            slp = 0
            k=lw + 1
           } # end of "if(winvm[k-1] > wd & winvm[k] > wd)"
            else{k=k+1}
           }#end of "while"
           if(k == lw) # if we have reached the end of the window without finding two consecutive epochs > wd 
           {
             wake = begwin + k - 1							
             if(wake-lstchng < nmin)
             {
               sleep = replace(sleep, lstchng:wake, "a")
               lstchng=wake
               if(lstchng > l){lstchng =l} #makes sure sleep vector doesn't too long
               sleep.start = sleep.start[1:length(sleep.start)-1]
             }
             else
             {
               if(wake-lstchng >= nmax) 
               {
                 sleep = replace(sleep, lstchng:wake, "br")
                 lstchng=wake
                 lc =  substr(ts[lstchng],1,19)
                 sleep.end = c (sleep.end, lc)
               }
               else   
               {
                 sleep = replace(sleep, lstchng:wake, "n")
                 sleep.start = sleep.start[1:length(sleep.start)-1]
                 lc =  substr(ts[lstchng],1,19)
                 lstchng=wake
                 if(lstchng > l){lstchng =l} #makes sure sleep vector doesn't too long
                 lc =  substr(ts[lstchng],1,19)
               }
             }
             slp = 0
           }
        }#end of "if(h_avg[m-1] < TH & h_avg[m] > TH )"
      }#end of "if(slp == 1)"
      
      if(slp == 0) {sleep = replace(sleep, lstchng:l, "a")} else {sleep = replace(sleep, lstchng:l, "br")}
    }# end of 'm' (hour) loop.
    bedrest = sleep
    
    dset_slp = data.frame(dset, bedrest)	
    
    # if(sleep.end[1] < sleep.start[1]){sleep.end= c(".", sleep.end)}
    if(length(sleep.end)>length(sleep.start)){sleep.start= c(sleep.start, ".")}
    if(length(sleep.end)<length(sleep.start)){sleep.end= c(sleep.end, ".")}  
    
    bedrest.start = sleep.start
    bedrest.end   = sleep.end
    
    subj_slp_sum = data.frame(bedrest.start, bedrest.end)
    result = paste(rstdr, "/", rstfl, ".csv", sep="")
    write.csv(subj_slp_sum, result, row.names = FALSE)  
    
    return = dset_slp
  } #end of else (not missing)
  
} #end of function markbr_y_wa
