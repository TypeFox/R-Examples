
markbedrest =  function(dset, age, loc, TS, cts, rstdr, rstfl, per, TH, wd, nvm, nmin, nmax, EC)
{
#  library("lubridate")
#  library("chron")
  vector.is.empty <- function(x) return(length(x) ==0 )
  
  
    if(!missing(age))
    {
    if(age=="Adult" | age=="ADULT" | age=="ad" | age=="AD"){
      age="adult"}
    }
    if(!missing(age))
    {
      if(age=="Youth" | age=="YOUTH" | age=="yo" | age=="YO"){
        age="youth"}
    }
    if(!missing(age))
    {
      if(age=="Preschool" | age=="PreSchool" |age=="PRESCHOOL" | age=="ps" | age=="PS" 
         | age=="Preschooler" | age=="PreSchooler" |age=="PRESCHOOLER"){
        age="preschool"}
    }
    
    if(!missing(loc))
    {
      if(loc=="Wrist" | loc=="WRIST" | loc=="wr" | loc=="WR" | loc=="Wr"){
        loc="wrist"}
    }
    if(!missing(loc))
    {
      if(loc=="Waist" | loc=="WAIST" | loc=="wa" | loc=="WA" | loc=="Wa"){
        loc="waist"}
    }
    
    if(missing(TS))  {TS ="TS"}  
    ts = dset[,TS]

    
# ERRORS that prevent it from running  

  if(missing(dset)){
    stop("You must specify the dset")}

# arg approach
  else if(missing(age) & missing(loc) & missing(TH)) {stop("Unless specifying age and loc, you must specify a TH")}
  else if(missing(age) & missing(loc) & missing(wd)) {stop("Unless specifying age and loc, you must specify a wd")}
  else if(missing(age) & missing(loc) & missing(nvm)) {stop("Unless specifying age and loc, you must specify a nvm")}

# age loc approach
  else if(!missing(age) & missing(loc)){
    stop("If you specify age, you must specify loc, where the monitor was worn")}
  else if(!missing(loc) & missing(age)){
    stop("If you specify loc, you must specify age (of wearer)")}
  


########################################
  else{
    if(missing(EC)) {EC=TRUE}  
    if(EC==TRUE){
        if(!missing(age)) {
          if(!(age=="adult" | age=="youth" | age=="preschool")){
            stop(paste0(paste0("'", age, "'"), ' is not a valid entry for age.  Only "preschool", "youth" and "adult" are valid entries'))
            }
#          if(age=="preschool" & !(as.POSIXct(dset$TS[2]) == as.POSIXct(dset$TS[1]) + dseconds(15))){
          if(age=="preschool" & !(dset$TS[2] == dset$TS[1] + dseconds(15))){
              stop('preschool data must be in 15 second epochs')
            }
#          if(age!="preschool" & !(as.POSIXct(dset$TS[2]) == as.POSIXct(dset$TS[1]) + dseconds(60))){
          if(age!="preschool" & !(dset$TS[2] == dset$TS[1] + dseconds(60))){
            stop('non-preschool data must be in 60 second epochs')
            }
          }
          
        if(!missing(age) & !missing(loc)){
          if(age=="preschool" & loc=="wrist"){
          stop("Sorry, we don't have validated cut-offs for that age loc combination")
          }
        
          if(!(loc=="waist" | loc=="wrist")){
          stop(paste0(paste0("'", loc, "'"), ' is not a valid entry for loc.  Only "waist" and "wrist" are valid entries'))
          }
        }
    }
    
# WARNINGS, but it can still run

    if(!(missing(age) & missing(loc))){ 
        if(!missing(TH)) 
        {
        print("WARNING! age and loc were specified thereby implying a TH, so TH specified will be ignored")
        }
        if(!missing(wd)){
        print("WARNING! age and loc were specified thereby implying a wd, so wd specified will be ignored")
        }
        if(!missing(nvm)){
        print("WARNING! age and loc were specified thereby implying a nvm, so nvm specified will be ignored")
        }
    }

#####################################
    
    if(missing(TS))  {TS ="TS"}
    if(missing(cts)) {cts = "Axis1"}
    
    if(missing(rstdr))  {rstdr = getwd()}
    if(missing(rstfl)) {rstfl = "subj_slp_sum"}
    
#    if(missing(TH))   {TH = 250}
#    if(missing(wd))   {wd = 3000}
    
    if(missing(nmin)) {nmin = 60}
    if(missing(nmax)) {nmax = 60}
#    if(missing(nvm))  {nvm = 50}
    if(missing(per))  {per  = 60}

if(!missing(age)){
      if(missing(nmin)) {if (age=="preschool") {nmin = 240} }
      if(missing(nmax)) {if (age=="preschool") {nmax = 240} }
      if(missing(per))  {if (age=="preschool") {per  = 240} }
      }

###################################    
    if(!missing(age) & !missing(loc))
      {
      if(age=="preschool" & loc=="waist"){
        TH=30
        wd=50
        nvm=10
        #_60
        }
        
        if(age=="youth" & loc=="waist"){
          TH=20
          wd=500
          nvm=50
          }
        if(age=="youth" & loc=="wrist"){
          TH=250
          wd=3000
          nvm=50
        }
        if(age=="adult" & loc=="waist"){
          TH=10
          wd=1000
          nvm=140
        }
        if(age=="adult" & loc=="wrist"){
          TH=400
          wd=1500
          nvm=150
          per=30
        }
      }
    # set length of window to be searched for waking (minutes) #
    lagwin = 2*per
    # set length of window to be searched for sleep (minutes) #
    lagslp = 2*per
        
    
    ts = as.vector(dset[,TS])
    vm = as.numeric(dset[,cts])
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
            if(falloff > l){falloff = l} #makes sure sleep vector doesn't get too long
            if(falloff <= lstchng)
            {falloff = lstchng
             sleep.end = sleep.end[1:length(sleep.end)-1]
             slp=1
             j=0
            }
            else {sleep = replace(sleep, lstchng:falloff, "a")
                  slp = 1
                  lstchng = falloff
                  if(lstchng > l){lstchng =l} #makes sure sleep vector doesn't get too long
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
              if(lstchng > l){lstchng =l} #makes sure sleep vector doesn't get too long
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
                if(lstchng > l){lstchng =l} #makes sure sleep vector doesn't get too long
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
    
    if(vector.is.empty(sleep.start)){sleep.start="."}
    if(vector.is.empty(sleep.end)){sleep.end="."}
    if(sleep.end[1] < sleep.start[1]){sleep.end= c(".", sleep.end)}
    if(length(sleep.end)>length(sleep.start)){sleep.start= c(sleep.start, ".")}
    if(length(sleep.end)<length(sleep.start)){sleep.end= c(sleep.end, ".")}  
    
    bedrest.start = sleep.start
    bedrest.end   = sleep.end
    
    subj_slp_sum = data.frame(bedrest.start, bedrest.end)
    result = paste(rstdr, "/", rstfl, ".csv", sep="")
    write.csv(subj_slp_sum, result, row.names = FALSE)  
    
    return = dset_slp
  } #end of else (not missing)
  
} #end of function markbedrest
