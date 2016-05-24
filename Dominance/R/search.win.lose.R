search.win.lose <-
function(data_sheet, ...)
{
    args = list(...)

if ("bits" %in% names(args))      
    bits <- args$bits
    else
    bits<- ""


if ("count_all" %in% names(args))      
    count_all <- args$count_all
    else
    count_all<- FALSE
 
     win_lose<-integer()
     max_items<- 0
	   max_items <- max(data_sheet$action.from,na.rm=TRUE)
	   count_max <- length(data_sheet$action.from)
	   if (max_items < max(data_sheet$action.to,na.rm=TRUE))
		    max_items <- max(data_sheet$action.to,na.rm=TRUE)
    if (max(data_sheet$item.number,na.rm=TRUE) < max_items)
	   {
		    print("Error max count of items does not match")
	      break;
      }
    max_items <-  max(data_sheet$item.number,na.rm=TRUE)       # because some items may not be used in actions

 # ------------------ checking whether bits is in a proper format
	if (bits != "0")
	{ 	
 	   max_actions <- max(data_sheet$action.number,na.rm=TRUE)

		detect_0 <- detect_bits(bits,FALSE)
		detect_1 <- detect_bits(bits,TRUE)
		length_bits <- 0
		if (detect_1[1] != -1) # -1: no bit 1 found
		length_bits <- length(detect_1)
		if (detect_0[1] != -1) # -1: no bit 1 found
		length_bits <- length_bits+length(detect_0)
		if (length_bits !=  max_actions )
			# Nr. of its == 1 + Nr. of bits = 0
		{
			print("Error: max count of actions does not length of bits")
		      break;
		}
      }

if ("weighting" %in% names(args))
    {
      if (length(args$weighting) != max(data_sheet$action.number,na.rm=TRUE))
       { 
        print(" length(weigting) !=  max(data_sheet$action.number)")
        break;
        } 
       weighting_internal <- args$weighting
     } else  
       weighting_internal <- data_sheet$weighting      
# -------- end checking whether bits is in a proper format

        temp_win <- integer()
        temp_lose <- integer()

            Data_temp = data.frame ("action.from"=data_sheet$action.from,"action.to"=data_sheet$action.to,"kind.of.action"=data_sheet$kind.of.action) #build an tem data frame to prevent missing dat in line 0 to max (actions)
            Data_temp = subset(Data_temp,data_sheet$action.from > 0 & data_sheet$action.to > 0 )  #delete all sctions from and to 0
            Data_temp = subset(Data_temp,match(Data_temp$kind.of.action,detect_1,nomatch=0)> 0 )  #delete all not with bytes set actions
            for (x in 1:max(data_sheet$action.number,na.rm=TRUE) )
            {  
           #  print(x)

             if (data_sheet$classification[x] ==2)
             {
              temp <- subset(Data_temp$action.from, x == Data_temp$kind.of.action & data_sheet$classification[x] ==2 )
              if ((length(temp) > 0))  
              {
#                print("A")
#                print(temp)
#                print(data$win)
                temp_win   <-  append(temp_win,temp)
              }
              temp <- subset(Data_temp$action.to,  x == Data_temp$kind.of.action & data_sheet$classification[x] ==2  )              
              if ((length(temp) > 0))  
              {
#                print("B")
#                print(temp)
#                print(data$lose)
                 temp_lose  <-  append(temp_lose,temp)
              }
             } 
             
             if (data_sheet$classification[x] ==1)
             {
              temp <- subset(Data_temp$action.from, x == Data_temp$kind.of.action & data_sheet$classification[x] ==1 )
              if ((length(temp) > 0))
              {
 #               print("C")
 #               print(temp)
 #               print(data$lose)
                 temp_lose  <-  append(temp_lose,temp)
                 rm(temp)
             }

              temp <- subset(Data_temp$action.to,  x == Data_temp$kind.of.action & data_sheet$classification[x] ==1  )              
              if ((length(temp) > 0))
              {
 #               print("D")
 #               print(temp)
 #               print(data$win)
                 temp_win   <-  append(temp_win,temp)
                 rm(temp)
              }  
             } 

             if (data_sheet$classification[x] > 2)
             break("error  Data_temp$classification > 2 not defined")
            } #  for (x in 1:max(Data_temp$action.number,na.rm=TRUE) )          
                                          
      win_lose <- data.frame("wins"=temp_lose,"loses"=temp_win)
      return(list(original=data_sheet,data.win.lose=win_lose,items=max_items))
}
