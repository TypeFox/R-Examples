match.list <- function (template.list,user.list)
{
  if ((length(template.list)==0)||(length(user.list)==0)) 
    return (template.list)
  else
  {
    attr1 <- attributes(template.list)[[1]]
    attr2 <- attributes(user.list)[[1]] # if there is at least one attribute in the list, then attr2 
                                        # will have the size of user.list, otherwise it will have size 0.
  
    if (length(attr2)==0) attr2 <- rep("", length(user.list)) 
    
    matches  <- match(attr2,attr1) # matches will have size of attr2 and user.list > 0 by now.
    
    copylist <- matches[!is.na(matches)]
    
    if (length(copylist)>0)
    {
      template.list[attr1[copylist]] <- user.list[attr1[copylist]]      
    }      
    
    guesslist <- which(attr2=="")
    
#   guesslist <- (1:length(matches))[is.na(matches)] # returnes indices of elements containing TRUE
        
    if (length(guesslist)>0)
    {
      misslist   <- setdiff(1:length(template.list), matches)
      
      if (length(misslist)>0)
      {
        # user.list can have more missed elements then 
        # the number of missed elements in the template.list
        copylength <- min(length(guesslist),length(misslist)) 
        template.list[misslist[1:copylength]] <- user.list[guesslist[1:copylength]]          
      }
    }
    return (template.list)    
  }
}