sls <-
function(l){
  ifelse(l %in% c("e","a","i","o","n","r","t","l","s","u"),1,
    ifelse(l %in% c("d","g"),2,
     ifelse(l %in% c("b","c","m","p"),3,
      ifelse(l %in% c("f","h","v","w","y"),4,
        ifelse(l %in% c("k"),5,
          ifelse(l %in% c("j","x"),8,
            ifelse(l %in% c("q","z"),10,0)))))))
}