USsenateClass <- function(x, senate=readUSsenate(),
   Office='Office', state='state',
   surname='surname', district='district',
                          senatePattern='^Senate') {
##
## 1.  subset x houseSenate
##
  hS <- which(names(x) == Office)
  if(length(hS) != 1)
      stop('failed to find a unique column of x matching ',
           Office)
#
  xSen <- grep(senatePattern, x[, hS])
  x.Sen <- x[xSen,]
  nx <- nrow(x.Sen)
  incumbent <- rep(FALSE, nx)
##
## 2.  Key
##
#  keyx <- with(x.Sen, paste(state, surname, sep=":"))
  keyx <- paste(x.Sen[, state], x.Sen[, surname], sep=':')
  keysen <- with(senate, paste(state, surname, sep=":"))
##
## 3.  find
##
  used <- (keysen %in% keyx)
  tabs <- table(keysen[used])
  if(any(tabs>1)){
#   State with both senators having the same last name
      warning("USsenateClass function may fail for a state where ",
              " both senators have the same last name:",
              names(tabs)[tabs>1][1])
  }
  senu <- senate[used,]
  keyu <- keysen[used]
  rownames(senu) <- keyu
#
  found <- (keyx %in% keyu)
  incumbent[found] <- TRUE
  Class <- senu[keyx[found], 'Class']
  names(Class) <- keyx[found]
#
  out <- character(nx)
  names(out) <- keyx
  out[keyx[found]] <- Class
# out2 <- character(nx)
#  out2[found] <- Class
#
  if(any(!found)){
      oops <- which(!found)
      senun <- senate[!used, ]
      for(o in oops){
          sto <- as.character(x.Sen[o, state])
#
          clo <- senun[senun$state==sto, 'Class']
          if(length(clo)>0){
              out[o] <- paste(sort(clo), collapse=' or ')
          } else {
#         Not in senun;  look in senu
              clo. <- senu[senu$state==sto, 'Class']
              if(length(clo.)>0){
                  out[o] <- paste(sort(clo.), collapse=' or ')
              } else {
                  print(x.Sen[o,])
                  warning('Could not find senator from state ',
                          sto)
              }
          }
      }
  }
##
## 4.  prepare for output
##
  Dist <- as.character(x[, district])
  Dist[xSen] <- out
  inc <- rep(NA, nrow(x))
  inc[xSen] <- incumbent
#
  Out <- data.frame(incumbent=inc, district=Dist,
                    stringsAsFactors=FALSE)
  Out
}

