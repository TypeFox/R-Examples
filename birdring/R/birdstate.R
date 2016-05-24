
# function to create a variable that indicates whether a bird was alive or dead
# at the time of reencounter
birdstate <- function(x){
  # x     condition (numeric code) as given by EURING
  livemat<-data.frame(condition=0:9, 
                      livestate=c("unknown", "dead", "dead", "dead", "sick", "sick", "alive", "alive", "alive", "alive"))
  state <- livemat$livestate[match(x, livemat$condition)]
  return(state)
}
