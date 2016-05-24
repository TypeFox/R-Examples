create.Lexis.matrix <-
function(HMD.dataset, Sex="Female", minage=50,
         maxage=100, minyear=1950, maxyear=2011) {
    
    reduced.dataset <- HMD.dataset[HMD.dataset$Year >= minyear &
                                   HMD.dataset$Year <= maxyear &
                                   HMD.dataset$Age >= minage &
                                   HMD.dataset$Age <= maxage,
                                   c("Year", "Age", Sex)]
    Lexis.matrix <- tapply(X=reduced.dataset[[Sex]],
                           INDEX=list(Age=reduced.dataset$Age,
                               Year=reduced.dataset$Year),
                           FUN=sum)
    return(Lexis.matrix)
}
