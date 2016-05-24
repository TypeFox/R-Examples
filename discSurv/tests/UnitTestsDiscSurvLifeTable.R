library(discSurv)

##################
# lifeTable

# Construct data set with 100 persons
TimeInput <- c(rep(1, 20), rep(2, 30), rep(3, 10), rep(4, 15), rep(5, 8), rep(6, 12), rep(7, 5))
EventInput <- c(rep(1, 10), rep(0, 10),
                rep(1, 12), rep(0, 18),
                rep(1, 7), rep(0, 3),
                rep(1, 15),
                rep(0, 8), 
                rep(1, 6), rep(0, 6),
                rep(1, 2), rep(0, 3))
genData <- data.frame(Time=TimeInput, Event=EventInput)
LifeTab <- lifeTable (dataSet=genData, timeColumn="Time", censColumn="Event")
LifeTab

# Checks
MaxTimeSeq <- 1:max(genData [, "Time"])+1
intervalBordersInput <- paste("[", c(0, MaxTimeSeq [-length(MaxTimeSeq)]), ", ", MaxTimeSeq, ")", sep = "")
LifeTab2 <- lifeTable (dataSet=genData, timeColumn="Time", censColumn="Event", intervalBorders=intervalBordersInput)
LifeTab2
stopifnot(all(rownames(LifeTab2)==intervalBordersInput [-length(intervalBordersInput)]))

# Additional checks
stopifnot(LifeTab$Output [1, 1]==dim(genData) [1])
SollN <- c(dim(genData) [1], 
           dim(genData) [1]-20, 
           dim(genData) [1]-50,
           dim(genData) [1]-60,
           dim(genData) [1]-75,
           dim(genData) [1]-83,
           dim(genData) [1]-95)
stopifnot(all(LifeTab$Output [, 1]==SollN))
stopifnot(identical(rownames(LifeTab$Output), paste("[", 0:6, ", ", 1:7, ")", sep="")))
stopifnot(LifeTab$Output [, "events"]==c(10, 12, 7, 15, 0, 6, 2))
stopifnot(LifeTab$Output [, "dropouts"]==c(10, 18, 3, 0, 8, 6, 3))
stopifnot(all(LifeTab$Output [, "atRisk"]==(SollN-c(10, 18, 3, 0, 8, 6, 3)/2)))
stopifnot(all.equal(LifeTab$Output [, "hazard"], c(10, 12, 7, 15, 0, 6, 2)/(SollN-c(10, 18, 3, 0, 8, 6, 3)/2)))
SollHaz <- c(10, 12, 7, 15, 0, 6, 2)/(SollN-c(10, 18, 3, 0, 8, 6, 3)/2)
SollSeHaz <- sqrt((SollHaz-SollHaz^2)/(SollN-c(10, 18, 3, 0, 8, 6, 3)/2))
stopifnot(all.equal(LifeTab$Output [, "seHazard"], SollSeHaz))
stopifnot(all.equal(LifeTab$Output [, "S"], cumprod(1-SollHaz)))
stopifnot(all.equal(LifeTab$Output [, "seS"], cumprod(1-SollHaz)*sqrt(cumsum(SollHaz/((1-SollHaz)*(SollN-c(10, 18, 3, 0, 8, 6, 3)/2)))), 4))
SollSeCumHaz <- sqrt(cumsum(c(10, 12, 7, 15, 0, 6, 2)/(SollN-c(10, 18, 3, 0, 8, 6, 3)/2)^2))
stopifnot(all.equal(LifeTab$Output [, "cumHazard"], cumsum(SollHaz)))
stopifnot(all.equal(LifeTab$Output [, "seCumHazard"], SollSeCumHaz))

# Additional Check if atRiskInput is corrected by last events, dropouts in no information intervals
# Check if also intervals before the first observation are generated
CheckDupRows <- data.frame(Time=rep(c(rep(10,9),10:24, 30:34, 36), 
                                    times=c(208, 185, 171,  92,  95,  20,  29, 209,  32,  
                                            25,  28,  17,  17,  20,  10,  11,  18,   7,   4,  24,   3,   3,   
                                            4,   5,   5,   5,   1,   4,   3,  15)), 
                           Censor=c(rep(0, 100), rep(1, 108), 
                                    rep(0, 7), rep(1, 178), 
                                    rep(0, 12), rep(1, 159), 
                                    rep(0, 3), rep(1, 89), 
                                    rep(0, 9), rep(1, 86),
                                    rep(1, 20),
                                    rep(0, 4), rep(1, 25),
                                    rep(0, 21), rep(1, 188),
                                    rep(0, 2), rep(1, 30),
                                    rep(0, 3), rep(1, 22),
                                    rep(0, 2), rep(1, 26),
                                    rep(0, 1), rep(1, 16),
                                    rep(0, 2), rep(1, 15),
                                    rep(0, 2), rep(1, 18),
                                    rep(0, 3), rep(1, 7),
                                    rep(0, 3), rep(1, 8),
                                    rep(0, 6), rep(1, 12),
                                    rep(0, 1), rep(1, 6),
                                    rep(1, 4),
                                    rep(0, 8), rep(1, 16),
                                    rep(1, 3),
                                    rep(0, 1), rep(1, 2),
                                    rep(0, 1), rep(1, 3),
                                    rep(0, 1), rep(1, 4),
                                    rep(0, 2), rep(1, 3),
                                    rep(0, 2), rep(1, 3),
                                    rep(1, 1),
                                    rep(0, 2), rep(1, 2),
                                    rep(1, 3), 
                                    rep(1, 15)))
CheckTab <- lifeTable(CheckDupRows, "Time", "Censor")
CheckIndicator <- all(CheckTab$Output [1:10, "n"] == 1270, CheckTab$Output [25:29, "n"] == 33, CheckTab$Output [35, "n"] == 15)
stopifnot(CheckIndicator)