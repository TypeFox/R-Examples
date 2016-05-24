### R code from vignette source 'Follow-up.rnw'

###################################################
### code chunk number 1: Follow-up.rnw:65-67
###################################################
library(Epi)
print( sessionInfo(), l=F )


###################################################
### code chunk number 2: Follow-up.rnw:146-153
###################################################
data( nickel )
nicL <- Lexis( entry = list( per=agein+dob,
                             age=agein,
                             tfh=agein-age1st ),
                exit = list( age=ageout ),
         exit.status = ( icd %in% c(162,163) )*1,
                data = nickel )


###################################################
### code chunk number 3: Follow-up.rnw:163-166
###################################################
str( nickel )
str( nicL )
head( nicL )


###################################################
### code chunk number 4: Follow-up.rnw:175-176
###################################################
summary( nicL )


###################################################
### code chunk number 5: nicL1
###################################################
plot( nicL )


###################################################
### code chunk number 6: nicL2
###################################################
par( mar=c(3,3,1,1), mgp=c(3,1,0)/1.6 )
plot( nicL, 1:2, lwd=1, col=c("blue","red")[(nicL$exp>0)+1],
      grid=TRUE, lty.grid=1, col.grid=gray(0.7),
      xlim=1900+c(0,90), xaxs="i",
      ylim=  10+c(0,90), yaxs="i", las=1 )
points( nicL, 1:2, pch=c(NA,3)[nicL$lex.Xst+1],
        col="lightgray", lwd=3, cex=1.5 )
points( nicL, 1:2, pch=c(NA,3)[nicL$lex.Xst+1],
        col=c("blue","red")[(nicL$exp>0)+1], lwd=1, cex=1.5 )


###################################################
### code chunk number 7: Follow-up.rnw:229-232
###################################################
nicS1 <- splitLexis( nicL, "age", breaks=seq(0,100,10) )
summary( nicL )
summary( nicS1 )


###################################################
### code chunk number 8: Follow-up.rnw:239-240
###################################################
round( subset( nicS1, id %in% 8:10 ), 2 )


###################################################
### code chunk number 9: Follow-up.rnw:245-247
###################################################
nicS2 <- splitLexis( nicS1, "tfh", breaks=c(0,1,5,10,20,30,100) )
round( subset( nicS2, id %in% 8:10 ), 2 )


###################################################
### code chunk number 10: Follow-up.rnw:253-258
###################################################
timeBand( nicS2, "age", "middle" )[1:20]
# For nice printing and column labelling use the data.frame() function:
data.frame( nicS2[,c("id","lex.id","per","age","tfh","lex.dur")],
            mid.age=timeBand( nicS2, "age", "middle" ),
            mid.tfh=timeBand( nicS2, "tfh", "middle" ) )[1:20,]


###################################################
### code chunk number 11: Follow-up.rnw:276-281
###################################################
subset( nicL, id %in% 8:10 )
agehi <- nicL$age1st + 50 / nicL$exposure
nicC <- cutLexis( data=nicL, cut=agehi, timescale="age",
                  new.state=2, precursor.states=0 )
subset( nicC, id %in% 8:10 )


###################################################
### code chunk number 12: Follow-up.rnw:287-292
###################################################
subset( nicS2, id %in% 8:10 )
agehi <- nicS2$age1st + 50 / nicS2$exposure
nicS2C <- cutLexis( data=nicS2, cut=agehi, timescale="age",
                    new.state=2, precursor.states=0 )
subset( nicS2C, id %in% 8:10 )


###################################################
### code chunk number 13: Follow-up.rnw:312-321
###################################################
data( nickel )
nicL <- Lexis( entry = list( per=agein+dob,
                             age=agein,
                             tfh=agein-age1st ),
                exit = list( age=ageout ),
         exit.status = ( icd > 0 ) + ( icd %in% c(162,163) ),
                data = nickel )
summary( nicL )
subset( nicL, id %in% 8:10 )


###################################################
### code chunk number 14: Follow-up.rnw:325-333
###################################################
nicL <- Lexis( entry = list( per=agein+dob,
                             age=agein,
                             tfh=agein-age1st ),
                exit = list( age=ageout ),
         exit.status = ( icd > 0 ) + ( icd %in% c(162,163) ),
                data = nickel,
              states = c("Alive","D.oth","D.lung") )
summary( nicL )


###################################################
### code chunk number 15: Follow-up.rnw:345-353
###################################################
nicL$agehi <- nicL$age1st + 50 / nicL$exposure
nicC <- cutLexis( data = nicL,
                   cut = nicL$agehi,
             timescale = "age",
             new.state = "HiExp",
      precursor.states = "Alive" )
subset( nicC, id %in% 8:10 )
summary( nicC, scale=1000 )


###################################################
### code chunk number 16: Follow-up.rnw:371-379
###################################################
nicC <- cutLexis( data = nicL,
                   cut = nicL$agehi,
             timescale = "age",
             new.state = "Hi",
           split.states=TRUE, new.scale=TRUE,
      precursor.states = "Alive" )
subset( nicC, id %in% 8:10 )
summary( nicC, scale=1000 )


