cal <- function(month, year) {

    yyy <- FALSE

    if(missing(year) && missing(month)) {  # no args, use current month
        tmp <- as.POSIXlt(Sys.time())
        year <- tmp$year+1900
        month <- tmp$mon+1
    } else if( missing(year) && is.numeric(month) && month > 12 ) { # switch month to year
        year <- month
        yyy <- TRUE
    } else if( missing(year) ) { # use current year
        tmp <- as.POSIXlt(Sys.time())
        year <- tmp$year+1900
    } else if( missing(month) ) { # no month do year
        yyy <- TRUE
    }

    if(yyy) {  # year calendar
        par(mfrow=c(4,3),oma=c(0,0,3.5,0))
        tmp <- seq( from=ISOdate(year,1,1), to=ISOdate(year,12,31), by='days' )
        tmp2 <- as.POSIXlt(tmp)
        wd <- tmp2$wd
        par(mar=c(1.5,1.5,2.5,1.5))
        for(i in 1:12){
            w <- (tmp2$mon+1) == i
            cs <- cumsum(wd[w]==0)
            if(cs[1] > 0) cs <- cs - 1
            nr <- max( cs ) + 1
            plot.new()
            plot.window( xlim=c(0,6), ylim=c(0,nr+1) )
            text( wd[w], nr - cs -0.5 , tmp2$mday[w] )
            title( main=month.name[i] )
            text( 0:6, nr+0.5, c('S','M','T','W','T','F','S') )
            mtext( year, outer=TRUE, line=1, cex=2 )
        }

    } else {  # month calendar
        if( is.character(month) ) {
            tmp <- pmatch( tolower(month), tolower(month.name) )
            if( is.na(tmp) ) {
                tmp <- pmatch( month, as.character(1:12))
            }
            if( is.na(tmp) ) {
                warning('Unable to match month, using current month')
                tmp <- as.POSIXlt(Sys.time())
                month <- tmp$mon+1
            } else {
                month <- tmp
            }
        }
        ld <- seq( from=ISOdate(year,month,1), length=2, by='months')[2]-86400
        days <- seq( from=ISOdate(year,month,1), to=ld, by='days')
        tmp <- as.POSIXlt(days)
        wd <- tmp$wday
        cs <- cumsum(wd == 0)
        if(cs[1] > 0) cs <- cs - 1
        nr <- max(cs) + 1
        par(oma=c(0.1,0.1,4.6,0.1))
        par(mfrow=c(nr,7))
        par(mar=c(0,0,0,0))
        for(i in seq_len(wd[1])){
            plot.new()
                                       # box()
        }
        day.name <- c('Sun','Mon','Tues','Wed','Thur','Fri','Sat')
        for(i in tmp$mday){
            plot.new()
            box()
            text(0,1, i, adj=c(0,1))
            if(i < 8) mtext( day.name[wd[i]+1], line=0.5,
                            at=grconvertX(0.5,to='ndc'), outer=TRUE )
        }
        mtext(month.name[month], line=2.5, at=0.5, cex=1.75, outer=TRUE)
                                        #box('inner') #optional

        invisible(function(day) {
            day <- day + wd[1] - 1
            rr <- day %/% 7 + 1
            cc <- day %%  7 + 1
            par(mfg=c(rr,cc))
        })
    }

}

### cal(10,2011)
### par(mfg=c(3,2))  # monday oct 10
### text(.5,.5, 'Some\nText', cex=2)
###
### par(mfg=c(2,3)) #Tues oct 4
### text(1,1, 'Top Right', adj=c(1,1))
###
### par(mfg=c(2,4)) # Wed oct 5
### text(0,0, 'Bottom Left', adj=c(0,0))
###
### par(mfg=c(6,2)) # oct 31
### tmp.x <- runif(25)
### tmp.y <- rnorm(25,tmp.x,.1)
### par(usr=c( range(tmp.x), range(tmp.y) ) )
### points(tmp.x,tmp.y)
###





