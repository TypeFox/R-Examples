`window_chron` <-
function(x, day1,hour1,day2,hour2, ...){
    window(x, start = chron(day1, hour1), end = chron(day2, hour2))
}

