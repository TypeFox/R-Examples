# Verified 1.3.18
diffmonths <-
function(date1, date2) {
        m1 = start(date1)
        m2 = start(date2)        
        if ( m1 > m2 ) { m = m2; m2=m1; m1=m }
        m1 = c(as.numeric(format(m1, format=c("%Y"))), as.numeric(format(m1, format=c("%m"))))
        m2 = c(as.numeric(format(m2, format=c("%Y"))), as.numeric(format(m2, format=c("%m"))))        
        y = (m2[1] - m1[1]) * 12
        if ( m2[2] < m1[2] ) { r = ( (y - 12) + (12 + m2[2] - m1[2]) - 1) }
        r = ( y + (m2[2] - m1[2]) )
        return(r)
}
