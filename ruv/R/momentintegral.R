momentintegral <-
function (d) 
{
    getf = function(d, x) {
        f = 1
        N = length(d)
        for (i in 1:N) {
            f = f/sqrt(1 + x * d[i])
        }
        return(f)
    }
    d = d/max(d)
    N = length(d)
    e = rep(0, N)
    x = c((0:1000)/1e+05, (100:1000)/10000, (100:1000)/1000, 
        (100:1000)/100, (100:1000)/10, (100:1000), (100:1000) * 
            10, (100:1000) * 100, (100:1000) * 1000, (100:1000) * 
            10000)
    y = getf(d, x)
    multiplier = 1
    int = (sum(y[1:1001]) - 0.5 * (y[1] + y[1001]))/1e+05
    nextsum = (sum(y[1002:1902]) - 0.5 * (y[1002] + y[1902]))/10000
    if (nextsum/int < 1e-06 & multiplier == 1) 
        multiplier = 1e+09
    else int = int + nextsum
    nextsum = (sum(y[1903:2803]) - 0.5 * (y[1903] + y[2803]))/1000
    if (nextsum/int < 1e-06 & multiplier == 1) 
        multiplier = 1e+08
    else int = int + nextsum
    nextsum = (sum(y[2804:3704]) - 0.5 * (y[2804] + y[3704]))/100
    if (nextsum/int < 1e-06 & multiplier == 1) 
        multiplier = 1e+07
    else int = int + nextsum
    nextsum = (sum(y[3705:4605]) - 0.5 * (y[3705] + y[4605]))/10
    if (nextsum/int < 1e-06 & multiplier == 1) 
        multiplier = 1e+06
    else int = int + nextsum
    nextsum = (sum(y[4606:5506]) - 0.5 * (y[4606] + y[5506]))
    if (nextsum/int < 1e-06 & multiplier == 1) 
        multiplier = 1e+05
    else int = int + nextsum
    nextsum = (sum(y[5507:6407]) - 0.5 * (y[5507] + y[6407])) * 
        10
    if (nextsum/int < 1e-06 & multiplier == 1) 
        multiplier = 10000
    else int = int + nextsum
    nextsum = (sum(y[6408:7308]) - 0.5 * (y[6408] + y[7308])) * 
        100
    if (nextsum/int < 1e-06 & multiplier == 1) 
        multiplier = 1000
    else int = int + nextsum
    nextsum = (sum(y[7309:8209]) - 0.5 * (y[7309] + y[8209])) * 
        1000
    if (nextsum/int < 1e-06 & multiplier == 1) 
        multiplier = 100
    else int = int + nextsum
    nextsum = (sum(y[8210:9110]) - 0.5 * (y[8210] + y[9110])) * 
        10000
    if (nextsum/int < 1e-06 & multiplier == 1) 
        multiplier = 10
    else int = int + nextsum
    d = d/multiplier
    x = c((0:10000)/1e+06, (1000:10000)/1e+05, (1000:10000)/10000, 
        (1000:10000)/1000, (1000:10000)/100, (1000:10000)/10, 
        (1000:10000), (1000:10000) * 10, (1000:10000) * 100, 
        (1000:10000) * 1000)
    y = getf(d, x)
    for (i in 1:N) {
        fully = y/(1 + x * d[i])
        int = (sum(fully[1:10001]) - 0.5 * (fully[1] + fully[10001]))/1e+06
        int = int + (sum(fully[10002:19002]) - 0.5 * (fully[10002] + 
            fully[19002]))/1e+05
        int = int + (sum(fully[19003:28003]) - 0.5 * (fully[19003] + 
            fully[28003]))/10000
        int = int + (sum(fully[28004:37004]) - 0.5 * (fully[28004] + 
            fully[37004]))/1000
        int = int + (sum(fully[37005:46005]) - 0.5 * (fully[37005] + 
            fully[46005]))/100
        int = int + (sum(fully[46006:55006]) - 0.5 * (fully[46006] + 
            fully[55006]))/10
        int = int + (sum(fully[55007:64007]) - 0.5 * (fully[55007] + 
            fully[64007]))
        int = int + (sum(fully[64008:73008]) - 0.5 * (fully[64008] + 
            fully[73008])) * 10
        int = int + (sum(fully[73009:82009]) - 0.5 * (fully[73009] + 
            fully[82009])) * 100
        int = int + (sum(fully[82010:91010]) - 0.5 * (fully[82010] + 
            fully[91010])) * 1000
        e[i] = int * d[i]
    }
    return(e/2)
}
