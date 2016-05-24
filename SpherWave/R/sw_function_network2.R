#Latitude (shown as a horizontal line) is the angular distance, in degrees, minutes, 
#and seconds of a point north or south of the Equator. Lines of latitude are often referred to as parallels. 
#Longitude (shown as a vertical line) is the angular distance, in degrees, minutes, 
#and seconds, of a point east or west of the Prime (Greenwich) Meridian. Lines of longitude are often referred 
#to as meridians. 


### reg.grid : Modified Gottlemann Regular Grid
# modnetlevel, modterritory, modregcenter
# distboun.comp, secdistboun.comp, multiselc.comp
"reg.grid" <-
function (x, latlon) 
{
    nlevel <- modnetlevel(x * pi/180)$nlevel
    terri <- modterritory(1)$mean/2
    netlab <- distboun.comp(1, latlon, modregcenter(1)$center, terri)$netlab
    boun <- modterritory(2)$min/4
    for (i in 2:nlevel) {
        netlab <- secdistboun.comp(i, latlon, modregcenter(i)$center, boun/2^(i - 2), netlab)$netlab
        netlab <- multiselc.comp(i, latlon, terri/2^(i - 1), netlab)$netlab
    }
    for (i in 1:length(netlab)) 
        if (netlab[i] == 0) 
            netlab[i] <- nlevel + 1
        
    list(netlab = netlab)
}

# Decide the number of levels
"modnetlevel" <-
function (angle) 
{
    l <- 1
    k <- seq(0, 2^l)
    i <- seq(0, 2^(l + 1))
    A <- k * (pi/2^l)
    B <- i * (pi/2^l)
    A <- A - (pi/2)
    B <- B - pi
    A <- 90 * A/(pi/2)
    B <- 180 * B/pi
    geod <- NULL
    for (j in 1:ceiling(length(A)/2)) {
        theta <- 90 - A[j]
        theta <- theta * pi/180
        phi <- B[1] * pi/180
        dot <- cos(theta) * cos((90 - A[j]) * pi/180)
        dot <- dot + sin(theta) * sin((90 - A[j]) * pi/180) * cos(phi - B[2] * pi/180)
        geod <- c(geod, acos(dot))
    }
    cvalue <- max(geod)
    while (cvalue > angle) {
        l <- l + 1
        k <- seq(0, 2^l)
        i <- seq(0, 2^(l + 1))
        A <- k * (pi/2^l)
        B <- i * (pi/2^l)
        A <- A - (pi/2)
        B <- B - pi
        A <- 90 * A/(pi/2)
        B <- 180 * B/pi
        geod <- NULL
        for (j in 1:ceiling(length(A)/2)) {
            theta <- 90 - A[j]
            theta <- theta * pi/180
            phi <- B[1] * pi/180
            dot <- cos(theta) * cos((90 - A[j]) * pi/180)
            dot <- dot + sin(theta) * sin((90 - A[j]) * pi/180) * cos(phi - B[2] * pi/180)
            geod <- c(geod, acos(dot))
        }
        cvalue <- max(geod)
    }
    nlevel <- l - 1
    
    list(nlevel = nlevel)
}

# Calculate Territories
"modterritory" <-
function (l) 
{
    k <- seq(0, 2^l)
    i <- seq(0, 2^(l + 1))
    A <- k * (pi/2^l)
    B <- i * (pi/2^l)
    A <- A - (pi/2)
    B <- B - pi
    A <- 90 * A/(pi/2)
    B <- 180 * B/pi
    dist <- NULL
    dd <- ceiling(length(A)/2)
    for (j in 2:dd) {
        theta <- 90 - A[j]
        theta <- theta * pi/180
        phi <- B[1] * pi/180
        dot <- cos(theta) * cos((90 - A[j]) * pi/180)
        dot <- dot + sin(theta) * sin((90 - A[j]) * pi/180) * cos(phi - B[2] * pi/180)
        dist <- c(dist, acos(dot))
    }

    list(dist = dist, max = max(dist), min = min(dist), mean = mean(dist), median = median(dist))
}

# Generate center points from regular grid
"modregcenter" <-
function (l) 
{
    center <- NULL
    for (k in 0:(2^l/2 - 1)) {
        a <- k * (pi/2^l)
        aa <- (k + 1) * (pi/2^l)
        b <- NULL
        for (i in 0:2^(l + 1)) 
            b <- c(b, i * (pi/2^l))
        
        for (j in 1:(length(b) - 1))
            center <- rbind(center, c((a + aa)/2, (b[j] + b[j + 1])/2))
    }
    center <- cbind(center[, 1] - (pi/2), center[, 2] - pi)
    center <- cbind(90 * center[, 1]/(pi/2), 180 * center[, 2]/pi)
    center <- cbind(c(-center[, 1], center[, 1]), c(center[, 2], center[, 2]))
    
    list(center = center)
}

# Select stations in initial stage
"distboun.comp" <-
function (l, xm, center, terri) 
{
    m <- length(xm[, 1])
    n <- length(center[, 1])
    netlab <- rep(0, m)
    dist <- rep(0, m)
    dot <- rep(0, m)
    theta <- 0
    phi <- 0
    min <- 0
    spot <- 0
    zz <- .Fortran("disboun", PACKAGE = "SpherWave", xm = as.double(xm), center = as.double(center), 
        dot = as.double(dot), dist = as.double(dist), min = as.double(min), 
        netlab = as.integer(netlab), theta = as.double(theta), 
        phi = as.double(phi), terri = as.double(terri), spot = as.integer(spot), 
        m = as.integer(m), n = as.integer(n), l = as.integer(l))
        
    list(netlab = zz$netlab)
}

# Select stations within some territories (L>2)
"secdistboun.comp" <-
function (l, xm, center, terri, netlab) 
{
    m <- length(xm[, 1])
    n <- length(center[, 1])
    dist <- rep(0, m)
    dot <- rep(0, m)
    theta <- 0
    phi <- 0
    min <- 0
    spot <- 0
    zz <- .Fortran("sedisboun", PACKAGE = "SpherWave", xm = as.double(xm), center = as.double(center), 
        dot = as.double(dot), dist = as.double(dist), min = as.double(min), 
        netlab = as.integer(netlab), theta = as.double(theta), 
        phi = as.double(phi), terri = as.double(terri), spot = as.integer(spot), 
        m = as.integer(m), n = as.integer(n), l = as.integer(l))
        
    list(netlab = zz$netlab)
}

# Delete stations which is close to stations selected in (L>2)
"multiselc.comp" <-
function (l, xm, bb, netlab) 
{
    for (j in 1:(l - 1)) {
        m <- length(xm[, 1])
        ym <- xm[netlab == l - j, ]
        n <- length(ym[, 1])
        dist <- rep(0, m)
        dot <- rep(0, m)
        theta <- 0
        phi <- 0
        zz <- .Fortran("seselc", PACKAGE = "SpherWave", xm = as.double(xm), ym = as.double(ym), 
            dot = as.double(dot), dist = as.double(dist), netlab = as.integer(netlab), 
            theta = as.double(theta), phi = as.double(phi), bb = as.double(bb), 
            m = as.integer(m), n = as.integer(n), l = as.integer(l))
        netlab <- zz$netlab
    }
    
    list(netlab = zz$netlab)
}

### red.grid : Modified Gottlemann Reduced Grid
# modnetlevel, modterritory, modredcenter
# distboun.comp, secdistboun.comp, multiselc.comp
"red.grid" <-
function (x, latlon) 
{
    nlevel <- modnetlevel(x * pi/180)$nlevel
    terri <- modterritory(1)$mean/2
    netlab <- distboun.comp(1, latlon, modredcenter(1)$center, terri)$netlab
    boun <- modterritory(2)$min/4
    for (i in 2:nlevel) {
        netlab <- secdistboun.comp(i, latlon, modredcenter(i)$center, boun/2^(i - 2), netlab)$netlab
        netlab <- multiselc.comp(i, latlon, terri/2^(i - 1), netlab)$netlab
    }
    for (i in 1:length(netlab)) 
        if (netlab[i] == 0) 
            netlab[i] <- nlevel + 1

    list(netlab = netlab)
}

# Generate center points from reduced grid
"modredcenter" <-
function (l) 
{
    center <- NULL
    for (k in 0:(2^l/2 - 1)) {
        a <- k * (pi/2^l)
        
        if (pi/4 <= a && a <= 3 * pi/4)
            r <- 0
        if (0 < a && a < pi/4)
            r <- floor(l - logb(pi * k, 2))
        if (3 * pi/4 < a && a < pi)
            r <- floor(l - logb(pi * (2^l - k), 2))
        if (a == 0 || a == pi)
            r <- l - 1

        rr <- (2^(l + 1) - 2^r)/2^r + 1
        
        b <- NULL
        for (i in 0:rr) 
            b <- c(b, i * 2^r * pi/2^l)

        aa <- (k + 1) * (pi/2^l)
        
        for (j in 1:(length(b) - 1))
            center <- rbind(center, c((a + aa)/2, (b[j] + b[j + 1])/2))

    }
    center <- cbind(center[, 1] - (pi/2), center[, 2] - pi)
    center <- cbind(90 * center[, 1]/(pi/2), 180 * center[, 2]/pi)
    center <- cbind(c(-center[, 1], center[, 1]), c(center[, 2], center[, 2]))
    
    list(center = center)
}

### gotreg.grid : Gottlemann Regular Grid
# gotnetlevel, newgotterritory, gotregcenter
# distboun.comp, secdistboun.comp, multiselc.comp
"gotreg.grid" <-
function (x, latlon) 
{
    nlevel <- gotnetlevel(x * pi/180)$nlevel
    terri <- newgotterritory(1)$mean/2
    netlab <- distboun.comp(1, latlon, gotregcenter(1)$center, terri)$netlab
    boun <- newgotterritory(2)$min/4
    for (i in 2:nlevel) {
        netlab <- secdistboun.comp(i, latlon, gotregcenter(i)$center, boun/2^(i - 2), netlab)$netlab
        netlab <- multiselc.comp(i, latlon, terri/2^(i - 1), netlab)$netlab
    }
    for (i in 1:length(netlab)) 
        if (netlab[i] == 0)
            netlab[i] <- nlevel + 1

    list(netlab = netlab)
}

# Decide the number of levels
"gotnetlevel" <-
function (angle) 
{
    l <- 1
    k <- seq(0, 2^(l + 1))
    i <- seq(0, 2^(l + 2))
    A <- k * (pi/2^(l + 1))
    B <- i * (pi/2^(l + 1))
    A <- A - (pi/2)
    B <- B - pi
    A <- 90 * A/(pi/2)
    B <- 180 * B/pi
    geod <- NULL
    for (j in 1:ceiling(length(A)/2)) {
        theta <- 90 - A[j]
        theta <- theta * pi/180
        phi <- B[1] * pi/180
        dot <- cos(theta) * cos((90 - A[j]) * pi/180)
        dot <- dot + sin(theta) * sin((90 - A[j]) * pi/180) * cos(phi - B[2] * pi/180)
        geod <- c(geod, acos(dot))
    }
    cvalue <- max(geod)
    while (cvalue > angle) {
        l <- l + 1
        k <- seq(0, 2^(l + 1))
        i <- seq(0, 2^(l + 2))
        A <- k * (pi/2^(l + 1))
        B <- i * (pi/2^(l + 1))
        A <- A - (pi/2)
        B <- B - pi
        A <- 90 * A/(pi/2)
        B <- 180 * B/pi
        geod <- NULL
        for (j in 1:ceiling(length(A)/2)) {
            theta <- 90 - A[j]
            theta <- theta * pi/180
            phi <- B[1] * pi/180
            dot <- cos(theta) * cos((90 - A[j]) * pi/180)
            dot <- dot + sin(theta) * sin((90 - A[j]) * pi/180) * cos(phi - B[2] * pi/180)
            geod <- c(geod, acos(dot))
        }
        cvalue <- max(geod)
    }
    nlevel <- l - 1
    
    list(nlevel=nlevel)
}

# Calculate Territories
"newgotterritory" <-
function (l) 
{
    out <- NULL
    k <- seq(0, 2^(l + 1))
    i <- seq(0, 2^(l + 2))
    A <- k * (pi/2^(l + 1))
    B <- i * (pi/2^(l + 1))
    A <- A - (pi/2)
    B <- B - pi
    A <- 90 * A/(pi/2)
    B <- 180 * B/pi
    dist <- NULL
    dd <- ceiling(length(A)/2)
    for (j in 2:dd) {
        theta <- 90 - A[j]
        theta <- theta * pi/180
        phi <- B[1] * pi/180
        dot <- cos(theta) * cos((90 - A[j]) * pi/180)
        dot <- dot + sin(theta) * sin((90 - A[j]) * pi/180) * cos(phi - B[2] * pi/180)
        dist <- c(dist, acos(dot))
    }
    
    list(dist = dist, max = max(dist), min = min(dist), mean = mean(dist), median = median(dist))
}

# Generate center points from regular grid
"gotregcenter" <-
function (l) 
{
    center <- NULL
    for (k in 0:(2^(l + 1)/2 - 1)) {
        a <- k * (pi/2^(l + 1))
        aa <- (k + 1) * (pi/2^(l + 1))
        b <- NULL
        for (i in 0:2^(l + 2)) 
            b <- c(b, i * (pi/2^(l + 1)))
        
        for (j in 1:(length(b) - 1)) 
            center <- rbind(center, c((a + aa)/2, (b[j] + b[j + 1])/2))
    }
    center <- cbind(center[, 1] - (pi/2), center[, 2] - pi)
    center <- cbind(90 * center[, 1]/(pi/2), 180 * center[, 2]/pi)
    center <- cbind(c(-center[, 1], center[, 1]), c(center[, 2], center[, 2]))
    
    list(center = center)
}

### gotred.grid : Gottlemann Reduced Grid
# gotnetlevel, newgotterritory, gotredcenter
# distboun.comp, secdistboun.comp, multiselc.comp
"gotred.grid" <-
function (x, latlon) 
{
    nlevel <- gotnetlevel(x * pi/180)$nlevel
    terri <- newgotterritory(1)$mean/2
    netlab <- distboun.comp(1, latlon, gotredcenter(1)$center, terri)$netlab
    boun <- newgotterritory(2)$min/4
    for (i in 2:nlevel) {
        netlab <- secdistboun.comp(i, latlon, gotredcenter(i)$center, boun/2^(i - 2), netlab)$netlab
        netlab <- multiselc.comp(i, latlon, terri/2^(i - 1), netlab)$netlab
    }
    for (i in 1:length(netlab))
        if (netlab[i] == 0) 
            netlab[i] <- nlevel + 1
    
    list(netlab = netlab)
}

# Generate center points from reduced grid
"gotredcenter" <-
function (l) 
{
    center <- NULL
    for (k in 0:(2^(l + 1)/2 - 1)) {
        a <- k * (pi/2^(l + 1))
        
        if (pi/4 <= a && a <= 3 * pi/4)
            r <- 0
        if (0 < a && a < pi/4) 
            r <- floor((l + 1) - logb(pi * k, 2))
        if (3 * pi/4 < a && a < pi) 
            r <- floor((l + 1) - logb(pi * (2^l - k), 2))
        if (a == 0 || a == pi) 
            r <- l - 1
        
        rr <- (2^(l + 2) - 2^r)/2^r + 1
        
        b <- NULL
        for (i in 0:rr) 
            b <- c(b, i * 2^r * pi/2^(l + 1))

        aa <- (k + 1) * (pi/2^(l + 1))
        
        for (j in 1:(length(b) - 1)) 
            center <- rbind(center, c((a + aa)/2, (b[j] + b[j + 1])/2))
    }
    center <- cbind(center[, 1] - (pi/2), center[, 2] - pi)
    center <- cbind(90 * center[, 1]/(pi/2), 180 * center[, 2]/pi)
    center <- cbind(c(-center[, 1], center[, 1]), c(center[, 2], center[, 2]))
    
    list(center = center)
}

### hsreg.grid : Oh Regular Grid
# netlevel, newterritory, regcenter
# distboun.comp, secdistboun.comp, multiselc.comp
"hsreg.grid" <-
function (x, latlon) 
{
    nlevel <- netlevel(x * pi/180)$nlevel
    terri <- newterritory(1)$mean/2
    netlab <- distboun.comp(1, latlon, regcenter(1)$center, terri)$netlab
    boun <- newterritory(2)$min/4
    for (i in 2:nlevel) {
        netlab <- secdistboun.comp(i, latlon, regcenter(i)$center, boun/2^(i - 2), netlab)$netlab
        netlab <- multiselc.comp(i, latlon, terri/2^(i - 1), netlab)$netlab
    }
    for (i in 1:length(netlab))
        if (netlab[i] == 0)
            netlab[i] <- nlevel + 1
    
    list(netlab = netlab)
}

# Decide the number of levels
"netlevel" <-
function (angle) 
{
    l <- 1
    A <- seq(-90, 90, length = 1 + 2^l)
    B <- seq(-180, 180, length = 1 + 2^l)
    geod <- NULL
    for (j in 1:ceiling(length(A)/2)) {
        theta <- 90 - A[j]
        theta <- theta * pi/180
        phi <- B[1] * pi/180
        dot <- cos(theta) * cos((90 - A[j]) * pi/180)
        dot <- dot + sin(theta) * sin((90 - A[j]) * pi/180) * cos(phi - B[2] * pi/180)
        geod <- c(geod, acos(dot))
    }
    cvalue <- max(geod)
    while (cvalue > angle) {
        l <- l + 1
        A <- seq(-90, 90, length = 1 + 2^l)
        B <- seq(-180, 180, length = 1 + 2^l)
        geod <- NULL
        for (j in 1:ceiling(length(A)/2)) {
            theta <- 90 - A[j]
            theta <- theta * pi/180
            phi <- B[1] * pi/180
            dot <- cos(theta) * cos((90 - A[j]) * pi/180)
            dot <- dot + sin(theta) * sin((90 - A[j]) * pi/180) * cos(phi - B[2] * pi/180)
            geod <- c(geod, acos(dot))
        }
        cvalue <- max(geod)
    }
    nlevel <- l - 1
    
    list(nlevel = nlevel)
}

# Calculate Territories
"newterritory" <-
function (l) 
{
    A <- seq(-90, 90, length = 1 + 2^l)
    B <- seq(-180, 180, length = 1 + 2^l)
    dist <- NULL
    dd <- ceiling(length(A)/2)
    for (j in 2:dd) {
        theta <- 90 - A[j]
        theta <- theta * pi/180
        phi <- B[1] * pi/180
        dot <- cos(theta) * cos((90 - A[j]) * pi/180)
        dot <- dot + sin(theta) * sin((90 - A[j]) * pi/180) * cos(phi - B[2] * pi/180)
        dist <- c(dist, acos(dot))
    }
    
    list(dist = dist, max = max(dist), min = min(dist), mean = mean(dist), median = median(dist))
}

# Generate center points from regular grid
"regcenter" <-
function (l) 
{
    center <- NULL
    for (k in 0:(2^l/2 - 1)) {
        a <- k * (pi/2^l)
        aa <- (k + 1) * (pi/2^l)
        b <- NULL
        for (i in 0:2^l)
            b <- c(b, i * (pi/2^(l - 1)))

        for (j in 1:(length(b) - 1))
            center <- rbind(center, c((a + aa)/2, (b[j] + b[j + 1])/2))
    }
    center <- cbind(center[, 1] - (pi/2), center[, 2] - pi)
    center <- cbind(90 * center[, 1]/(pi/2), 180 * center[, 2]/pi)
    center <- cbind(c(-center[, 1], center[, 1]), c(center[, 2], center[, 2]))
    
    list(center = center)
}

### hsred.grid : Oh Reduced Grid
# netlevel, newterritory, redcenter
# distboun.comp, secdistboun.comp, multiselc.comp
"hsred.grid" <-
function (x, latlon) 
{
    nlevel <- netlevel(x * pi/180)$nlevel
    terri <- newterritory(1)$mean/2
    netlab <- distboun.comp(1, latlon, redcenter(1)$center, terri)$netlab
    boun <- newterritory(2)$min/4
    for (i in 2:nlevel) {
        netlab <- secdistboun.comp(i, latlon, redcenter(i)$center, boun/2^(i - 2), netlab)$netlab
        netlab <- multiselc.comp(i, latlon, terri/2^(i - 1), netlab)$netlab
    }
    for (i in 1:length(netlab))
        if (netlab[i] == 0)
            netlab[i] <- nlevel + 1

    list(netlab = netlab)
}

# Generate center points from reduced grid
"redcenter" <-
function (l) 
{
    center <- NULL
    for (k in 0:(2^l/2 - 1)) {
        a <- k * (pi/2^l)
        
        if (pi/4 <= a && a <= 3 * pi/4) 
            r <- 0
        if (0 < a && a < pi/4) 
            r <- round(l - logb(pi * k, 2))     
        if (3 * pi/4 < a && a < pi)
            r <- round(l - logb(pi * (2^l - k), 2))
        if (a == 0 || a == pi)
            r <- l - 1
        
        rr <- (2^l - 2^r)/2^r + 1
        
        b <- NULL
        for (i in 0:rr)
            b <- c(b, i * 2^r * pi/2^(l - 1))

        aa <- (k + 1) * (pi/2^l)
        
        for (j in 1:(length(b) - 1))
            center <- rbind(center, c((a + aa)/2, (b[j] + b[j + 1])/2))
    }
    center <- cbind(center[, 1] - (pi/2), center[, 2] - pi)
    center <- cbind(90 * center[, 1]/(pi/2), 180 * center[, 2]/pi)
    center <- cbind(c(-center[, 1], center[, 1]), c(center[, 2], center[, 2]))
        
    list(center = center)
}

"network.design" <- 
function (latlon, method="Oh", type="reduce", nlevel, x) 
{
    if (method != "cover" && method != "Gottlemann" && method != "ModifyGottlemann" && method != "Oh")
        stop("Only cover, Gottlemann and Oh's grid are provided.")
    if ((method == "ModifyGottlemann" || method == "Gottlemann" || method =="Oh") && (type != "reduce" && type !="regular"))
         stop("Specify type of the grid, reduce or regular.")
    if (method == "cover" && nrow(latlon) != sum(nlevel))
        stop("The number of data is not the same to the total number of network label.")
    if (method == "cover" && !all(diff(nlevel) < 0))
        stop("The number of data is not the same to the total number of network label.")
    
    if (method == "cover") {
        nlevels <- length(nlevel)
        netlab <- rep(1, nrow(latlon))
        assign(paste("m", nlevels, sep=""),
            cover.design(cbind(latlon[, 2], latlon[, 1]), nlevel[nlevels], DIST=rdist.earth))       
                  
        for (i in 2:(nlevels - 1))
            assign(paste("m", nlevels + 1 - i, sep=""),
                cover.design(cbind(latlon[, 2], latlon[, 1]), nlevel[nlevels + 1 - i], DIST=rdist.earth, 
                    fixed=get(paste("m", nlevels + 2 - i, sep=""))$best.id))      
        for (i in 2:nlevels)
            netlab[get(paste("m", i, sep=""))$best.id] <- i
    }
    if (method == "Gottlemann" && type == "regular") {
        netlab <- gotreg.grid(x, latlon)$netlab
        netlab <- (max(netlab) + 1) - netlab
    }
    if (method == "Gottlemann" && type == "reduce") {
        netlab <- gotred.grid(x, latlon)$netlab
        netlab <- (max(netlab) + 1) - netlab
    }        
    if (method == "ModifyGottlemann" && type == "regular") {
        netlab <- reg.grid(x, latlon)$netlab
        netlab <- (max(netlab) + 1) - netlab
    }
    if (method == "ModifyGottlemann" && type == "reduce") {
        netlab <- red.grid(x, latlon)$netlab
        netlab <- (max(netlab) + 1) - netlab
    }   
    if (method == "Oh" && type == "regular") {
        netlab <- hsreg.grid(x, latlon)$netlab
        netlab <- (max(netlab) + 1) - netlab
    }
    if (method == "Oh" && type == "reduce") {
        netlab <- hsred.grid(x, latlon)$netlab
        netlab <- (max(netlab) + 1) - netlab
    }      
    netlab
}

#=================================
# Generating Grid
#  Note A=latitute and B=longitude
#  Rectangular grid - reggrid, redgrid
#  Modified Gottelmann grid - modreggrid, modredgrid
#  Gottelmann grid - gotreggrid, gotredgrid
#=================================

"reggrid" <-
function (l) 
{
    k <- seq(0, 2^l)
    i <- k
    A <- k * (pi/2^l)
    B <- i * (pi/2^(l - 1))
    A <- A - (pi/2)
    B <- B - pi
    A <- 90 * A/(pi/2)
    B <- 180 * B/pi
    grid <- NULL
    for (i in 1:length(A))
        for (j in 1:length(B)) 
            grid <- rbind(grid, c(A[i], B[j]))

    list(grid = grid)
}

"redgrid" <-
function (l) 
{
    points <- NULL
    for (k in 0:2^l) {
        a <- k * (pi/2^l)
        if (pi/4 <= a && a <= 3 * pi/4)
            r <- 0
        
        if (0 < a && a < pi/4) 
            r <- round(l - logb(pi * k, 2))

        if (3 * pi/4 < a && a < pi) 
            r <- round(l - logb(pi * (2^l - k), 2))

        if (a == 0 || a == pi)
            r <- (l - 1)

        rr <- (2^l - 2^r)/2^r + 1
        for (i in 0:rr) {
            b <- i * 2^r * pi/2^(l - 1)
            points <- rbind(points, c(a, b))
        }
    }
    grid <- cbind(points[, 1] - (pi/2), points[, 2] - pi)
    grid <- cbind(90 * grid[, 1]/(pi/2), 180 * grid[, 2]/pi)
    
    list(grid = grid)
}

"modreggrid" <-
function (l) 
{
    k <- seq(0, 2^l)
    i <- seq(0, 2^(l + 1))
    A <- k * (pi/2^l)
    B <- i * (pi/2^l)
    A <- A - (pi/2)
    B <- B - pi
    A <- 90 * A/(pi/2)
    B <- 180 * B/pi
    grid <- NULL
    for (i in 1:length(A))
        for (j in 1:length(B))
            grid <- rbind(grid, c(A[i], B[j]))

    list(grid = grid)
}

"modredgrid" <-
function (l) 
{
    points <- NULL
    for (k in 0:2^l) {
        a <- k * (pi/2^l)
        if (pi/4 <= a && a <= 3 * pi/4)
            r <- 0

        if (0 < a && a < pi/4) 
            r <- floor((l) - logb(pi * k, 2))

        if (3 * pi/4 < a && a < pi) 
            r <- floor((l) - logb(pi * (2^l - k), 2))

        if (a == 0 || a == pi)
            r <- (l - 1)

        rr <- (2^(l + 1) - 2^r)/2^r + 1
        for (i in 0:rr) {
            b <- i * 2^r * pi/2^l
            points <- rbind(points, c(a, b))
        }
    }
    grid <- cbind(points[, 1] - (pi/2), points[, 2] - pi)
    grid <- cbind(90 * grid[, 1]/(pi/2), 180 * grid[, 2]/pi)
    
    list(grid = grid)
}

"gotreggrid" <-
function (l) 
{
    k <- seq(0, 2^(l + 1))
    i <- seq(0, 2^(l + 2))
    A <- k * (pi/2^(l + 1))
    B <- i * (pi/2^(l + 1))
    A <- A - (pi/2)
    B <- B - pi
    A <- 90 * A/(pi/2)
    B <- 180 * B/pi
    grid <- NULL
    for (i in 1:length(A)) 
        for (j in 1:length(B)) 
            grid <- rbind(grid, c(A[i], B[j]))

    list(grid = grid)
}

"gotredgrid" <-
function (l) 
{
    points <- NULL
    for (k in 0:2^(l + 1)) {
        a <- k * (pi/2^(l + 1))
        if (pi/4 <= a && a <= 3 * pi/4) 
            r <- 0

        if (0 < a && a < pi/4)
            r <- floor((l + 1) - logb(pi * k, 2))

        if (3 * pi/4 < a && a < pi) 
            r <- floor((l + 1) - logb(pi * (2^(l + 1) - k), 2))

        if (a == 0 || a == pi)
            r <- (l - 1)

        rr <- (2^(l + 2) - 2^r)/2^r + 1
        for (i in 0:rr) {
            b <- i * 2^r * pi/2^(l + 1)
            points <- rbind(points, c(a, b))
        }
    }
    grid <- cbind(points[, 1] - (pi/2), points[, 2] - pi)
    grid <- cbind(90 * grid[, 1]/(pi/2), 180 * grid[, 2]/pi)
    
    list(grid = grid)
}
