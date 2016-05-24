library('maSAE')
message('## ## ## unclustered')
message('## ## non-exhaustive')
message('## load s2 data')
data('s2'); str(s2)
message('## load s1 data')
data('s1'); str(s1)
message('## add sample indicators to s2')
s2$s1 <- s2$s2 <- TRUE
message('## add sample indicators to s1')
s1$s1 <- TRUE
s1$s2 <- FALSE
message('## prepare s1 data')
eval(parse(text=(paste('s1$', setdiff(names(s2), names(s1)), ' <- NA' , sep = ''))))
message('## union s1 and s2 data')
s12 <- rbind(s1, s2)
message('## create object')
saeO <- saObj(data = s12, f = y ~x1 + x2 + x3 | g, s2 = 's2')
message('## small area estimation')
predict(saeO)


message('## ## three-phase')
message('## load s0 data')
data('s0')
message('## s0 has all 3 potential predictors, we keep one only')
s0$x1 <- s0$x3 <- NULL
message('## add sample indicators to s1')
s0$s1 <- s0$s2 <- FALSE
message('## prepare s1 data')
eval(parse(text=(paste('s0$', setdiff(names(s12), names(s0)), ' <- NA' , sep = ''))))
message('## union s12 and s0 data')
s012 <- rbind(s0, s12)
message('## create object')
saeO <- saObj(data = s012,  f = y ~x1 + x2 + x3 | g, s1 = 's1', s2 = 's2')
message('## small area estimation')
predict(saeO)

message('## ## partially exhaustive')
tm <- as.data.frame(tapply(s012$x2 , s012$g, mean));tm$g <- row.names(tm); names(tm) <- c('x2', 'g')
predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g,  s2 = 's2', smallAreaMeans = tm))

message('## ## exhaustive')
message('## true Means as sample means from s12')
preds <- paste('x',1:3, sep='')
tm <- as.data.frame(
    rbind(
        colMeans(subset(s12, g =='a')[, preds])
        , colMeans(subset(s12, g =='b')[, preds])
        )
    ); tm$g=c('a', 'b')
predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g,  s2 = 's2', smallAreaMeans = tm))

message('## ## ## clustered')
message('## ## non-exhaustive')
message('## load s2 data')
data('s2'); str(s2)
message('## load s1 data')
data('s1'); str(s1)
message('## add sample indicators to s2')
s2$s1 <- s2$s2 <- TRUE
message('## add sample indicators to s1')
s1$s1 <- TRUE
s1$s2 <- FALSE
message('## prepare s1 data')
eval(parse(text=(paste('s1$', setdiff(names(s2), names(s1)), ' <- NA' , sep = ''))))
message('## union s1 and s2 data')
s12 <- rbind(s1, s2)
message('## create object')
saeO <- saObj(data = s12, f = y ~x1 + x2 + x3 | g, s2 = 's2', cluster = 'clustid')
message('## small area estimation')
predict(saeO)


message('## ## three-phase')
message('## load s0 data')
data('s0')
message('## s0 has all 3 potential predictors, we keep one only')
s0$x1 <- s0$x3 <- NULL
message('## add sample indicators to s1')
s0$s1 <- s0$s2 <- FALSE
message('## prepare s1 data')
eval(parse(text=(paste('s0$', setdiff(names(s12), names(s0)), ' <- NA' , sep = ''))))
message('## union s12 and s0 data')
s012 <- rbind(s0, s12)
message('## create object')
saeO <- saObj(data = s012,  f = y ~x1 + x2 + x3 | g, s1 = 's1', s2 = 's2', cluster = 'clustid')
message('## small area estimation')
predict(saeO)

message('## ## partially exhaustive')
tm <- as.data.frame(tapply(s012$x2 , s012$g, mean));tm$g <- row.names(tm); names(tm) <- c('x2', 'g')
predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g,  s2 = 's2', cluster = 'clustid', smallAreaMeans = tm))

message('## ## exhaustive')
message('## true Means as sample means from s12')
preds <- paste('x',1:3, sep='')
tm <- as.data.frame(
    rbind(
        colMeans(subset(s12, g =='a')[, preds])
        , colMeans(subset(s12, g =='b')[, preds])
        )
    ); tm$g=c('a', 'b')
predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g,  s2 = 's2', cluster = 'clustid', smallAreaMeans = tm))

