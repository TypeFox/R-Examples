require("dissUtils")

make.example.data <- function(N = 250,
                              n.groups = 5,
                              scale.k = 12,
                              p.count = 20){

    fake.data <- data.frame(ID=factor(sample(1:n.groups, N, replace = T)));

    fake.data$radius <- runif(nrow(fake.data), 0, 1);

    fake.data$angle <- as.numeric(fake.data$ID) * pi / n.groups;

    fake.data$X <- scale.k*fake.data$radius * cos(fake.data$angle) + rnorm(N);

    fake.data$Y <- scale.k*fake.data$radius * sin(fake.data$angle) + rnorm(N);

    return(invisible(fake.data));
}

make.density.data <- function(Y, ids, p){

    Y <- as.matrix(Y);

    ids <- as.factor(ids);

    n.groups <- nlevels(ids);

    require(dissUtils);

    hoops <- data.frame(ID = factor(levels(ids)),
                        p = rep(p,n.groups),
                        delta = numeric(n.groups),
                        ix = integer(n.groups),
                        X = numeric(n.groups),
                        Y = numeric(n.groups));

    dense <- groupwise.density(Y, ids, p.neighbors = p);

    tmp <- apply(scale(dense,
                       rep(0,ncol(dense)),
                       apply(dense,2,sum)),
                 1,
                 mean);

    mix <- match(max(tmp),tmp);

    all.ix <- 1:nrow(Y);

    for(j in 1:nlevels(ids)){

        lvl <- levels(ids)[j];

        f <- ids == lvl;

        tmp <- diss(matrix(Y[mix,],nrow = 1),
                    Y[f,]);

        hoops$delta[j] <- min(tmp);

        ix <- match(hoops$delta[j], tmp);

        hoops$ix[j] <- all.ix[f][ix];

        hoops[j,5:6] <- Y[f,][ix,];
    }

    return(invisible(hoops));
}

extract <- function(X,nom){as.numeric(t(sapply(X, function(x){x[[nom]]})))};

whole.hoops <- list();
whole.data <- list();

for(lab in c(letters,LETTERS)){

    fake.data <- make.example.data();

    p.count <- nlevels(fake.data$ID);

    hoops <- make.density.data(fake.data[,4:5],
                               fake.data$ID,
                               1);

    for(i in (p.count-1):1){

        temp <- make.density.data(fake.data[,4:5], fake.data$ID, i / p.count);

        hoops <- rbind(temp, hoops);
    }

    whole.hoops[[lab]] <- hoops;

    whole.data[[lab]] <- fake.data;
}


all.data <- data.frame(ID = factor(extract(whole.data, "ID")));

for(n in names(whole.data$A)[-1]){

    all.data[[n]] <- extract(whole.data, n);
}

all.data$ID <- factor(all.data$ID);

all.densest <- data.frame(ID = factor(extract(whole.hoops, "ID")));

for(n in names(whole.hoops$A)[-1]){

    all.densest[[n]] <- extract(whole.hoops, n);
}

temp <- aggregate(all.densest$delta, as.list(all.densest[,1:2]), mean);

temp$sd <- aggregate(all.densest$delta, as.list(all.densest[,1:2]), sd)$x;

temp$N <- aggregate(all.densest$delta, as.list(all.densest[,1:2]), length)$x;

temp$lo <- temp$x - temp$sd / sqrt(temp$N);

temp$hi <- temp$x + temp$sd / sqrt(temp$N);


dev.new(title="all scatter");

par(mar = c(10,10,2,2)/3,
    mgp = 2 * c(3,1,0)/3,
    bg = "black",
    fg = "white", col.axis = "white", col.lab = "white");

plot(Y ~ X, all.data, pch ='.', col = rainbow(5,alpha = 2/3)[ID], asp = 1);

points(Y~X, all.densest, cex = 1/2, col = rainbow(5)[ID]);

dev.new(title="all trends")

par(mar = c(10,10,2,2)/3, mgp = 2 * c(3,1,0)/3)

plot(1,10, type = 'n',
     xlab = "proportion of neighbors examined",
     ylab = "distance to group densest point",
     xlim = c(0,1),
     ylim = range(temp[,6:7]),
     log = 'y');

for(j in 1:5){

    f <- temp$ID == j;

    segments(temp$p[f], temp$lo[f],
             temp$p[f], temp$hi[f],
             col = rainbow(5)[j]);

    points(x ~ p, temp[f,], col = rainbow(5)[j], type = 'o');
}
