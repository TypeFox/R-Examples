ghpoints3 <-
function(K) {
    tmp = ghpoints(K);
    points = tmp[,1];
    weights = tmp[,2];
    points1 = rep(points,each=K^2);
    w1 = rep(weights,each=K^2);
    points2 = rep(rep(points,each=K),K);
    w2 = rep(rep(weights,each=K),K);
    points3 = rep(points,K^2);
    w3 = rep(weights,K^2);
    Points = rbind(points1, points2, points3);
    W = w1*w2*w3;  
    return(list(Points=Points, Weights=W));
}
