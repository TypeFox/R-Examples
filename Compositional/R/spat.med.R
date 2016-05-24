################################
#### Spatial median 
#### Tsagris Michail 10/2014  
#### References: Jyrki Mottonen, Klaus Nordhausen and Hannu Oja (2010)
#### Asymptotic theory of the spatial median 
#### In Nonparametrics and Robustness in Modern Statistical Inference and Time Series
#### Analysis: A Festschrift in honor of Professor Jana Jureckova 
#### On computation of spatial median for robust data mining (2005)
#### T. Karkkaminen and S. Ayramo
#### Evolutionary and Deterministic Methods for Design, Optimization
#### and Control with Applications to Industrial and Societal Problems 
#### EUROGEN 2005
#### R. Schilling, W.Haase, J. Periaux, H. Baier, G. Bugeda (Eds)
#### FLM, Munich, 2005
#### http://users.jyu.fi/~samiayr/pdf/ayramo_eurogen05.pdf
#### mtsagris@yahoo.gr
################################

spat.med <- function(x) {
  x <- as.matrix(x)
  p <- ncol(x)
  s <- diag(p)
  u <- matrix(nrow = 10000, ncol = p)
  u[1, ] <- apply(x, 2,median)
  ww <- sqrt( mahalanobis(x, u[1, ], s ) )
  u[2, ] <- colSums(x / ww) / sum(1 / ww)
  i <- 2
  while ( sum( abs(u[i, ] - u[i - 1, ]) ) > 1e-10 ) {
    i <- i +1
    ww <- sqrt( mahalanobis(x, u[i - 1, ], s ) )
    u[i, ] <- colSums(x / ww) / sum(1 / ww)
  }
  u[i, ]
} 
