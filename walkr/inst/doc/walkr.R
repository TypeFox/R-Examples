## ----example0, eval = TRUE, echo = FALSE, fig.cap = "We use the two sampling algorithms, hit-and-run and Dikin walk, on the 2D, 5D, and 10D-simplex. We show the histograms for the first parameter, $x_1$, because the distribution for every parameter should be the same (there is nothing special about $x_1$). Every histogram contains 10000 sampled points. The 2D-simplex is the line segment described by $x_1+x_2=1$ and $x_i$ greater equal to $0$. For the 2D-simplex, we see that hit-and-run samples uniformly across [0,1], while Dikin walk concentrates in regions away from the center. For higher dimensions (5D and 10D histograms), consider the 3D-simplex analogy. The 3D-simplex (Figure 1) is a triangle in three dimensional space. If we look at the distribution for $x_1$, that is equivalent to projecting this triangle onto the $x_1$ axis. As the samples are drawn uniformly from the 3D-simplex, there are more points near $0$ than near $1$. Therefore, we see this downward sloping distribution. In the 5D and 10D-simplex cases, we see that hit-and-run samples uniformly, while Dikin again concentrates in regions away from the edges.", message= FALSE, cache = TRUE, fig.show = 'hold', out.width = '1\\linewidth', out.height = '12cm', fig.align = 'center'----

library(walkr)

set.seed(314)
## initialize matrix

result <- matrix(1, ncol = 6, nrow = 10000)

## 2D simplex
A <- matrix(1, ncol = 2)
b <- 1 

## hitandrun and dikin
result[,1] <- walkr(A = A, b = b, points = 10000, thin = 1, burn = 0, method = "hit-and-run")[1,]
result[,2] <- walkr(A = A, b = b, points = 10000, thin = 1, burn = 0, method = "dikin")[1,]

## 5D simplex
A <- matrix(1, ncol = 5)
b <- 1 

## hitandrun and dikin
result[,3] <- walkr(A = A, b = b, points = 10000, thin = 1, burn = 0, method = "hit-and-run")[1,]
result[,4] <- walkr(A = A, b = b, points = 10000, thin = 1, burn = 0, method = "dikin")[1,]

## 10D simplex
A <- matrix(1, ncol = 10)
b <- 1 

## hitandrun and dikin
result[,5] <- walkr(A = A, b = b, points = 10000, thin = 1, burn = 0, method = "hit-and-run")[1,]
result[,6] <- walkr(A = A, b = b, points = 10000, thin = 1, burn = 0, method = "dikin")[1,]


df <- as.data.frame(result)
colnames(df) <- c("har2", 
                  "dikin2",
                  "har5", 
                  "dikin5",
                  "har10", 
                  "dikin10")
library(grid)
library(gridExtra)

m1 <- ggplot(df, aes(x=har2)) 
m2 <- ggplot(df, aes(x=dikin2))
m3 <- ggplot(df, aes(x=har5))
m4 <- ggplot(df, aes(x=dikin5))
m5 <- ggplot(df, aes(x=har10))
m6 <- ggplot(df, aes(x=dikin10))

m1 <- m1 + geom_histogram(binwidth = 0.02, colour = "darkgreen", fill = "black") + 
      xlab(expression(paste("Value of ", x[1]))) + 
      ggtitle("Hit-and-run: 2D-simplex") + scale_y_continuous(limits = c(0,500))
m2 <- m2 + geom_histogram(binwidth = 0.02, colour = "darkgreen", fill = "black") + 
      xlab(expression(paste("Value of ", x[1]))) + 
      ggtitle("Dikin walk: 2D-simplex") + scale_y_continuous(limits = c(0,500))
m3 <- m3 + geom_histogram(binwidth = 0.02, colour = "darkgreen", fill = "black") + 
      xlab(expression(paste("Value of ", x[1]))) + 
      ggtitle("Hit-and-run: 5D-simplex") + scale_y_continuous(limits = c(0,1000))
m4 <- m4 + geom_histogram(binwidth = 0.02, colour = "darkgreen", fill = "black") + 
      xlab(expression(paste("Value of ", x[1]))) + 
      ggtitle("Dikin walk: 5D-simplex") + scale_y_continuous(limits = c(0,1000))
m5 <- m5 + geom_histogram(binwidth = 0.02, colour = "darkgreen", fill = "black") + 
      xlab(expression(paste("Value of ", x[1]))) + 
      ggtitle("Hit-and-run: 10D-simplex") + scale_y_continuous(limits = c(0,2000))
m6 <- m6+ geom_histogram(binwidth = 0.02, colour = "darkgreen", fill = "black") + 
      xlab(expression(paste("Value of ", x[1]))) + 
      ggtitle("Dikin walk: 10D-simplex") + scale_y_continuous(limits = c(0,2000))

## plot it 3 by 2 

grid.arrange(m1,m2,m3,m4,m5,m6, ncol = 2)


## ----example1, eval = TRUE, cache = TRUE---------------------------------
A <- matrix(1, ncol = 3)
b <- 1
sampled_points <- walkr(A = A, b = b, points = 1000, 
                        method = "hit-and-run", chains = 5, ret.format = "matrix")

## ----example3, eval = TRUE, cache = TRUE---------------------------------
A <- matrix(sample(c(0,1,2), 40, replace = TRUE), ncol = 20)
b <- c(0.5, 0.3)
sampled_points <- walkr(A = A, b = b, points = 10000, chains = 5, 
                        method = "hit-and-run", ret.format = "list")


## ----example4, eval = TRUE, cache = TRUE---------------------------------
sampled_points <- walkr(A = A, b = b, points = 1000, chains = 5, thin = 500, 
                        method = "hit-and-run", ret.format = "list")         

## ----example5, eval = TRUE, cache = TRUE---------------------------------
sampled_points <- walkr(A = A, b = b, points = 1000, chains = 5, thin = 10,  
                        method = "dikin", ret.format = "list")

