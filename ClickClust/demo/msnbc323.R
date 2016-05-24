
set.seed(1234)

data(msnbc323)

n <- length(msnbc323)

C <- click.read(msnbc323)

# DATA ANALYSIS USING THE EM ALGORITHM (with beta's)

M1 <- click.EM(X = C$X, y = C$y, K = 1, iter = 10, r = 250, scale.const = 2)
M2 <- click.EM(X = C$X, y = C$y, K = 2, iter = 10, r = 250, scale.const = 2)
M3 <- click.EM(X = C$X, y = C$y, K = 3, iter = 10, r = 250, scale.const = 2)
M4 <- click.EM(X = C$X, y = C$y, K = 4, iter = 10, r = 250, scale.const = 2)
M5 <- click.EM(X = C$X, y = C$y, K = 5, iter = 10, r = 250, scale.const = 2)

# DATA ANALYSIS USING THE EM ALGORITHM (without beta's)

N1 <- click.EM(X = C$X, K = 1, iter = 10, r = 250, scale.const = 2)
N2 <- click.EM(X = C$X, K = 2, iter = 10, r = 250, scale.const = 2)
N3 <- click.EM(X = C$X, K = 3, iter = 10, r = 250, scale.const = 2)
N4 <- click.EM(X = C$X, K = 4, iter = 10, r = 250, scale.const = 2)
N5 <- click.EM(X = C$X, K = 5, iter = 10, r = 250, scale.const = 2)

rbind(c(M1$BIC, M2$BIC, M3$BIC, M4$BIC, M5$BIC),
      c(N1$BIC, N2$BIC, N3$BIC, N4$BIC, N5$BIC))

state.names <- c("frontpage", "news", "tech", "local", "opinion", "on-air",
  "misc", "weather", "msn-news", "health", "living", "business", "msn-sports",
  "sports", "summary", "bbs", "travel")

dev.new(width = 11, height = 11)
click.plot(X = C$X, y = C$y, id = M3$id, col.levels = 10, top.srt = 90,
           colors = c("lightyellow", "red", "darkred"), font.cex = 1.5,
           marg = 4, states = state.names)