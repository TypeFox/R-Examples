tangentLine <-
function(c1, c2, r1, r2, k = 1)
{
dx <- c1[1] - c2[1]
dy <- c1[2] - c2[2]
dr <- r1 - r2

d <- sqrt(dx ^ 2 + dy ^ 2)
X <- dx / d
Y <- dy / d
R <- dr / d

if (abs(R) > 1) return (c(a=NA, b=NA, c=NA))
mR <- sqrt(1 - R ^ 2)

a <- R * X - k * Y * mR
b <- R * Y + k * X * mR
c <- r1 - (a * c1[1] + b * c1[2])
c(a=a, b=b, c=c)
}
