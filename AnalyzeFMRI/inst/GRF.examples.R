
library(fmri, lib.loc = "./tmp/")
## GRF simulations and checking of Bonferroni and GRF thresholds

ksize <- 9
d <- c(64, 64, 21)
FWHM <- 9
sigma <- diag(FWHM^2, 3) / (8 * log(2))
voxdim <- c(2, 2, 4)

filtermat <- GaussSmoothKernel(voxdim, ksize, sigma)
image(filtermat[(ksize + 1) / 2, , ])

mask <- f.read.analyze.volume(file = "/users/marchini/fmri/fmri_data/jm2.null.mc.mask.img")[, , , 1]
num.vox <- sum(mask)

ans <- 1:10000
for(i in 1:10000) {
    print(i)
    ans[i] <- Sim.3D.GRF(d = d, voxdim = voxdim, sigma = sigma, ksize = ksize, mask = mask, type = "max")$max
}

thr.bonf <- Threshold.Bonferroni(p.val = 0.05, n = num.vox)
thr.grf <- Threshold.GRF(p.val = 0.05, sigma = sigma, voxdim = voxdim, num.vox = num.vox)

thr.bonf
thr.grf

sum(ans > thr.bonf)
sum(ans > thr.grf)

a <- sort(ans)
p <- 1:10000 / 10000
plot(a, p, typ = "l")

p1 <- 1 - EC.3D(a, sigma = sigma, voxdim = voxdim, num.vox = sum(mask))
lines(a, p1, col = 2)



