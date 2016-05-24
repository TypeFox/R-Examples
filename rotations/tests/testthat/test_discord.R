Rs <- ruars(20,rcayley)
Qs <- as.Q4(Rs)

so3Hi <- discord(Rs, type='i')
qsHi <- discord(Qs, type='i')
so3HiE <- discord(Rs, type='e')
qsHiE <- discord(Qs, type='e')

context("Discord Invariance")
expect_equal(so3Hi, qsHi)
expect_equal(so3HiE, qsHiE)

expect_error(discord(Qs), "type needs to be one of 'intrinsic' or 'extrinsic'")
expect_error(discord(Rs), "type needs to be one of 'intrinsic' or 'extrinsic'")
