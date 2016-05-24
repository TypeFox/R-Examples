library("partykit")
set.seed(1)

dat <- data.frame(v1 = as.double(1:100))

sv1 <- partysplit(as.integer(1), breaks = as.double(50))
character_split(sv1, dat)
stopifnot(all(kidids_split(sv1, dat) == ((dat$v1 > 50) + 1)))

sv1 <- partysplit(as.integer(1), breaks = as.double(50), 
                index = as.integer(c(2, 1)))
character_split(sv1, dat)
stopifnot(all(kidids_split(sv1, dat) == ((dat$v1 <= 50) + 1)))

sv1 <- partysplit(as.integer(1), breaks = as.double(50), right = FALSE)
character_split(sv1, dat)
stopifnot(all(kidids_split(sv1, dat) == ((dat$v1 >= 50) + 1)))

sv1 <- partysplit(as.integer(1), breaks = as.double(50), 
                index = as.integer(c(2, 1)), right = FALSE)
character_split(sv1, dat)
stopifnot(all(kidids_split(sv1, dat) == ((dat$v1 < 50) + 1)))

sv1 <- partysplit(as.integer(1), breaks = as.double(c(25, 75)))
character_split(sv1, dat)
stopifnot(all(kidids_split(sv1, dat) == 
              as.integer(cut(dat$v1, c(-Inf, 25, 75, Inf)))))

sv1 <- partysplit(as.integer(1), breaks = as.double(c(25, 75)), right = FALSE)
character_split(sv1, dat)
stopifnot(all(kidids_split(sv1, dat) == 
              as.integer(cut(dat$v1, c(-Inf, c(25, 75), Inf), right = FALSE))))

sv1 <- partysplit(as.integer(1), breaks = as.double(c(25, 75)), 
                index = as.integer(3:1), right = FALSE)
character_split(sv1, dat)
stopifnot(all(kidids_split(sv1, dat) == 
              (3:1)[as.integer(cut(dat$v1, c(-Inf, c(25, 75), Inf), right = FALSE))]))


dat$v2 <- gl(4, 25)

sv2 <- partysplit(as.integer(2), index = as.integer(c(1, 2, 1, 2)))
character_split(sv2, dat)
kidids_split(sv2, dat)

sv2 <- partysplit(as.integer(2), breaks = as.integer(c(1, 3)))
character_split(sv2, dat)
kidids_split(sv2, dat)


dat <- data.frame(x = gl(3, 30, labels = LETTERS[1:3]), y = rnorm(90), 
                  z = gl(9, 10, labels = LETTERS[1:9], ordered = TRUE))
csp <- partysplit(as.integer(1), index = as.integer(c(1, 2, 1)))
kidids_split(csp, dat)
#kidids_node(list(csp), dat)

nsp <- partysplit(as.integer(2), breaks = c(-1, 0, 1), index = as.integer(c(1, 2, 1, 3)))
kidids_split(nsp, dat)

osp <- partysplit(as.integer(3), breaks = as.integer(c(3, 6)), index = as.integer(c(2, 1, 2)))
kidids_split(osp, dat)

nadat <- dat
nadat$x[1:10] <- NA
nadat$y[11:20] <- NA
#kidids_node(list(csp, nsp, osp), nadat)

character_split(csp, dat)
character_split(nsp, dat)
character_split(osp, dat)

