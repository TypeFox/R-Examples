theta <- acos(runif(1, -1, 1))
phi <- runif(1, -pi, pi)
u<- c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
r<-rvmises(1)

context("Conversions")
expect_equal(as.Q4(as.SO3(u,r)),as.Q4(u,r))
expect_equal(as.SO3(as.Q4(u,r)),as.SO3(u,r))
expect_true(is.SO3(as.SO3(u,r)))
expect_true(is.Q4(as.Q4(u,r)))

#context("Project to SO3")
#Rs<-ruars(5,rcayley)
#Rs<-rbind(Rs,rnorm(9))

#txt<-"Row(s)"
#txt<-paste(txt,6)
#txt<-paste(txt,"was(were) not in SO(3).")

#expect_message(as.SO3(Rs),cat(txt,"\n"))
