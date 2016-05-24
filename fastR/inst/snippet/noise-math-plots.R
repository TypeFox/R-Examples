noise.interaction <- xyplot(score~noise,groups=group,data=mathnoise,
    type='a',auto.key=list(lines=T,points=F,columns=2))
noise.xy <- xyplot(score~noise,groups=group,data=mathnoise,
                auto.key=list(lines=F,points=T,columns=2))
noise.interaction2 <- xyplot(score~group,groups=noise,data=mathnoise,
    type='a',auto.key=list(lines=T,points=F,columns=2))
noise.xy2 <- xyplot(score~group,groups=noise,data=mathnoise,
                auto.key=list(lines=F,points=T,columns=2))
