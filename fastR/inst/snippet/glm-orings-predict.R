predict(orings.model,newdata=data.frame(temp=31)) -> r; r
ilogit(r)                        # inverse logit transformation
predict(orings.model,newdata=data.frame(temp=31),type='response')->p; p
