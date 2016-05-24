caFactorialDesign<-function(data, type="null", cards=NA)
{
        set.seed(123)
        num<-data.frame(data.matrix(data))
        vars.number<-length(num)
        levels.number<-0
        for (i in 1:length(num)) levels.number<-levels.number+max(num[i])
        ca.number<-levels.number-vars.number+1
        aca.number=3*(levels.number-vars.number+1)-levels.number
        profiles.number<-1
        for (i in 1:length(num)) profiles.number<-profiles.number*max(num[i])
        if (type=="full") return(data)
        if (type=="null") temp.design<-optFederov(~., data)
        if (type=="fractional") 
        {
                if (is.na(cards)==TRUE) temp.design<-optFederov(~.,data, approximate=FALSE, nullify=1) 
                else temp.design<-optFederov(~., data, nTrials=cards)
        }
        if (type=="ca") temp.design<-optFederov(~., data, nTrials=ca.number)
        if (type=="aca") temp.design<-optFederov(~., data, nTrials=aca.number)
        if (type=="orthogonal")
        {
                for (i in ca.number: profiles.number)
                {
                        temp.design<-optFederov(~., data, nTrials=i, approximate=FALSE, nRepeats=50)
                        test.design<-temp.design$design
                        if (det(cor(data.matrix(test.design)))==1) break
                }
        }
        return(temp.design$design)
}
