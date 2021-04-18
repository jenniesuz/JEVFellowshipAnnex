library(plyr)
library(bbmle)
#******************************likelihood functions*******************

#****************single FOI************
nll.binom <- function(logfoi,dat){ 
  num <- dat$totN - dat$immN                # number individuals positive
  denom <- dat$totN                         # number individuals screened
  
  foi <- exp(logfoi)
  
  ll <- sum(dbinom(dat$immN,prob=(1-exp(-foi*dat$meanAge)),size=dat$totN,log=T))
  return(-ll)
}
#**************two FOIS**************
nll.binom2 <- function(logfoi1,logfoi2,dat){ 
  num <- dat$totN - dat$immN                #number individuals positive
  denom <- dat$totN                         #number individuals screened
  
  dat1 <- dat[(dat$v %in% "1")|(dat$v %in% "2"),]
  dat2 <- dat[dat$v %in% "3",]

  foi1 <- exp(logfoi1)
  foi2 <- exp(logfoi2)
  
  db1 <- sum(dbinom(dat1$immN,prob=(1-exp(-foi1*dat1$meanAge)),size=dat1$totN,log=T))
  db2 <- sum(dbinom(dat2$immN,prob=(1-exp(-foi2*dat2$meanAge)),size=dat2$totN,log=T))

  ll <- sum(db1,db2)
  return(-ll)
}
#***************three FOIS***************
nll.binom3 <- function(logfoi1,logfoi2,logfoi3,dat){ 
  num <- dat$totN - dat$immN                #number individuals positive
  denom <- dat$totN                         #number individuals screened
  
  dat1 <- dat[dat$v %in% "1",]
  dat2 <- dat[dat$v %in% "2",]
  dat3 <- dat[dat$v %in% "3",]
  
  foi1 <- exp(logfoi1)
  foi2 <- exp(logfoi2)
  foi3 <- exp(logfoi3)
  
  db1 <- sum(dbinom(dat1$immN,prob=(1-exp(-foi1*dat1$meanAge)),size=dat1$totN,log=T))
  db2 <- sum(dbinom(dat2$immN,prob=(1-exp(-foi2*dat2$meanAge)),size=dat2$totN,log=T))
  db3 <- sum(dbinom(dat3$immN,prob=(1-exp(-foi3*dat3$meanAge)),size=dat3$totN,log=T))
  
  ll <- sum(db1,db2,db3)
  return(-ll)
}
#********************************************************

#************************************************************
# constant FOI, proportion of hosts aged a, at time t that remain susceptible:
# can be estimated as, lambda=foi
propSusFunc <- function(lambda,age){
  sus <- exp(-lambda*age)
  inf <- 1 - sus
  return(c(sus,inf))
}
propSusFunc.V <- Vectorize(propSusFunc)
#***********************************************************


#*****************simulate data*******************
# Single FOI

# simulateDataFunc <- function(trueFoi,numHosts,ageGroups,maxAge){
# 
#   randomAge <- runif(numHosts,1,maxAge)
#   randomAgeGroups <- cut(randomAge,ageGroups,labels=c(1:7))
# 
#   ages <- cbind.data.frame(randomAge,randomAgeGroups)
#   meanAges <- ddply(ages,.(randomAgeGroups),transform,meanAge=round(mean(randomAge),0))
# 
#   #*******true underlying probability of having had infection
#   trueProps <- propSusFunc.V(trueFoi,meanAges$randomAge)
#   meanAges$trueSusProb <- trueProps[1,]
#   meanAges$trueImmProb <- 1 - trueProps[1,]
# 
#   #****randomly assign to whether have had or not based on true prob****
#   meanAges$randomImm <- sapply(1:length(meanAges[,1]),function(x){
#     temp <- meanAges[x,]
#     return(rbinom(1,1,prob=temp$trueImmProb))
#   })
# 
#   # sum by age 
#   propInf <- ddply(meanAges,.(meanAge),summarise,totN=length(randomImm),immN=sum(randomImm))
#   propInf$propImm <- propInf$immN/propInf$totN
# 
#   mleFit <- mle2(function(par1){nll.binom(par1,propInf)}             
#                  ,start=list(par1=log(0.04)))
#   
#   biasFOI <- trueFoi - exp(coef(mleFit))
# 
# #return(list(meanAges,propInf,biasFOI,estFOI=exp(coef(mleFit))))
#  return(cbind.data.frame(biasFOI,estFOI=exp(coef(mleFit)),numHosts))
# }

#*****************************three FOIs func*******************************

simulateDataFunc <- function(trueFoi1,trueFoi2,trueFoi3,numHosts,ageGroups,maxAge){
  
  randomAge1 <- runif(numHosts,6,maxAge)
  randomAgeGroups1 <- cut(randomAge1,ageGroups,labels=c(1:7))
  ages1 <- cbind.data.frame(randomAge1,randomAgeGroups1)
  meanAges1 <- ddply(ages1,.(randomAgeGroups1),transform,meanAge=round(mean(randomAge1),0))
  
  randomAge2 <- runif(numHosts,6,maxAge)
  randomAgeGroups2 <- cut(randomAge2,ageGroups,labels=c(1:7))
  ages2 <- cbind.data.frame(randomAge2,randomAgeGroups2)
  meanAges2 <- ddply(ages2,.(randomAgeGroups2),transform,meanAge=round(mean(randomAge2),0))
  
  randomAge3 <- runif(numHosts,6,maxAge)
  randomAgeGroups3 <- cut(randomAge3,ageGroups,labels=c(1:7))
  ages3 <- cbind.data.frame(randomAge3,randomAgeGroups3)
  meanAges3 <- ddply(ages3,.(randomAgeGroups3),transform,meanAge=round(mean(randomAge3),0))
  
  #*******true underlying probability of having had infection
  trueProps1 <- propSusFunc.V(trueFoi1,meanAges1$randomAge1)
  meanAges1$trueSusProb <- trueProps1[1,]
  meanAges1$trueImmProb <- 1 - trueProps1[1,]
  
  trueProps2 <- propSusFunc.V(trueFoi2,meanAges2$randomAge2)
  meanAges2$trueSusProb <- trueProps2[1,]
  meanAges2$trueImmProb <- 1 - trueProps2[1,]
  
  trueProps3 <- propSusFunc.V(trueFoi3,meanAges3$randomAge3)
  meanAges3$trueSusProb <- trueProps3[1,]
  meanAges3$trueImmProb <- 1 - trueProps3[1,]
  
  #****randomly assign to whether have had or not based on true prob****
  meanAges1$randomImm <- sapply(1:length(meanAges1[,1]),function(x){
    temp <- meanAges1[x,]
    return(rbinom(1,1,prob=temp$trueImmProb))
  })
  meanAges2$randomImm <- sapply(1:length(meanAges2[,1]),function(x){
    temp <- meanAges2[x,]
    return(rbinom(1,1,prob=temp$trueImmProb))
  })
  meanAges3$randomImm <- sapply(1:length(meanAges3[,1]),function(x){
    temp <- meanAges3[x,]
    return(rbinom(1,1,prob=temp$trueImmProb))
  })
  
  # sum by age 
  propInf1 <- ddply(meanAges1,.(meanAge),summarise,totN=length(randomImm),immN=sum(randomImm))
  propInf1$propImm <- propInf1$immN/propInf1$totN
  propInf2 <- ddply(meanAges2,.(meanAge),summarise,totN=length(randomImm),immN=sum(randomImm))
  propInf2$propImm <- propInf2$immN/propInf2$totN
  propInf3 <- ddply(meanAges3,.(meanAge),summarise,totN=length(randomImm),immN=sum(randomImm))
  propInf3$propImm <- propInf3$immN/propInf3$totN
  
  propInf1$v <- "1"
  propInf2$v <- "2"
  propInf3$v <- "3"
  
  propInf <- rbind.data.frame(propInf1,propInf2,propInf3)
  
  mleFit3 <- mle2(function(par1,par2,par3){nll.binom3(par1,par2,par3,propInf)}             
                  ,start=list(par1=log(0.01),par2=log(0.05),par3=log(0.1)))
  
  mleFit2 <- mle2(function(par1,par2){nll.binom2(par1,par2,propInf)}             
                  ,start=list(par1=log(0.01),par2=log(0.05)))
  
  mleFit <- mle2(function(par1){nll.binom(par1,propInf)}             
                 ,start=list(par1=log(0.05)))
  
  aov1 <- anova(mleFit,mleFit3)
  prob1<-aov1[,5]
  p1 <- if(prob1[2]<=0.05){
    p1 <- 1
  }else{p1<-0}
  
  aov2 <- anova(mleFit2,mleFit3)
  prob2<-aov2[,5]
  p2 <- if(prob2[2]<=0.05){
    p2 <- 1
  }else{p2<-0}
  
  return(c(p1=p1,p2=p2))
}

#*******************************************************************


