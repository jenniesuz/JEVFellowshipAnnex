library(lhs)
require(sensitivity)
library(ggplot2)
library(reshape)
library(grid)
library(gridExtra)
library(zoo)
library(plyr)
library(parallel)
library(deSolve) 

source("r_odeThreeHostVarMoz.R")             # read in model code
source("r_odeThreeHostVarMozSummaryFuncs.R") # read in functions to summarise model output

#*************parameters to vary**********
eip <- c(1/14,1/5) 
prob2moz <- c(0.01,0.8)
prob2host <- c(0.5,1)
maxMoz <- c(10^3.3,10^4)
birdPref <- c(0.01,0.1)

min <- c(eip[1],prob2moz[1],prob2host[1],maxMoz[1],birdPref[1])
max <- c(eip[2],prob2moz[2],prob2host[2],maxMoz[2],birdPref[2])

params <- c("eip"
            ,"prob2moz" 
            ,"prob2host"
            ,"maxMoz" 
            ,"birdPref"
            )

params <- cbind.data.frame(params,min,max)
# # select random sets of parameter values within parameter value ranges
r <- randomLHS(500,length(params[,1]))
parmVals <- lapply(1:length(params[,1]),function(x){
  temp <- params[x,]
  randomSample <- runif(r[,x],min=temp$min,max=temp$max)
  
})
parmVals <- do.call(cbind.data.frame,parmVals)
names(parmVals) <- params$params

randSnd <- parmVals
randSnd$run <- 1:500
#*************************************
#********cattle:bird scenarios********
numHosts <- c(1100)         # number of hosts/ village
numCattle <- c(numHosts/100*5
               ,numHosts/100*10
               ,numHosts/100*15
               ,numHosts/100*20
               ,numHosts/100*25
               ,numHosts/100*30
               ,numHosts/100*35
               ,numHosts/100*40
               ,numHosts/100*45
               ,numHosts/100*50) 
numBirds <- numHosts - numCattle  # then assign numbers of birds
numBirds/(numCattle+numBirds)
#*****************************************************
randSnd$numBirds <- numBirds[4]
randSnd$numCattle <- numCattle[4]
#*************run scenario assuming c 80% birds*****************
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

sc1Sims <- parRapply(cl,randSnd,function(x){
  library(plyr)
  library(zoo)
  library(binom)
  source("r_odeThreeHostVarMoz.R")
  source("r_odeThreeHostVarMozSummaryFuncs.R")
  x<-as.numeric(x)
  initial <- c(Sv=2000        # susceptible vectors
               ,Ev1=0         # fed and exposed/ infected vectors
               ,Ev2=0
               ,Ev3=0
               ,Iv=1         # infectious vectors
               ,Sh1=x[7]     # susceptible hosts 
               ,Ih1=0        # infected hosts
               ,Rh1=0        # recovered hosts
               ,Sh2=1
               ,Ih2=0
               ,Rh2=0
               ,Sh3=x[8]
               ,Ih3=0
               ,Rh3=0
  ) 
  yrs <- 4                   # How long to run model for
  times <- seq(1,365*yrs,1)  # Timesteps for reporting
  
  simData <-  as.data.frame(lsoda(y=initial
                                  , times=times
                                  , func=mod
                                  , parms=params(
                                    years=yrs
                                    # mosquito population dynamics parameters
                                    , maxRecruit = x[4]
                                    ,baselineRecruit = 0.2    #
                                    ,interval = 1             #
                                    ,firstRecruit = 0.4       #
                                    ,firstMean = 130          #
                                    ,firstSd = 5              #
                                    ,secondMean = 150         #
                                    ,secondSd = 5   
                                    ,mu_v=1/20
                                    ,sigma_v = x[1]           # EIP - relatively long
                                    # transmission parameters to change
                                    ,p_h1v = x[2]         # prob host to vector transmission
                                    ,p_h2v = x[2]
                                    ,p_h3v = 0
                                    ,mu_h1 = 1/90         # host mortality rate - birds
                                    ,mu_h2 = 1/365        # pigs
                                    ,mu_h3 = 1/(365*5)    # cattle
                                    ,beta_h1  = 1/90      # host1 birth rate
                                    ,beta_h2 = 1/365
                                    ,beta_h3 = 1/(365*5)
                                    ,pref_h1 = x[5]        # preference of vector for host 1 (birds)
                                    ,pref_h2 = 0         # pigs - currently no pigs
                                    ,pref_h3 = 1-x[5]        # cattle
                                    ,p_vh1 = x[3]          # prob vector to host1 transmission
                                    ,p_vh2 = x[3]
                                    ,p_vh3 = x[3]
                                    ,phi_h1 = 0.2         # host recovery rate
                                    ,phi_h2 = 0.2
                                    ,phi_h3 = 0.2
                                  )
  ))
  
  simData$run <- x[6]
  
  # take the last year
  oneYr <- simData[1095:(1095+364),]   #*NB this is hard coded - need to adjust so can simulate different times
  
  oneYr$newTimes <- c(1:length(oneYr[,1]))
  oneYr$date <- as.Date(oneYr$newTimes, origin = "2019-01-01")      # assign dates so easier to view and summarise
  oneYr$Yrmon <- as.yearmon(oneYr$date)
  
  oneYr$totMozInf <- oneYr$Iv+oneYr$Ev1+oneYr$Ev2+oneYr$Ev3         # infected mosquitoes over time
  oneYr$totMoz <- oneYr$Iv+oneYr$Ev1+oneYr$Ev2+oneYr$Ev3+oneYr$Sv   # total mosquitoes over time
  
  
  byMonthMozPrev <- ddply(oneYr,.(Yrmon),summarise,meanP=(sum(totMozInf)/sum(totMoz))*1000,totMoz=mean(totMoz)) # rough monthly prevalence
  
  hosts <- sum(x[2],x[3])
 
  
  byMonthMozPrev$moz2h <- byMonthMozPrev$totMoz/hosts               # monthly mosquito to host ratio
  byMonthMozPrev$run <- x[6]                                        # assign unique number for each combination of params
  
  summaryFOIc <- ddply(oneYr,.(Yrmon),summarise,monFOI=sum(FOIc))   # return force of infection by month
  mean(summaryFOIc[,2])  # average FOI
  
  byMonthMozPrev$FOIc <- summaryFOIc$monFOI                         # add to table to be returned
  return(byMonthMozPrev)
  
})

stopCluster(cl)

#*****************************************************************
test <- do.call(rbind,sc1Sims)
test$run <- as.factor(test$run)

dat <- cbind.data.frame(randSnd,test)


dat<-dat[,-6]


dat <- melt(dat,measure.vars=c("eip"
                               ,"prob2moz" 
                               ,"prob2host"
                               ,"maxMoz" 
                               ,"birdPref"
))



dat$variable2 <- factor(dat$variable
                        , levels = c("eip"
                                     ,"prob2moz" 
                                     ,"prob2host"
                                     ,"maxMoz" 
                                     ,"birdPref"
                        )
                        , ordered = TRUE
                        , labels=c("Extrinsic incubation period"
                                   ,"Probability of host to vector transmission"
                                   ,"Probability of mosquito to host transmission"
                                   ,"Maximum mosquito density"
                                   ,"Preference for birds"
                        ))


ggplot(dat) +
  geom_point(aes(x=value,y=meanP,col=as.factor(Yrmon)),size=0.01,alpha=0.2)  +
  labs( y= "Mean mosquitoes infected/1000"
        , x="Parameter value") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'grey')
        ,text=element_text(size=9)
        ,panel.spacing = unit(1, "lines")
        #,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.9,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=9)
        ,legend.position =c(0.9,0.1)
        ,legend.title = element_text(size=12)
        ,strip.background = element_rect(fill="white",color="white")
  ) +
  facet_wrap(~variable,scales="free_x",labeller = label_parsed)


ggplot(dat) +
  geom_point(aes(x=value,y=FOIc,col=as.factor(Yrmon)),size=0.01,alpha=0.2)  +
  labs( y= "Monthly FOI"
        , x="Parameter value") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'grey')
        ,text=element_text(size=9)
        ,panel.spacing = unit(1, "lines")
        #,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.9,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=9)
        ,legend.position =c(0.9,0.1)
        ,legend.title = element_text(size=12)
        ,strip.background = element_rect(fill="white",color="white")
  ) +
  facet_wrap(~variable,scales="free_x",labeller = label_parsed)


ggplot(dat) +
  geom_point(aes(y=meanP,x=Yrmon),size=0.01,alpha=0.2)  +
  labs( y= "Mean mosquitoes infected/1000"
        , x="Parameter value") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'grey')
        ,text=element_text(size=9)
        ,panel.spacing = unit(1, "lines")
        #,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.9,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=9)
        ,legend.position =c(0.9,0.1)
        ,legend.title = element_text(size=12)
        ,strip.background = element_rect(fill="white",color="white")
  ) +
  facet_wrap(~variable,scales="free_x",labeller = label_parsed)

#average the FOI over the year and remove those for which average is > 0.5
dat1 <- ddply(dat,.(run,variable),transform,meanFOI=mean(FOIc))
dat2 <- dat1[(dat1$meanFOI<=5)&(dat1$meanFOI>=0.2),]

# then want to know min and max MLEP
dat3 <- ddply(dat2,.(run,variable),summarise,minMLE=min(meanP),maxMLE=max(meanP))
round(dat3[,c(3,4)],2)

# check this seems reasonable

ggplot(dat2) +
  geom_point(aes(y=meanP,x=Yrmon),size=0.01,alpha=0.2)  +
  labs( y= "Mean monthly mosquitoes infected/1000"
        , x="Parameter value") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'grey')
        ,text=element_text(size=9)
        ,panel.spacing = unit(1, "lines")
        #,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.9,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=9)
        ,legend.position =c(0.9,0.1)
        ,legend.title = element_text(size=12)
        ,strip.background = element_rect(fill="white",color="white")
  ) +
  facet_wrap(~variable,scales="free_x",labeller = label_parsed)

# Mean monthly MLE p - c. 1 to a peak of 20
ddply(dat2,.(Yrmon),summarise,meanMLE=round(mean(meanP),2))
