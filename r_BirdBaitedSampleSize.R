library(plyr)
library(matrixStats)
library(MASS)
# variance from light trap sampling in households - at the household level over three nights
# from village 4 Cx t minimum 1 maximum 158, mean 51, variance 2627, sd 51

# if mean 50 and var 2600, neg bin distribution: var = mu + mu^2/ size
# var - mu = mu^2/size
# size = mu^2/(var-mu)
50^2 / (2600 - 50) # 0.98
rnbinom(100,mu=50,size=0.98)

# bird-baited traps 'true' differences
# scenario 1 in rice fields
riceBird1 <- c(1:100)   # potential number of mosquitoes caught per day/ evening per trap
riceBird2 <- riceBird1*2
riceBird3 <- riceBird1/2

# scenario 2 in households
houseBird1 <- riceBird1*2
houseBird2 <- riceBird2*2
houseBird3 <- riceBird3*2





# lets say three trials each of each bird x location

# function to simulate mosquito catches from three birds at two locations
sampleFunc <- function(nTrials   # number of samples for each bird/ location
                       ,bird1mean # mean number expected from bird 1 in location 1
                       ,bird2p    # multiplier for bird 2 in location 1
                       ,bird3p    # multiplier for bird 3 in location 1
                       ,houseMult){  # multiplier for all birds in location 2
  
  riceBird1 <- rnbinom(nTrials,mu=bird1mean,size=bird1mean^2/(bird1mean^2 - bird1mean)  )
  riceBird2 <- rnbinom(nTrials,mu=bird1mean*bird2p,size=(bird1mean*bird2p)^2/ ( (bird1mean*bird2p)^2 - (bird1mean*bird2p)))
  riceBird3 <- rnbinom(nTrials,mu=bird1mean*bird3p,size=(bird1mean*bird3p)^2/ ( (bird1mean*bird3p)^2 - (bird1mean*bird3p)) )
  houseBird1 <- rnbinom(nTrials,mu=bird1mean*houseMult,size=(bird1mean*houseMult)^2/ ( (bird1mean*houseMult)^2 - (bird1mean*houseMult)) )
  houseBird2 <- rnbinom(nTrials,mu=bird1mean*bird2p*houseMult,size=(bird1mean*bird2p*houseMult)^2/ ( (bird1mean*bird2p*houseMult)^2 - (bird1mean*bird2p*houseMult)) )
  houseBird3 <- rnbinom(nTrials,mu=bird1mean*bird3p*houseMult,size= (bird1mean*bird3p*houseMult)^2/ ( (bird1mean*bird3p*houseMult)^2 - (bird1mean*bird3p*houseMult)))
  
counts <- c(riceBird1,riceBird2,riceBird3,houseBird1,houseBird2,houseBird3)
location <- c(rep("rice",nTrials*3),rep("house",nTrials*3))
bird <- c(rep("Bird1",nTrials),rep("Bird2",nTrials),rep("Bird3",nTrials),rep("Bird1",nTrials),rep("Bird2",nTrials),rep("Bird3",nTrials))
dat <- cbind.data.frame(counts,location,bird)
return(dat)
}




# now use the function to generate some data each for increasing number of trials for bird x location
# fit a glm to the data

t1 <- 100
t2 <- 2
t3 <- 0.5
t4 <- 1.5

nSims <- 500

# here 1 to 100 trials, mean of 100 for bird 1 in location 1, multipliers of 2, 0.5 & 2
mu100 <- lapply(c(10,20,30,40,50,60,70,80,90,100),function(x){
  
  reps <- rdply(nSims,sampleFunc(x,t1,t2,t3,t4) )
  #randomSample <- sampleFunc(x,100,2,0.5,2) 
  modCoefs <- sapply(unique(reps$.n),function(y){
    repsSub <- reps[reps$.n %in% y,]
    mod <- glm.nb(counts ~ bird+location,data=repsSub)
    coefs <- exp(coefficients(mod))
    
    s <- summary(mod)
    pVals <- as.numeric(coefficients(s)[,4]) < 0.05
    if (sum(pVals) == 4) {
      pVals <- 1
    }else{
      pVals <- 0
    }
    
    strrs <- broom::tidy(mod)$std.error
    ci <- exp(1.96 * strrs)

    return(c(coefs,ci,pVals))
  })
  
  trueVals <- c(t1,t2,t3,1/t4)
  modEst <- modCoefs[1:4,]
  modCi <- modCoefs[5:8,]
  modUpper <- modEst + modCi
  modLower <- modEst - modCi
  
  rMeans <- rowMeans(modEst)
  
  # bias
  bias <- rMeans - trueVals
  # variance
  rVar <- rowVars(modEst)
  # CI width
  
  # mse
  mse <- bias^2 + rVar
  # coverage
  cov1 <- (trueVals[1] < modUpper[1,]) & (trueVals[1] > modLower[1,])
  cov2 <- (trueVals[2] < modUpper[2,]) & (trueVals[2] > modLower[2,])
  cov3 <- (trueVals[3] < modUpper[3,]) & (trueVals[3] > modLower[3,])
  cov4 <- (trueVals[4] < modUpper[4,]) & (trueVals[4] > modLower[4,])
  
  cov1 <- length(cov1[cov1==T]) / length(cov1)
  cov2 <- length(cov2[cov2==T]) / length(cov2)
  cov3 <- length(cov3[cov3==T]) / length(cov3)
  cov4 <- length(cov4[cov4==T]) / length(cov4)
  
  cov <- c(cov1,cov2,cov3,cov4)
  
  power <- (sum(modCoefs[9,])/length(modCoefs[9,]) )*100
  
  summary <- cbind.data.frame(trueVals,power,rMeans,bias,rVar,mse,cov,nTrials=rep(x,4),coef=c("intercept","Bird 2 coefficient","Bird 3 coefficient","Location coefficient"))
  
  return(summary)
})


cs1 <- do.call(rbind.data.frame,mu100)   # this returns the fitted model coefficients for each experiment
# differing in number of trials

# bias - expected difference between the estimate and the true value
# of the parameter
# variance - variability of point estimates around their mean value
# confidence interval width - width of CI in absolute terms or as a proportion of the estimated value
# mean squared error -  bias^2 + variance. Total variation around the true value
# Coverage - proportion of simulations in which CI actually include the true value

#cs1 <- cs[!cs$coef %in% "intercept",]

library(ggplot2)
library(grid)
library(gridExtra)

ggplot(cs1) +
  geom_point(aes(x=nTrials,y=power))

biasp <- ggplot(cs1) +
  geom_point(aes(x=nTrials,y=bias)) + 
  theme_set(theme_bw()) +
  labs(x=" ",y="Bias") +
  theme(panel.border=element_blank()
        ,strip.background = element_rect(colour=NA, fill=NA)
        ,axis.line = element_line(color = 'grey')
        ,text=element_text(size=11)
        ,panel.spacing = unit(1, "lines")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.6,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,legend.position =c(0.89,0.7)
        ,legend.title = element_text(size=10)
  ) +
  facet_wrap(~coef,scales="free")

varp <- ggplot(cs1) +
  geom_point(aes(x=nTrials,y=rVar)) +
  labs(x=" ",y="Variance") +
  theme_set(theme_bw()) +
  theme(panel.border=element_blank()
        ,strip.background = element_rect(colour=NA, fill=NA)
        ,axis.line = element_line(color = 'grey')
        ,text=element_text(size=11)
        ,panel.spacing = unit(1, "lines")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.6,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,legend.position =c(0.89,0.7)
        ,legend.title = element_text(size=10)
  ) +
  facet_wrap(~coef,scales="free")

msep <- ggplot(cs1) +
  geom_point(aes(x=nTrials,y=mse)) + 
  labs(x="Sample size",y="Mean squared error") +
  theme_set(theme_bw()) +
  theme(panel.border=element_blank()
        ,strip.background = element_rect(colour=NA, fill=NA)
        ,axis.line = element_line(color = 'grey')
        ,text=element_text(size=20)
        ,panel.spacing = unit(1, "lines")
        ,axis.text=element_text(size=15)
        ,legend.key.size = unit(0.6,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=20)
        ,legend.position =c(0.89,0.7)
        ,legend.title = element_text(size=10)
  ) +
  facet_wrap(~coef,scales="free")

msep

tiff("AnnexFig1.tiff", height = 4.5, width =6, units = 'in', compression="lzw", res=400)
msep
dev.off()
