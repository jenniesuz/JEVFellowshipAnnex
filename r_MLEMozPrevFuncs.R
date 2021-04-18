library(binom)
library(bbmle)
library(plyr)
library(ggplot2)
library(parallel)
# simulating MLE of mosquito infection prevalence

#**********************Fit functions*******************
#*#*************One sample************
nll.binom <- function(logPropInf,dat){ 
  num <- dat$num                #number pools positive
  denom <- dat$denom            #number of pools screened
  
  denomZero <- denom[num %in% 0] # select from denom those that were all negative
  denomOne <- denom[num %in% 1]
  
  p <- exp(logPropInf)
  
  ll <- sum( log((1-p)^(denomZero)) ) + sum ( log(1 - (1-p)^(denomOne)) )
  
  return(-ll)
}
#************two samples****************
nll.binom2 <- function(logPropInf1,logPropInf2,dat){ 
  dat1 <- dat[dat$prev %in% "low",]
  dat2 <- dat[dat$prev %in% "high",]
  
  num1 <- dat1$num                #number pools positive
  denom1 <- dat1$denom            #number of pools screened
  num2 <- dat2$num
  denom2 <- dat2$denom
  
  denomZero1 <- denom1[num1 %in% 0] # select from denom those that were all negative
  denomOne1 <- denom1[num1 %in% 1]
  
  denomZero2 <- denom2[num2 %in% 0] # select from denom those that were all negative
  denomOne2 <- denom2[num2 %in% 1]
  
  
  p1 <- exp(logPropInf1)
  p2 <- exp(logPropInf2)
  
  ll1 <- sum( log((1-p1)^(denomZero1)) ) + sum ( log(1 - (1-p1)^(denomOne1)) )
  ll2 <- sum( log((1-p2)^(denomZero2)) ) + sum ( log(1 - (1-p2)^(denomOne2)) )
  ll <- ll1 + ll2
  
  return(-ll)
}
#*****************************************

#******************Function to generate random samples (single prev)*************
randSampleFunc <- function(tp=truePrev,numberPools){
  
  denom <- sample(20:50,numberPools,replace=T)
  
  num <- sapply(denom,function(x){
    randomBinom <- rbinom(x,1,prob=tp)
    randomBinom <- max(randomBinom)
    if (randomBinom >=1 ){return(1)}else{
      return(0)
    }
  })
  
  dat <- cbind.data.frame(denom,num)
  
  totMoz <- sum(denom)
  
  mleFit <- mle2(function(par1){nll.binom(par1,dat)}             
                 ,start=list(par1=log(0.05)))
  
  est <- as.numeric(exp(coef(mleFit)))
  ci <- exp(as.numeric(confint(mleFit)))
  ciWidth <- ci[2] - ci[1]
  bias <- est - tp
  return(c(est=est,ci=ci,bias=bias,ciWidth=ciWidth,totMoz=totMoz,truePrev=tp,numberPools=numberPools))
}

#******************Function to generate random samples (two prev)*************
randSampleFunc2 <- function(tp1=truePrev1,tp2=truePrev2,numberPools){
  
  denom1 <- sample(20:50,numberPools,replace=T) # low prevalence
  denom2 <- sample(20:50,numberPools,replace=T) # high prevalence
  
  num1 <- sapply(denom1,function(x){            # low
    randomBinom <- rbinom(x,1,prob=tp1)
    randomBinom <- max(randomBinom)
    if (randomBinom >=1 ){return(1)}else{
      return(0)
    }
  })
  num2 <- sapply(denom2,function(x){            # high
    randomBinom <- rbinom(x,1,prob=tp2)
    randomBinom <- max(randomBinom)
    if (randomBinom >=1 ){return(1)}else{
      return(0)
    }
  })
  

  dat1 <- cbind.data.frame(denom1,num1)
  names(dat1) <- c("denom","num")
  dat2 <- cbind.data.frame(denom2,num2)
  names(dat2) <- c("denom","num")
  dat1$prev <- "low"
  dat2$prev <- "high"
  
  dat <- rbind.data.frame(dat1,dat2)
  
  totMoz <- sum(denom1)
  
  mleFit1 <- mle2(function(par1){nll.binom(par1,dat)}  # assume single prev             
                 ,start=list(par1=log(0.05)))
  mleFit2 <- mle2(function(par1,par2){nll.binom2(par1,par2,dat)}  # assume two prev             
                  ,start=list(par1=log(0.05),par2=log(0.05)))
  
  p <- 0
  anovaRes <- anova(mleFit1,mleFit2)
  anovaRes <- anovaRes[2,5]
  if(anovaRes<=0.05){
    p<-1
  }
  return(c(p,totMoz,numberPools))
}

#********Function to calculate power********
powerFunc <- function(dat){
  numMoz <- mean(dat$totMoz)
  variance <- var(dat$est)
  meanBias <- mean(dat$bias)
  mse <- meanBias^2 + variance
  ciw <- mean(dat$ciWidth)
  cov <- sapply(1:length(dat[,1]),function(x){
    temp <- dat[x,]
    if(temp$truePrev > temp$ci1 & temp$truePrev < temp$ci2){
      return(1)}else{
        return(0)
      }
  })
  coverage <- sum(cov)/length(cov)
  
  return(c(variance=variance,meanBias=meanBias,mse=mse,coverage=coverage,numMoz=numMoz,ciw=ciw))
}