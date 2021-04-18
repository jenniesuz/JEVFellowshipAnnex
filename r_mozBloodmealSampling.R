
library("ggplot2")
library("GLMMmisc")
library("parallel")
library("binom")
library("plyr")
#***************************************************************8
#*BLOODMEALS

intercept<- -7.5
coef <- log(1.062)

odds1 <- exp((intercept + coef*66))/(1+exp((intercept + coef*66)))
odds2 <- exp((intercept + coef*65))/(1+exp((intercept + coef*65)))
odds3 <- exp((intercept + coef*95))/(1+exp((intercept + coef*95)))

odds1/(1+odds1)
odds3/(1+odds3)

oddsRatio <- odds1/odds2
oddsRatio





inds <- 1000
birdRatio <- c(65,75,85,95)
households <- LETTERS[1:10]

moz.data <- expand.grid(inds=inds,birdRatio=birdRatio,households=households)
birdRatio <- sapply(1:length(moz.data$birdRatio),function(x){
  temp <- moz.data[x,]
  return(round(runif(1,temp$birdRatio-5,temp$birdRatio+5)))
})
moz.data$birdRatio <- birdRatio
moz.data$row.id <- factor(paste("row", 1:nrow(moz.data), sep = ""))
moz.data$response <- inds
moz.data$n <- moz.data$response

moz.data<-
  sim.glmm(
    design.data = moz.data,
    fixed.eff =
      list(
        intercept = -7.5,      
        birdRatio =
          log(1.062)),
    rand.V =NULL
    ,distribution = "binomial")     

#****test****#
tbin <- binom.confint(moz.data$response,moz.data$n,methods="exact")
moz.data$mean <- tbin$mean
ggplot(moz.data) +
geom_point(aes(x=birdRatio,y=mean),alpha=0.5,col="red")


#******add in random effect of household
moz.data<-
  sim.glmm(
    design.data = moz.data,
    fixed.eff =
      list(
        intercept = -7.5,      
        birdRatio =
          log(1.062)),
    rand.V =
    inv.mor(
      c(households = 2)) # there is variation between households (MRR=1.5)
    ,distribution = "binomial")    


tbin <- binom.confint(moz.data$response,moz.data$n,methods="exact")
moz.data$mean <- tbin$mean
ggplot(moz.data) +
  geom_point(aes(x=birdRatio,y=mean),alpha=0.5,col="red")

moz.bin <-   lme4::glmer(cbind(response, n - response) ~ birdRatio + (1 | households) 
                         ,  family = "binomial", data = moz.data)

summary(moz.bin)

#********Power analysis assuming no random variation between households********
pwrFunc <- function(...){
  library("GLMMmisc")
  
  inds <- 30
  birdRatio <- c(65,75,85,95)
  households <- LETTERS[1:10]
  
  moz.data <- expand.grid(inds=inds,birdRatio=birdRatio,households=households)
  birdRatio <- sapply(1:length(moz.data$birdRatio),function(x){
    temp <- moz.data[x,]
    return(round(runif(1,temp$birdRatio-5,temp$birdRatio+5)))
  })
  moz.data$birdRatio <- birdRatio
  moz.data$row.id <- factor(paste("row", 1:nrow(moz.data), sep = ""))
  moz.data$response <- inds
  moz.data$n <- moz.data$response
  
  moz.data<-
    sim.glmm(
      design.data = moz.data,
      fixed.eff =
        list(
          intercept = -7.5,      
          birdRatio =
            log(1.062)),
      rand.V =NULL
       ,distribution = "binomial")     
  
  moz.bin <- glm(cbind(response,n-response)~birdRatio,family="binomial",data=moz.data)
  output <- summary(moz.bin)
  cfs<-output$coefficients
  prbirds<-cfs[2,4]
  pvalbirds <- 0
  
  if(length(output$optinfo$conv$lme4$messages)==0){
    if(prbirds<=0.05){pvalbirds<-1}
  }
  return(c(pvalbirds))
}


sim.res50 <- mclapply(1:1000, pwrFunc, mc.cores =1)

l<-do.call(rbind,sim.res50)
l2<-rowSums(l)
sum(l2)






#********Power analysis including some variation between households********
pwrFunc <- function(...){
  library("GLMMmisc")
  inds <- 30
  birdRatio <- c(65,75,85,95)
  households <- LETTERS[1:10]
  
  moz.data <- expand.grid(inds=inds,birdRatio=birdRatio,households=households)
  birdRatio <- sapply(1:length(moz.data$birdRatio),function(x){
    temp <- moz.data[x,]
    return(round(runif(1,temp$birdRatio-5,temp$birdRatio+5)))
  })
  moz.data$birdRatio <- birdRatio
  moz.data$row.id <- factor(paste("row", 1:nrow(moz.data), sep = ""))
  moz.data$response <- inds
  moz.data$n <- moz.data$response
  
  moz.data<-
    sim.glmm(
      design.data = moz.data,
      fixed.eff =
        list(
          intercept = -7.5,      
          birdRatio =
            log(1.062)),
      rand.V =
      inv.mor(
        c(households = 1.5)) 
      ,distribution = "binomial")     
  
  moz.bin <-   lme4::glmer(cbind(response, n - response) ~ birdRatio + (1 | households) ,
  family = "binomial", data = moz.data)
  output <- summary(moz.bin)
  cfs<-output$coefficients
  prbirds<-cfs[2,4]
  pvalbirds <- 0
  
  if(length(output$optinfo$conv$lme4$messages)==0){
    if(prbirds<=0.05){pvalbirds<-1}
  }
  return(c(pvalbirds))
}


sim.res <- mclapply(1:1000, pwrFunc, mc.cores =1)

l<-do.call(rbind,sim.res)
l2<-rowSums(l)
sum(l2)
# 814