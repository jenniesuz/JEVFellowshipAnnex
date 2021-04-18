# Simulations to assess numbers of traps required to detect differences 
# in mosquito abundance between bird types
library("ggplot2")
library("GLMMmisc")
library("parallel")
library("plyr")
library("MASS")

lighttraps <- read.csv("RajshahiCxt.csv")
byDay <- ddply(lighttraps,.(Date),varCount=var(Females))
range(byDay[,2]) # 1 to 196
mean(byDay[,2]) # c. 100

lme4::glmer(Females ~ (1 | Date),
            family = "poisson", data =lighttraps)
#1.58

#***************************************************************
numTraps<-5
traps <- 1:numTraps
bird <- c("chicken","duck","pigeon")
habitat <- c("household","rice")

moz.data <-
  expand.grid(traps=traps,habitat=habitat,bird=bird)

moz.data$row.id <- factor(paste("row", 1:nrow(moz.data), sep = ""))
moz.data<-
  sim.glmm(
    design.data = moz.data,
    fixed.eff =
      list(
        intercept = log(20),      
        habitat =
          log(                     
            c(household = 1,               
              rice = 0.5)),
        bird=log(
          c(pigeon=1,
            chicken=0.5,
            duck=1.5)
        )),
    rand.V =NULL
    ,distribution = "negbinomial"
    ,theta=0.5
  )      # we are simulating 



moz.pois <-
  glm.nb(response ~ bird + habitat, data = moz.data)
output <- summary(moz.pois)
cfs<-output$coefficients

trueVals<-c(20,1.5,0.5,0.5)
exp(cfs[,1])

trueVals - exp(cfs[,1])

#*****************Simulations to determine number of traps - poisson************
pwrFunc <- function(...){
  library("GLMMmisc")
  
  traps <- 1:40
  bird <- c("chicken","duck","pigeon")
  habitat <- c("household","rice")
  
  moz.data <-
    expand.grid(traps=traps,habitat=habitat,bird=bird)
  
  moz.data$row.id <- factor(paste("row", 1:nrow(moz.data), sep = ""))
  
  moz.data<-
    sim.glmm(
      design.data = moz.data,
      fixed.eff =
        list(
          intercept = log(20),      
          habitat =
            log(                     
              c(household = 1,               
                rice = 0.5)),
          bird=log(
            c(pigeon=1,
              chicken=0.5,
              duck=1.5)
          )),
      rand.V =NULL
      ,distribution = "poisson"
 )      # we are simulating 
  
  
  moz.pois <-
    glm(response ~ bird + habitat, data = moz.data, family="poisson")
  output <- summary(moz.pois)
  cfs<-output$coefficients
  prduck<-cfs[2,4]
  prpigeon<-cfs[3,4]
  prhabitat<-cfs[4,4]
  pvalduck <- 0
  pvalpigeon <- 0
  pvalhabitat <- 0
  
  if(length(output$optinfo$conv$lme4$messages)==0){
    if(prduck<=0.05){pvalduck<-1}
    if(prpigeon<=0.05){pvalpigeon<-1}
    if(prhabitat<=0.05){pvalhabitat<-1}
  }
  return(c(pvalduck,pvalpigeon,pvalhabitat))
}

sim.res <- mclapply(1:1000, pwrFunc, mc.cores =1)
l<-do.call(rbind,sim.res)
l2<-rowSums(l)
length(l2[l2==3])
# 1000

#*****************Simulations to determine number of traps - negative binomial************
pwrFunc <- function(...){
  library("GLMMmisc")
  
  traps <- 1:40
  bird <- c("chicken","duck","pigeon")
  habitat <- c("household","rice")
  
  moz.data <-
    expand.grid(traps=traps,habitat=habitat,bird=bird)
  
  moz.data$row.id <- factor(paste("row", 1:nrow(moz.data), sep = ""))
  
  moz.data<-
    sim.glmm(
      design.data = moz.data,
      fixed.eff =
        list(
          intercept = log(20),      
          habitat =
            log(                     
              c(household = 1,               
                rice = 0.5)),
          bird=log(
            c(pigeon=1,
              chicken=0.5,
              duck=1.5)
          )),
      rand.V =NULL
      ,distribution = "negbinomial"
      ,theta=0.5
    )      # we are simulating 
  
  
  moz.pois <-
    glm.nb(response ~ bird + habitat, data = moz.data)
  output <- summary(moz.pois)
  cfs<-output$coefficients
  prduck<-cfs[2,4]
  prpigeon<-cfs[3,4]
  prhabitat<-cfs[4,4]
  pvalduck <- 0
  pvalpigeon <- 0
  pvalhabitat <- 0
  
  if(length(output$optinfo$conv$lme4$messages)==0){
    if(prduck<=0.05){pvalduck<-1}
    if(prpigeon<=0.05){pvalpigeon<-1}
    if(prhabitat<=0.05){pvalhabitat<-1}
  }
  return(c(pvalduck,pvalpigeon,pvalhabitat))
}

sim.res <- mclapply(1:1000, pwrFunc, mc.cores =1)
l<-do.call(rbind,sim.res)
l2<-rowSums(l)
length(l2[l2==3])

