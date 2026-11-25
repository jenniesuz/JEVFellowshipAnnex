
# Simulations to assess the power to detect differences in mosquito abundance between traps
# containing either duck, chicken or pigeon, informing study design

library("ggplot2")    # packages required
library("GLMMmisc")
library("parallel")
library("plyr")
library("MASS")
library("GLMMmisc")
#library(devtools)
#devtools::install_github("pcdjohnson/GLMMmisc")
#library(GLMMmisc)

#*******************Single simulation**********************************
# Study design
numTraps <- 40                                                    # number of traps per bird and area
traps <- 1:numTraps   
bird <- c("chicken","duck","pigeon")                              # bird type
habitat <- c("household","rice")                                  # area

moz.data <- expand.grid(traps=traps,habitat=habitat,bird=bird)    # create table
moz.data$row.id <- factor(paste("row", 1:nrow(moz.data), sep = ""))

# Simulate data
moz.data <-
  sim.glmm(
    design.data = moz.data,
    fixed.eff =                  # assumed fixed effects
      list(
        intercept = log(20),      
        habitat =
          log(c(household = 1,                
              rice = 0.5)),
        bird=log(
          c(pigeon = 1,
            chicken = 0.5,
            duck = 1.5)
        )),
    rand.V = NULL
    ,distribution = "negbinomial"
    ,theta = 0.5                     # shape
  )     


# Plot
ggplot(moz.data) +
  geom_boxplot((aes(x=bird,y=response))) +
  ylab("Number of mosquitoes per trap per day") +
  facet_wrap(~habitat)

# Fit model
moz.mod <-
  glm.nb(response ~ bird + habitat, data = moz.data)
output <- summary(moz.mod)
output
cfs <- output$coefficients


#*****************Simulations to determine number of traps - negative binomial************
pwrFunc <- function(...){
  library("GLMMmisc")
  numTraps <- 40                                                    # number of traps per bird and area
  traps <- 1:numTraps   
  bird <- c("chicken","duck","pigeon")                              # bird type
  habitat <- c("household","rice")                                  # area
  
  moz.data <- expand.grid(traps=traps,habitat=habitat,bird=bird)    # create table
  moz.data$row.id <- factor(paste("row", 1:nrow(moz.data), sep = ""))
  
  moz.data <-
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
    )      
  
  
  moz.mod <- glm.nb(response ~ bird + habitat, data = moz.data)
  output <- summary(moz.mod)
  cfs <- output$coefficients
  prduck <- cfs[2,4]
  prpigeon <- cfs[3,4]
  prhabitat <- cfs[4,4]
  pvalduck <- 0
  pvalpigeon <- 0

  if(length(output$optinfo$conv$lme4$messages)==0){
    if(prduck<=0.05){pvalduck<-1}
    if(prpigeon<=0.05){pvalpigeon<-1}
  }
  return(c(pvalduck,pvalpigeon))
}

# Assess power to detect a difference between bird types
sim.res <- mclapply(1:1000, pwrFunc, mc.cores=1)
l <- do.call(rbind,sim.res)
l2 <- rowSums(l)
length(l2[l2==2])/1000*100

