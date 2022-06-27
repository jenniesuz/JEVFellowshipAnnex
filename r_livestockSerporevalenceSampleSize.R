# script to estimate power to compare JEV seroprevalence between predicted high and low risk areas

start.time <- Sys.time()
# load packages
library(GLMMmisc) # available via devtools::install_github("pcdjohnson/GLMMmisc")
library(lme4)
library(parallel)
library(binom)
library(ggplot2)
library(plyr)

#par.tab <- read.csv("parameter.estimates.csv", row.names = 1)

# simulate JEV serology data in each species
dat <- expand.grid(hh = 1:n.hh, village = 1:n.village, n = n)
print(sum(dat$n))
# allocate villages to high and low prevalence in 1:1 ratio 
dat$risk.level <- dat$village %% 2 - 0.5
# simulate seropositives
simdata <-
  sim.glmm(
    design.data = dat, 
    fixed.eff = 
      list(
        intercept = -3.53, #par.tab["(Intercept)",1],
        risk.level = log(OR)),
    distribution = "binomial",
    rand.V = c(hh = 0.4 #par.tab["barcode_hh",1], 
               ,village = 0.6 #par.tab["village",1]
               ))


#************function to visualise data************

vizDatFunc <- function(simdat=simdata){
summariseVillage <- ddply(simdat,.(village,risk.level),summarise,denom=sum(n),num=sum(response))
binomVillage <- binom.confint(summariseVillage$num
                              ,summariseVillage$denom,methods="exact")
binomVillage$villageNum <- 1:length(binomVillage[,1])
binomVillage$risk.level <- as.factor(summariseVillage$risk.level)

ggplot(binomVillage) +
  geom_point(aes(x=villageNum,y=mean*100,col=risk.level)) +
  geom_errorbar(aes(x=villageNum,ymin=lower*100,ymax=upper*100))


summariseRisk <- ddply(simdat,.(risk.level),summarise
                       ,denom=sum(n)
                       ,num=sum(response))
binomRisk <- binom.confint(summariseRisk$num,summariseRisk$denom,methods="exact")

ggplot(binomRisk) +
  geom_point(aes(x=1:2,y=mean*100)) +
  geom_errorbar(aes(x=1:2,ymin=lower*100,ymax=upper*100))
}
#**************************************

#**********************************
# function to simulate data and estimate p-value for null hypothesis 
# that high and low risk areas have the same seroprevalence
res.tab.fn <- function(...) {
  # create template data set
  inter <- -2.197
  dat <- expand.grid(hh = 1:n.hh, village = 1:n.village, n = n)
  # allocate villages to high and low prevalence in 1:1 ratio 
  dat$risk.level <- dat$village %% 2 - 0.5
  # simulate seropositives
  simdat <-
    sim.glmm(
      design.data = dat, 
      fixed.eff = 
        list(
          intercept = inter,
          risk.level = log(OR)),
      distribution = "binomial",
      rand.V = c(hh = 0.4, 
                 village = 0.6))
  
  fit <- glmer(cbind(response, n - response) ~ risk.level + (1 | hh) +(1 | village), family = binomial, data = simdat)
  fit0 <- update(fit, ~ . - risk.level)
  anova(fit, fit0)[2, "Pr(>Chisq)"]
}
#*****************************************


n.village <- 36 #30
n.hh <- 40 #20         # per village
n <- 1             # no of animals per household
OR <- 3
inter <- -2.197

# simulate RVFV serology data in each species
dat <- expand.grid(hh = 1:n.hh, village = 1:n.village, n = n)
print(sum(dat$n))
# allocate villages to high and low prevalence in 1:1 ratio 
dat$risk.level <- dat$village %% 2 - 0.5
# simulate seropositives
simdata <-
  sim.glmm(
    design.data = dat, 
    fixed.eff = 
      list(
        intercept = inter,
        risk.level = log(OR)),
    distribution = "binomial",
    rand.V = c(hh = 0.4, #par.tab["barcode_hh",1], 
               village = 0.6))


vizDatFunc()



#************assume c. 20 - 40% prevalence in high risk areas***********
# and three fold difference in risk
# and OR of 3
# repeat simulations many times and calculate p-value
nsim<-1000
n.village <- 36 #24
n.hh <- 40 #20         # per village
n <- 1 #3             # no of animals per household
# 1800 animals for 30 villages or 1440 for 24 villages
inter <- -2.197

sim.res <- mclapply(1:nsim, res.tab.fn)

# estimate power
apply(do.call("rbind", sim.res) < 0.05, 2, mean)
# 0.849
# 0.948

