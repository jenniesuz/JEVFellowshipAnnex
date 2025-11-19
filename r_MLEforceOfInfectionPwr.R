source("r_MLEforceofInfectionFuncs.R")
library(parallel)
# 
# repSimsFunc <- function(trueFoi=tf,numHosts=nh,ageGroups=ag,maxAge=ma){
#   no_cores <- detectCores() - 1
#   cl <- makeCluster(no_cores)
#   clusterExport(cl,varlist=c("trueFoi","numHosts","ageGroups","maxAge"))
#   
#   simFits <- parLapply(cl,1:500,function(i){
#     pacman::p_load(bbmle, plyr) 
#     source("r_MLEforceofInfectionFuncs.R")
#     ests <- simulateDataFunc(trueFoi=trueFoi,numHosts=numHosts,ageGroups=ageGroups,maxAge=maxAge)
#     return(ests)
#   })
#   
#   stopCluster(cl)
#   
#   simFits1 <- do.call(rbind,simFits)
#   return(simFits1)
# }
# 

repSimsFunc <- function(trueFoi1=trueFoi1,trueFoi2=trueFoi2,trueFoi3=trueFoi3,numHosts=numHosts,ageGroups=ageGroups,maxAge=maxAge){
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  clusterExport(cl,varlist=c("trueFoi1","trueFoi2","trueFoi3","numHosts","ageGroups","maxAge"))
  
  simFits <- parLapply(cl,1:1000,function(i){
    source("r_MLEforceofInfectionFuncs.R")
    ests <- simulateDataFunc(trueFoi1=trueFoi1,trueFoi2=trueFoi2,trueFoi3=trueFoi3,numHosts=numHosts,ageGroups=ageGroups,maxAge=maxAge)
    return(ests)
  })
  
  stopCluster(cl)
  
  simFits1 <- do.call(rbind,simFits)
  return(simFits1)
}

#***************Parameters***************************
trueFoi1 <- 0.1
trueFoi2 <- 0.2
trueFoi3 <- 0.5
numHosts <- 60
ageGroups <- 7
maxAge <- 20
#************************************************************
hosts <- repSimsFunc(trueFoi1=trueFoi1,trueFoi2=trueFoi2,trueFoi3=trueFoi3
                        ,numHosts=numHosts
                        ,ageGroups=ageGroups
                        ,maxAge=maxAge)
h <- data.frame(hosts)
length(h$p1[h$p1>0]) 
length(h$p2[h$p2>0])

