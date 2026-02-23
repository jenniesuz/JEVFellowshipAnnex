library(deSolve) 

params <- function(   # Parameters
  maxRecruit = 10^2.5       # mosquito population dynamics model parameters - see separate script of this model for details
  ,baselineRecruit = 0.2    #
  ,interval = 100           #
  ,firstRecruit = 0.4       #
  ,firstMean = 130          #
  ,firstSd = 5              #
  ,secondMean = 150         #
  ,secondSd = 5             #
  ,mu_v= 1/30            # mosquito death rate
  ,years=3               # number of years to run model    
  ,alpha = 1/3           # biting rate
  ,sigma_v = 1/14        # 1/ incubation period
  #************HOST-RELATED PARAMETERS************
  ,p_h1v = 0.05          # prob host to vector transmission
  ,p_h2v = 0.05
  ,p_h3v = 0
  ,mu_h1 = 1/90         # host mortality rate
  ,mu_h2 = 1/365        #
  ,mu_h3 = 1/(365*3)
  ,beta_h1  = 1/90      # host1 birth rate
  ,beta_h2 = 1/365
  ,beta_h3 = 1/(365*3)
  ,pref_h1 = 0.05       # preference of vector for host 1 
  ,pref_h2 = 0.15
  ,pref_h3 = 0.8
  ,p_vh1 = 0.8          # prob vector to host1 transmission
  ,p_vh2 = 0.8
  ,p_vh3 = 0.8
  ,phi_h1 = 0.2         # host recovery rate
  ,phi_h2 = 0.2
  ,phi_h3 = 0.2
)
return(as.list(environment()))

times <- seq(1,365*params()$years,1) # Timesteps for reporting



#********************transmission model***********************************

mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  # total mosquitoes and hosts
  Nv <- Sv + Ev1 + Ev2 + Ev3 + Iv
  Nh1 <- Sh1 + Ih1 + Rh1
  Nh2 <- Sh2 + Ih2 + Rh2
  Nh3 <- Sh3 + Ih3 + Rh3
  # proportion of bloodmeals on each host species
  rho_h1 <- (pref_h1*Nh1)/(pref_h1*Nh1+pref_h2*Nh2+pref_h3*Nh3)
  rho_h2 <- pref_h2*Nh2/(pref_h1*Nh1+pref_h2*Nh2+pref_h3*Nh3)
  rho_h3 <- pref_h3*Nh3/(pref_h1*Nh1+pref_h2*Nh2+pref_h3*Nh3)
  # mosquito population dynamics parameters
  index <- tt
  allTime <- rep(1:365,years)
  time <- allTime[index]
  secondRecruit = 1 - (baselineRecruit + firstRecruit)
  delta1 <- (1/2*pi*firstSd)*exp( -(firstMean - time)^2 / ((2*firstSd)^2) ) 
  delta2 <- (1/2*pi*secondSd)*exp( -(secondMean - time)^2 / ((2*secondSd)^2) ) 
  recruit <- (baselineRecruit*maxRecruit) /interval + firstRecruit*maxRecruit*delta1 + secondRecruit*maxRecruit*delta2
  # ODEs
  # vector
  deriv <- rep(NA,14) 
  deriv[1] <- recruit - alpha*(p_h1v*rho_h1*Ih1/Nh1 + p_h2v*rho_h2*Ih2/Nh2)*Sv - mu_v*Sv  
  deriv[2] <-  alpha*(p_h1v*rho_h1*Ih1/Nh1 + p_h2v*rho_h2*Ih2/Nh2)*Sv - (mu_v + sigma_v*3)*Ev1   # Exposed vectors
  deriv[3] <- (sigma_v*3)*Ev1 - (mu_v + sigma_v*3)*Ev2 
  deriv[4] <- (sigma_v*3)*Ev2 - (mu_v + sigma_v*3)*Ev3
  deriv[5] <- (sigma_v*3)*Ev3 - mu_v*Iv 
  # host 1
  deriv[6] <- beta_h1*Nh1 - alpha*rho_h1*p_vh1*Iv*Sh1/Nh1 - mu_h1*Sh1             # S hosts
  deriv[7] <- alpha*rho_h1*p_vh1*Iv*Sh1/Nh1 - phi_h1*Ih1 - mu_h1*Ih1              # I hosts
  deriv[8] <- phi_h1*Ih1 - mu_h1*Rh1                                              # R hosts
  # host 2 
  deriv[9] <- beta_h2*Nh2 - alpha*rho_h2*p_vh2*Iv*Sh2/Nh2 - mu_h2*Sh2      
  deriv[10] <-  alpha*rho_h2*p_vh2*Iv*Sh2/Nh2 - phi_h2*Ih2 - mu_h2*Ih2  
  deriv[11] <- phi_h2*Ih2 - mu_h2*Rh2 
  # host 2 can become infected and develop antibody response but don't transmit virus to mosquitoes
  deriv[12] <- beta_h3*Nh3 - alpha*rho_h3*p_vh3*Iv*Sh3/Nh3 - mu_h3*Sh3     
  deriv[13] <-  alpha*rho_h3*p_vh3*Iv*Sh3/Nh3 - phi_h3*Ih3 - mu_h3*Ih3  
  deriv[14] <- phi_h3*Ih3 - mu_h3*Rh3 
  
  FOIc<-alpha*rho_h3*p_vh3*Iv/Nh3
  return(list(deriv,Nv=Nv,Nh1=Nh1,Nh2=Nh2,Nh3=Nh3,FOIc=FOIc))
})
#*******************Function to run model**********************************
sim <- function(init=initial, tseq = times, modFunction=mod
                    , parms = params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  return(simDat)
}


simV <- Vectorize(sim)

sim()
