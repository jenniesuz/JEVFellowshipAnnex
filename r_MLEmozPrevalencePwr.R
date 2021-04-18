source("r_MLEMozPrevFuncs.R")

# count sims for which >4 AIC difference between model with one prevalence
# and model with two
sims1 <- matrix(nrow=1000,ncol=3)
for(i in 1:length(sims1[,1])){
  sims1[i,] <- randSampleFunc2(tp1=0.001,tp2=0.02,numberPools=80)
}
sum(sims1[,1])


sims2 <- matrix(nrow=1000,ncol=3)
for(i in 1:length(sims2[,1])){
  sims2[i,] <- randSampleFunc2(tp1=0.01,tp2=0.02,numberPools=100)
}
sum(sims2[,1]) # 862
mean(sims2[,2]) # total mosquitoes in sample
