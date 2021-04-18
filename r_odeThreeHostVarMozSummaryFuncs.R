
plotSimFunc <- function(simDat){

  mozPerHost <- simDat$Nv/ ( simDat$Nh1 + simDat$Nh2 + simDat$Nh3 )

  hostProps <- cbind.data.frame(
    time=rep(simDat$time,9)
  
    ,host=c(rep("Host1",length(simDat$Ih1)*3)
          ,rep("Host2",length(simDat$Ih1)*3)
          ,rep("Host3",length(simDat$Ih1)*3))
  
    ,class=c(rep(c(rep("Susceptible",length(simDat$Ih1))
                 ,rep("Infected",length(simDat$Ih1))
                 ,rep("Resistant",length(simDat$Ih1))),3)
    )
  
    ,value=c((simDat$Sh1/simDat$Nh1)*100
           ,(simDat$Ih1/simDat$Nh1)*100
           ,(simDat$Rh1/simDat$Nh1)*100
           ,(simDat$Sh2/simDat$Nh2)*100
           ,(simDat$Ih2/simDat$Nh2)*100
           ,(simDat$Rh2/simDat$Nh2)*100
           ,(simDat$Sh3/simDat$Nh3)*100
           ,(simDat$Ih3/simDat$Nh3)*100
           ,(simDat$Rh3/simDat$Nh3)*100
    )
  )

  mozProps <- cbind.data.frame(time=simDat$time
                             ,pIm=(simDat$Ev1+simDat$Ev2+simDat$Ev3+simDat$Iv)/(simDat$Nv)
                             ,n1000=(simDat$Ev1+simDat$Ev2+simDat$Ev3+simDat$Iv)/(simDat$Nv)*1000
                             ,mph= mozPerHost
  )


  mozPerHost <- ggplot(mozProps) +
    geom_line(aes(x=time,y=mph)) +
    xlab("Time") +
    ylab("Mosquitoes per host") + 
    # xlim(730,1095) +
    theme_set(theme_bw()) +
    theme(
      axis.line = element_line(color = 'black')
      ,text=element_text(size=10)
      ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
      ,axis.text=element_text(size=8)
      ,legend.key.size = unit(0.8,"line")
      ,legend.background = element_blank()
      ,legend.text=element_text(size=9)
      ,legend.position =c(0.8,0.2)
      ,legend.title = element_blank()
      ,strip.background = element_rect(colour="white", fill="white")
      ,panel.border = element_blank()
  )

  hostPlot <- ggplot(hostProps) +
    geom_line(aes(x=time,y=value,col=class,linetype=host)) +
    xlab("Time") +
    ylab("Percentage of hosts") + 
    # xlim(730,1095) +
    theme_set(theme_bw()) +
    theme(
      axis.line = element_line(color = 'black')
      ,text=element_text(size=10)
      ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
      ,axis.text=element_text(size=8)
      ,legend.key.size = unit(0.8,"line")
      ,legend.background = element_blank()
      ,legend.text=element_text(size=9)
      ,legend.position =c(0.8,0.4)
      ,legend.title = element_blank()
      ,strip.background = element_rect(colour="white", fill="white")
      ,panel.border = element_blank()
    )

  mozPlot <- ggplot(mozProps) +
    geom_line(aes(x=time,y=n1000)) +
    xlab("Time") +
    # xlim(730,1095) +
    ylab("Proportion of infected mosquitoes * 1000") + 
    theme_set(theme_bw()) +
    theme(
      axis.line = element_line(color = 'black')
      ,text=element_text(size=10)
      ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
      ,axis.text=element_text(size=8)
      ,legend.key.size = unit(0.8,"line")
      ,legend.background = element_blank()
      ,legend.text=element_text(size=9)
      ,legend.position =c(0.8,0.2)
      ,legend.title = element_blank()
      ,strip.background = element_rect(colour="white", fill="white")
      ,panel.border = element_blank()
    )


  foiCPlot <- ggplot(simDat) +
   geom_line(aes(x=time,y=FOIc))  +
   xlab("Time") +
    ylim(0,0.5) +
    ylab("Force of infection from mosquitoes to cattle") + 
    theme_set(theme_bw()) +
    theme(
      axis.line = element_line(color = 'black')
      ,text=element_text(size=10)
      ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
      ,axis.text=element_text(size=8)
      ,legend.key.size = unit(0.8,"line")
      ,legend.background = element_blank()
      ,legend.text=element_text(size=9)
      ,legend.position =c(0.8,0.2)
      ,legend.title = element_blank()
      ,strip.background = element_rect(colour="white", fill="white")
      ,panel.border = element_blank()
  )


  grid.arrange(mozPerHost,hostPlot,mozPlot,foiCPlot,ncol=2)
  return(mozProps)
}


