### Ellie Broadman
### December 2020

##Part 1: Calculating the amount of groundwater inflow (m^3/year) at steady state, without uncertainty, using the model of Krabbenhoft (1990)

#the model requires the following 9 inputs:
  #L is the volume of lake water in m^3;
  #dL is the d18O of lake water in per mille VSMOW;
  #SA is the lake's surface area in m^2;
  #P is the amount of on-lake precipitation in m;
  #dP is the d18O of average annual precip in per mille VSMOW;
  #TA is the air temperature in °C;
  #TL is the lake water temperature in °C;
  #RH is relative humidity in %;
  #and dG is the d18O of groundwater entering the lake in per mille VSMOW.

library(gridExtra)
#if we use Matanuska instead (USE THIS):
KrabbenhoftB = function(dL,SA,E,P,dP,TA,TL,RH,dG){
  #calculate on-lake precip and evap (in m^3)
  precip = P*SA
  evap = E*SA
  #E = 0.7*365*SA #259970.04
  #empirically derive d18O of evaporated water vapor (per mille VSMOW)
  TA = TA+273.15
  TL = TL+273.15
  esa = 6.108*(2.71828^((17.27*TA)/(TA+237.7)))
  lnalpha = (0.35041*(10^6/TA^3)) - (1.6664*(10^3/TA^2)) + (6.7123*(1/TA)) - (7.685*10^-3)
  alpha = exp(1)^lnalpha 
  h = RH*(TA/TL)
  epsilonP = 1000*(1-alpha)
  dA = dP - epsilonP
  epsilonK = 14.2*(1-h)
  epsilon = epsilonK + epsilonP
  dE = ((alpha*dL) - (h*dA) - epsilon) / (1-h+(0.001*epsilonK))
  #and finally calculate the amount of groundwater inflow:
  X = (precip*(dL-dP)+evap*(dE-dL))/(dG-dL)
  return(X)
}

#for annual:
KrabbenhoftB(dL = -13.9, SA = 590841, E = 0.44, P=0.59 ,dP = -17.7 ,TA = 3,TL = 9.5, RH = 89, dG = -16.5)
KrabbenhoftB(dL = -13.9, SA = 590841, E = 0.35, P=0.64 ,dP = -17.7 ,TA = 3,TL = 9.5, RH = 89, dG = -16.5)
590841*0.44
590841*0.59

#for summer: (also need to change 365 to 184 in function)
#Krabbenhoft(dL = -13.9, SA = 590841, P=0.32 ,dP = -15.4 ,TA = 10,TAmin = -6, TAmax = 21,TL = 15, RH = 89, dG = -16.5, z = 94,u = 4, A = 65.5)


##Part 2: Sensitivity testing for the Early and Middle/Late Holocene (steady state)

#first, to derive d18Olake from d18Odiatom:
#we will test using the modern relationships which best represents this system
#Crespin et al 2010
#t (°C) =245.3 - 6.25 (d18Oopal - d18Owater)
d18Oopal = 25.31
dL = -13.9

Twater = 245.3 - 6.25*(d18Oopal-dL)
#if we use 4°C, and average the most recent 3 diatom isotopes samples (avg 24.6), we get:
Twater = 4
d18Oopal = 24.6 
dL_modern = (((Twater - 245.3)/-6.25)-d18Oopal)*-1

#for the early Holocene:
#using an average of the last 3 samples (9.7-9.2 ka)
d18Oopal = 22.9
dL_earliest = (((Twater - 245.3)/-6.25)-d18Oopal)*-1
#using the average of 9.7-7.3 ka
#d18Oopal = 24.5
#dL_early = (((Twater - 245.3)/-6.25)-d18Oopal)*-1

#for the middle to late Holocene:
#using the average of many samples (7.3-modernish)
d18Oopal = 27.6
dL_midlate = (((Twater - 245.3)/-6.25)-d18Oopal)*-1

#derive +/- values from the SD of the diatom d18O values
#incorporate analytical uncertainty??

#second, calculate lake surface area and volume based on a lower lake level
#minimum lake level in the early Holocene would have barely covered the site for Core 4, but not core 3 (~0.25 water depth at core 4 site, or 3.75 m lower than present)
#maximum lake level would be 1 m lower than present (maximum Chara depth at that site)

SA_early_min = 371697
#V_early_min = 1695867.563
SA_early_max = 532847
#V_early_max = 2980612.906
#the mean of the min and max:
SA_early = (371697+532847)/2

#assume modern lake level for middle/late Holocene ± 1 m?
SA_late_max = 607460
SA_late_min = 590841
SA_late = (607460+590841)/2

## FOR EARLY HOLOCENE
## grid cell for Kelly Lake is 30 (lon) and 62.5 (lat)
#gotta turn all that into a function gurl!!!

#USE THIS???
Krabbenhoft5 = function(dL,SA,E,dP,P,TA,TL,RH){
  #convert P and E to volume:
  P = c(P*1.25, P, P*0.75)
  E = c(E*1.25, E, E*0.75)
  precip = P*SA
  evap = E*SA
  #empirically derive d18O of evaporated water vapor (per mille VSMOW), for 3 different d18O scenarios
  TA = TA+273.15
  TL = TL+273.15
  esa = 6.108*(2.71828^((17.27*TA)/(TA+237.7)))
  lnalpha = (0.35041*(10^6/TA^3)) - (1.6664*(10^3/TA^2)) + (6.7123*(1/TA)) - (7.685*10^-3)
  alpha = exp(1)^lnalpha 
  h = RH*(TA/TL)
  epsilonP = 1000*(1-alpha)
  epsilonK = 14.2*(1-h)
  epsilon = epsilonK + epsilonP
  dA = dP - epsilonP
  dE = ((alpha*dL) - (h*dA) - epsilon) / (1-h+(0.001*epsilonK))
  #groundwater as related to precip
  dG = c(dP-1.5, dP, dP+1.5)
  #calculate G based off ranges of dP and dG
  G1 = matrix(NA, length(dG),length(evap))
  G2 = matrix(NA, length(dG),length(evap))
  G3 = matrix(NA, length(dG),length(evap))
  for (g in 1:3) {
    for (e in 1:3) {
    G1[g,e] = ((precip[1]*(dL-dP)+evap[e]*(dE-dL))/(dG[g]-dL))
    }
  }
  for (g in 1:3) {
    for (e in 1:3) {
      G2[g,e] = ((precip[2]*(dL-dP)+evap[e]*(dE-dL))/(dG[g]-dL))
    }
  }
  for (g in 1:3) {
    for (e in 1:3) {
      G3[g,e] = ((precip[3]*(dL-dP)+evap[e]*(dE-dL))/(dG[g]-dL))
    }
  }
  G = cbind(G1,G2,G3)
  row.names(G) = dG
  colnames(G) = c(rep(precip[1],times = 3), rep(precip[2],times = 3),rep(precip[3],times = 3))
  G = as.data.frame(cbind(dG,G))
  library(ggplot2)
  library(scales)
  plot = ggplot()+
    geom_ribbon(aes(x=G[,1], ymax=G[,2], ymin = G[,4],colour = "PMax", fill = "Pmax"), size =1,alpha = 0.1, linetype = 2)+
    geom_line(aes(x=G[,1], y=G[,3], colour = "PMax"),size =1)+
    #geom_line(aes(x=G[,1], y=G[,4], colour = "PMax"), linetype = 2)+
    geom_ribbon(aes(x=G[,1], ymax=G[,5], ymin = G[,7],colour = "P", fill = "P"), size =1,alpha = 0.1, linetype = 2)+
    geom_line(aes(x=G[,1], y=G[,6], colour = "P"),size =1)+
    #geom_line(aes(x=G[,1], y=G[,7], colour = "P"), linetype = 2)+
    geom_ribbon(aes(x=G[,1], ymax=G[,8], ymin = G[,10],colour = "Pmin", fill = "Pmin"), size =1,alpha = 0.1, linetype = 2)+
    geom_line(aes(x=G[,1], y=G[,9], colour = "Pmin"),size =1)+
    #geom_line(aes(x=G[,1], y=G[,10], colour = "Pmin"), linetype = 2)+
    geom_hline(aes(yintercept = 0), colour = "black")+
    geom_hline(aes(yintercept = 118000), linetype = "dotted", size = 1, colour = "black")+
    xlab("Groundwater d18O")+
    ylab("Groundwater inflow")+
    #scale_y_continuous(labels = comma, breaks = c(0,100000,200000,300000,400000,500000))+
    scale_color_discrete(name = "Precipitation rate", labels = c("P","PMax","PMin"))+
    coord_cartesian(ylim = c(0, 500000)) +
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.position = "none", legend.text = element_text(size = 18), legend.title = element_text(size = 18), panel.background = element_rect(fill = 'white', colour = 'black'))
  return(plot)
}

#450000*0.5
#600000*0.61
#590841*0.59
early_summer = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -14)
early_1 = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -15)
early_2 = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -16)
early_same = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -17)
early_3 = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -18)
early_4 = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -19)
early_winter = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -20)
early_5 = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -21)
early_6 = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -22)
#if dP and dG were much lower, there would be less groundwater throughflow
#if 

grid.arrange(modern_summer,early_summer,late_summer,
             modern_1, early_1, late_1,
             modern_2, early_2, late_2,
             modern_same, early_same, late_same,
             modern_3, early_3, late_3,
             modern_4, early_4, late_4,
             modern_winter, early_winter, late_winter,
             modern_5, early_5, late_5,
             modern_6, early_6, late_6, ncol = 3)



grid.arrange(modern_summer, modern_same, late_same,early_summer,early_same,early_winter,late_summer,late_same,late_winter, nrow = 3)


late_summer = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -14)
late_1 = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -15)
late_2 = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -16)
late_same = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -17)
late_3 = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -18)
late_4 = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -19)
late_winter = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -20)
late_5 = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -21)
late_6 = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -22)

#modern:
modern_summer = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -14)
modern_1 = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -15)
modern_2 = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -16)
modern_same = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -17)
modern_3 = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -18)
modern_4 = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -19)
modern_winter = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -20)
modern_5 = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -21)
modern_6 = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -22)


library(gridExtra)
grid.arrange(modern_summer, modern_same, late_same,early_summer,early_same,early_winter,late_summer,late_same,late_winter, nrow = 3)
?grid.arrange

P9k= 1.74798*10^-5
P4k= 2.09953*10^-5
P0k= 2.02376*10^-5
(2.02376*10^-5)*31536000/1000
0.64
P4k*31536000/1000
P9k*31536000/1000


E9k= 0.652746
E4k= 0.990383
E0k= 0.967294

E9k/1000*365
E4k/1000*365
E0k/1000*365

T0k = -1.9
T4k = -1.9
T9k = -2.4
#4k is 2 degrees cooler than recent years; 9k is 0.5 cooler than 4k Holocene

(P0k - P9k)/P0k
#14% less in the early holocene
(P4k - P0k)/P0k
#4% more in the late holocene

(E0k - E9k)/E0k
#33% less in the early holocene
(E4k - E0k)/E0k
#2% more in the late holocene

dPmin = c(-14,-17,-20)
dGmin = c(dPmin)
for (i in 1:length(dP)*3) {
  dG[i,i+1,i+2] = c(dP[i]-2, dP[i], dP[i]+2)
}



#scratch

precip = 0.2
dL = -16
dP = -17
evap = 0.2
dE = -20
dG = -17
X = (precip*(dL-dP)+evap*(dE-dL))/(dG-dL)
X

precip = 0.4
dL = -14
dP = -12
evap = 0.5
dE = -20
dG = -11
X = (precip*(dL-dP)+evap*(dE-dL))/(dG-dL)
X

#dE is always lower than dL
#model assumes dL is higher than dG, which is a reasonable assumption because X

Krabbenhoft5 = function(dL,SA,E,P,TA,TL,RH){
  #convert P and E to volume:
  P = c(P*1.25, P, P*0.75)
  E = c(E*1.25, E, E*0.75)
  precip = P*SA
  evap = E*SA
  #pick values for dP based on dL
  dP = c(dL-0.5, dL-2.5, dL-4.5)
  #empirically derive d18O of evaporated water vapor (per mille VSMOW), for 3 different d18O scenarios
for (p in 1:3) {
  TA = TA+273.15
  TL = TL+273.15
  esa = 6.108*(2.71828^((17.27*TA)/(TA+237.7)))
  lnalpha = (0.35041*(10^6/TA^3)) - (1.6664*(10^3/TA^2)) + (6.7123*(1/TA)) - (7.685*10^-3)
  alpha = exp(1)^lnalpha 
  h = RH*(TA/TL)
  epsilonP = 1000*(1-alpha)
  epsilonK = 14.2*(1-h)
  epsilon = epsilonK + epsilonP
  dA = dP - epsilonP
  dE = ((alpha*dL) - (h*dA) - epsilon) / (1-h+(0.001*epsilonK))
  #groundwater as related to precip
  dG = c(dP-1.5, dP, dP+1.5)
  #calculate G based off ranges of dP and dG
  G1 = matrix(NA, length(dG),length(evap))
  G2 = matrix(NA, length(dG),length(evap))
  G3 = matrix(NA, length(dG),length(evap))
  for (g in 1:3) {
    for (e in 1:3) {
      G1[g,e] = ((precip[1]*(dL-dP)+evap[e]*(dE-dL))/(dG[g]-dL))
    }
  }
  for (g in 1:3) {
    for (e in 1:3) {
      G2[g,e] = ((precip[2]*(dL-dP)+evap[e]*(dE-dL))/(dG[g]-dL))
    }
  }
  for (g in 1:3) {
    for (e in 1:3) {
      G3[g,e] = ((precip[3]*(dL-dP)+evap[e]*(dE-dL))/(dG[g]-dL))
    }
  }
  G = cbind(G1,G2,G3)
  row.names(G) = dG
  colnames(G) = c(rep(precip[1],times = 3), rep(precip[2],times = 3),rep(precip[3],times = 3))
}
  G = as.data.frame(cbind(dG,G))
  library(ggplot2)
  library(scales)
  plot = ggplot()+
    geom_ribbon(aes(x=G[,1], ymax=G[,2], ymin = G[,4],colour = "PMax", fill = "Pmax"), size =1,alpha = 0.1, linetype = 2)+
    geom_line(aes(x=G[,1], y=G[,3], colour = "PMax"),size =1)+
    #geom_line(aes(x=G[,1], y=G[,4], colour = "PMax"), linetype = 2)+
    geom_ribbon(aes(x=G[,1], ymax=G[,5], ymin = G[,7],colour = "P", fill = "P"), size =1,alpha = 0.1, linetype = 2)+
    geom_line(aes(x=G[,1], y=G[,6], colour = "P"),size =1)+
    #geom_line(aes(x=G[,1], y=G[,7], colour = "P"), linetype = 2)+
    geom_ribbon(aes(x=G[,1], ymax=G[,8], ymin = G[,10],colour = "Pmin", fill = "Pmin"), size =1,alpha = 0.1, linetype = 2)+
    geom_line(aes(x=G[,1], y=G[,9], colour = "Pmin"),size =1)+
    #geom_line(aes(x=G[,1], y=G[,10], colour = "Pmin"), linetype = 2)+
    geom_hline(aes(yintercept = 0), colour = "black")+
    geom_hline(aes(yintercept = 118000), linetype = "dotted", size = 1, colour = "black")+
    xlab("Groundwater d18O")+
    ylab("Groundwater inflow")+
    #scale_y_continuous(labels = comma, breaks = c(0,100000,200000,300000,400000,500000))+
    scale_color_discrete(name = "Precipitation rate", labels = c("P","PMax","PMin"))+
    coord_cartesian(ylim = c(0, 500000)) +
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.position = "none", legend.text = element_text(size = 18), legend.title = element_text(size = 18), panel.background = element_rect(fill = 'white', colour = 'black'))
  return(plot)
}

#450000*0.5
#600000*0.61
#590841*0.59
early_summer = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -14.7)
early_1 = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -15.7)
early_2 = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -16.7)
early_same = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -17.7)
early_3 = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -18.7)
early_4 = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -19.7)
early_winter = Krabbenhoft5(dL = -15.7,SA = 450000,E = 0.32, P = 0.5,TA = 0.5,TL = 7, RH = 91.5, dP = -20.7)
#if dP and dG were much lower, there would be less groundwater throughflow
#if 

grid.arrange(early_summer,early_1,early_2, early_same, early_3, early_4, early_winter, nrow = 1)


grid.arrange(modern_summer, modern_same, late_same,early_summer,early_same,early_winter,late_summer,late_same,late_winter, nrow = 3)


late_summer = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -14.7)
late_same = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -17.7)
late_winter = Krabbenhoft5(dL = -11,SA = 600000,E = 0.45,P = 0.61,TA = 1,TL = 7.5, RH = 88.7, dP = -20.7)

#modern:
modern_summer = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -16.4)
modern_same = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -17.7)
modern_winter = Krabbenhoft5(dL = -13.9,SA = 590841,E = 0.44, P = 0.59,TA = 3,TL = 9.5, RH = 89, dP = -20.7)
-15.7+1.5
-15.7-1.5
-17.7+1.5
library(gridExtra)
grid.arrange(modern_summer, modern_same, late_same,early_summer,early_same,early_winter,late_summer,late_same,late_winter, nrow = 3)
?grid.arrange

P9k= 1.74798*10^-5
P4k= 2.09953*10^-5
P0k= 2.02376*10^-5
(2.02376*10^-5)*31536000/1000
0.64
P4k*31536000/1000
P9k*31536000/1000


E9k= 0.652746
E4k= 0.990383
E0k= 0.967294

E9k/1000*365
E4k/1000*365
E0k/1000*365

T0k = -1.9
T4k = -1.9
T9k = -2.4
#4k is 2 degrees cooler than recent years; 9k is 0.5 cooler than 4k Holocene

(P0k - P9k)/P0k
#14% less in the early holocene
(P4k - P0k)/P0k
#4% more in the late holocene

(E0k - E9k)/E0k
#33% less in the early holocene
(E4k - E0k)/E0k
#2% more in the late holocene

dPmin = c(-14,-17,-20)
dGmin = c(dPmin)
for (i in 1:length(dP)*3) {
  dG[i,i+1,i+2] = c(dP[i]-2, dP[i], dP[i]+2)
}



#scratch

precip = 0.2
dL = -16
dP = -17
evap = 0.2
dE = -20
dG = -17
X = (precip*(dL-dP)+evap*(dE-dL))/(dG-dL)
X

precip = 0.4
dL = -14
dP = -12
evap = 0.5
dE = -20
dG = -11
X = (precip*(dL-dP)+evap*(dE-dL))/(dG-dL)
X

#dE is always lower than dL
#model assumes dL is higher than dG, which is a reasonable assumption because X


Krabbenhoft5 = function(dL,dP,SA,E,P,TA,TL,RH){
  #convert P and E to volume:
  P = c(P*1.25, P, P*0.75)
  E = c(E*1.25, E, E*0.75)
  precip = P*SA
  evap = E*SA
  #empirically derive d18O of evaporated water vapor (per mille VSMOW), for 3 different d18O scenarios, based on multiple possible dP values
for (i in 1:length(dP)) {
  TA = TA+273.15
  TL = TL+273.15
  esa = 6.108*(2.71828^((17.27*TA)/(TA+237.7)))
  lnalpha = (0.35041*(10^6/TA^3)) - (1.6664*(10^3/TA^2)) + (6.7123*(1/TA)) - (7.685*10^-3)
  alpha = exp(1)^lnalpha 
  h = RH*(TA/TL)
  epsilonP = 1000*(1-alpha)
  epsilonK = 14.2*(1-h)
  epsilon = epsilonK + epsilonP
  dA = dP[i] - epsilonP
  dE = ((alpha*dL) - (h*dA) - epsilon) / (1-h+(0.001*epsilonK))
  #groundwater as related to precip
  dG = c(dP[i]-1.5, dP[i], dP[i]+1.5)
  #calculate G based off ranges of dP and dG
  G1 = matrix(NA, length(dG),length(evap))
  G2 = matrix(NA, length(dG),length(evap))
  G3 = matrix(NA, length(dG),length(evap))
  for (g in 1:3) {
    for (e in 1:3) {
      G1[g,e] = ((precip[1]*(dL-dP[i])+evap[e]*(dE-dL))/(dG[g]-dL))
    }
  }
  for (g in 1:3) {
    for (e in 1:3) {
      G2[g,e] = ((precip[2]*(dL-dP[i])+evap[e]*(dE-dL))/(dG[g]-dL))
    }
  }
  for (g in 1:3) {
    for (e in 1:3) {
      G3[g,e] = ((precip[3]*(dL-dP[i])+evap[e]*(dE-dL))/(dG[g]-dL))
    }
  }
  G = cbind(G1,G2,G3)
  row.names(G) = dG
  colnames(G) = c(rep(precip[1],times = 3), rep(precip[2],times = 3),rep(precip[3],times = 3))
  G = as.data.frame(cbind(dG,G))
}
  library(ggplot2)
  library(scales)
  plot = ggplot()+
    geom_ribbon(aes(x=G[,1], ymax=G[,2], ymin = G[,4],colour = "PMax", fill = "Pmax"), size =1,alpha = 0.1, linetype = 2)+
    geom_line(aes(x=G[,1], y=G[,3], colour = "PMax"),size =1)+
    #geom_line(aes(x=G[,1], y=G[,4], colour = "PMax"), linetype = 2)+
    geom_ribbon(aes(x=G[,1], ymax=G[,5], ymin = G[,7],colour = "P", fill = "P"), size =1,alpha = 0.1, linetype = 2)+
    geom_line(aes(x=G[,1], y=G[,6], colour = "P"),size =1)+
    #geom_line(aes(x=G[,1], y=G[,7], colour = "P"), linetype = 2)+
    geom_ribbon(aes(x=G[,1], ymax=G[,8], ymin = G[,10],colour = "Pmin", fill = "Pmin"), size =1,alpha = 0.1, linetype = 2)+
    geom_line(aes(x=G[,1], y=G[,9], colour = "Pmin"),size =1)+
    #geom_line(aes(x=G[,1], y=G[,10], colour = "Pmin"), linetype = 2)+
    geom_hline(aes(yintercept = 0), colour = "black")+
    geom_hline(aes(yintercept = 118000), linetype = "dotted", size = 1, colour = "black")+
    xlab("Groundwater d18O")+
    ylab("Groundwater inflow")+
    #scale_y_continuous(labels = comma, breaks = c(0,100000,200000,300000,400000,500000))+
    scale_color_discrete(name = "Precipitation rate", labels = c("P","PMax","PMin"))+
    coord_cartesian(ylim = c(0, 500000)) +
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.position = "none", legend.text = element_text(size = 18), legend.title = element_text(size = 18), panel.background = element_rect(fill = 'white', colour = 'black'))
  return(plot)
}
