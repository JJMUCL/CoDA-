##################################################################################################################
#
# LR_ILR: ilr transformation for linear regression
# based on Hron et al 2009
#
# Author S.Chastin
# May 2017
#
# Please cite 
# Chastin SFM, Palarea-Albaladejo J, Dontje ML, Skelton DA (2015) Combined Effects of Time Spent in Physical Activity, Sedentary Behaviors and Sleep on Obesity and Cardio-Metabolic Health Markers: A Novel Compositional Data Analysis Approach. PLoS ONE 10(10): e0139984. https://doi.org/10.1371/journal.pone.0139984
# AND Hron et al. 2009 Journal of Applied Statistics
# and appropriate reference therein and within these codes
# if using this or part of this code.
#
# compo is a four part composition 
# base_rot is the basis rotation  1 compo[,1],compo[,2],compo[,3],compo[,4]
#                                 2 compo[,2],compo[,3],compo[,4],compo[,1]
#                                 etc....
##################################################################################################################


LR_ILR <- function(compo,base_rot=1){
D<-4 #four part composition  
### rotate through the basis
  #for (j in 1:ncol(X)) {
  if (base_rot ==1){
    data<-compo
  }
  if (base_rot ==2){
    data<-cbind(compo[,2:4],compo[,1])
  
  }
  if (base_rot ==3){
    data<-cbind(compo[,3:4],compo[,1:2])
  }
  if (base_rot ==4){
    data<-cbind(compo[,4],compo[1:3])
  }
#print(data)
z<-data[,1:3] #just to initialise
# z1
z[,1]<-sqrt(D-1/(D-1+1))*log(data[,1]/(data[,2]*data[,3]*data[,4])^(1/(D-1)))
#Z2
z[,2]<-sqrt(D-2/(D-2+1))*log(data[,2]/(data[,3]*data[,4])^(1/(D-2)))
#z3
z[,3]<-sqrt(D-3/(D-3+1))*log(data[,3]/(data[,4])^(1/(D-3)))
#names(Z)<-c("Z1","Z2","Z3")

      return(z)    
}
  