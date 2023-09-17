rm(list=ls())

#parameter values (standard case)
theta1nc=c(8.6957,8.0664,12.0566,0) 

THETAHATC<-c(8.680591, 8.0018861, 3.24841556) # estimates derived from the minimization procedure in Toptim_nu0_standard
lambdac=1   



theta1c=c(7.2976,4.3860,2.5821,1) 

THETAHATNC<-c(6.193258,  0.001000,  9.099974) # estimates derived from the minimization procedure in Toptim_nu1_standard
lambdanc=0   


#--------------------------------------------------------------------------------------------------------------
# Note that the following designs have been resulted from separate programs!
# T optimum designs, Table 5, Standard case 

dest1<-matrix(c(30.00,0,3.214,0,30.00,21.413,5.625,11.152),nrow=4,byrow=TRUE)
wt1<-c(0.063,0.063,0.310,0.564)


dest2<-matrix(c(30.00,0,3.348,0,30.00,22.082,5.759,11.375),nrow=4,byrow=TRUE)
wt2<-c(0.058,0.189,0.260,0.493)


dest3<-matrix(c(30.00,0,2.678,0,30.00,25.874,6.562,8.922),nrow=4,byrow=TRUE)
wt3<-c(0.052,0.137,0.330,0.481)


dest4<-matrix(c(30.00,0,3.080,0,30.00,22.751,5.491,11.598),nrow=4,byrow=TRUE)
wt4<-c(0.060,0.250,0.250,0.440)
#--------------------------------------------------------------------------------------------------------------

#delta_0 (lack of fit sum of squares when the non competitive model is fixed)
DELTA0=function(TETA,XXX1,XXX2,WW){
  ETA1=(theta1nc[1]*XXX1)/((theta1nc[2]*(1+(XXX2/theta1nc[3])))+(XXX1*(1+((XXX2*(1-theta1nc[4]))/theta1nc[3]))))
  ETA2=(TETA[1]*XXX1)/((TETA[2]*(1+(XXX2/TETA[3])))+(XXX1*(1+((XXX2*(1-lambdac))/TETA[3]))))
  sum(WW*((ETA1-ETA2)^2))
}

#delta_1 (lack of fit sum of squares when the competitive model is fixed)
DELTA1=function(THETA,XX1,XX2,WW1){
  ETA1=(theta1c[1]*XX1)/((theta1c[2]*(1+(XX2/theta1c[3])))+(XX1*(1+((XX2*(1-theta1c[4]))/theta1c[3]))))
  ETA2=(THETA[1]*XX1)/((THETA[2]*(1+(XX2/THETA[3])))+(XX1*(1+((XX2*(1-lambdanc))/THETA[3]))))
  sum(WW1*((ETA1-ETA2)^2))
}

#--------------------------------------------------------------------------------------------------------------
#T_efficiency formulation

g1<-DELTA0(THETAHATC,dest2[,1],dest2[,2],wt2)
g2<-DELTA0(THETAHATC,dest3[,1],dest3[,2],wt3)
g3<-DELTA0(THETAHATC,dest4[,1],dest4[,2],wt4)
gg<-DELTA0(THETAHATC,dest1[,1],dest1[,2],wt1)
T11<-gg/gg
T21<-g1/gg
T31<-g2/gg
T41<-g3/gg
T1<-c(T11,T21,T31,T41)
T1


g11<-DELTA1(THETAHATNC,dest1[,1],dest1[,2],wt1)
g22<-DELTA1(THETAHATNC,dest2[,1],dest2[,2],wt2)
g33<-DELTA1(THETAHATNC,dest3[,1],dest3[,2],wt3)
ggg<-DELTA1(THETAHATNC,dest4[,1],dest4[,2],wt4)
T12<-g11/ggg
T22<-g22/ggg
T32<-g33/ggg
T42<-ggg/ggg
T2<-c(T12,T22,T32,T42)
T2

round(cbind(T1,T2)*100, digits=2)
