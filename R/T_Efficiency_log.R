rm(list=ls())

#parameter values (log case)
theta1nc=c(12.0125, 8.5359, 5.6638,0) 

MEANTH.HATC<-c(8.025290, 5.520947, 3.575972) # estimates derived from the minimization procedure in Toptim_A1_log
lambdac=1   


theta1c=c(6.0645, 3.2799, 3.3153,1) 

MEANTH.HATNC<-c(14.34820, 13.86322, 8.715987) # estimates derived from the minimization procedure in Toptim_A4_log
lambdanc=0   


#--------------------------------------------------------------------------------------------------------------
# Note that the following designs have been resulted from separate programs!
# T optimum designs, Table 5, Log case 

dest1<-matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE)
wt1<-c(0.0095,0.1402,0.36,0.4903)

dest2<-matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE)
wt2<-c(0.1688,0.1818,0.3002,0.3492)

dest3<-matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE)
wt3<-c(0.1634,0.2189,0.2811,0.3366)

dest4<-matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE)
wt4<-c(0.2500,0.2502,0.2500,0.2498)
#--------------------------------------------------------------------------------------------------------------
#delta_0 (lack of fit sum of squares of the optimum designs, when the non competitive log model is fixed)

DELTA0=function(TETA,XXX1,XXX2,WW){   
  ETA10=log((theta1nc[1]*XXX1)/((theta1nc[2]*(1+(XXX2/theta1nc[3])))+(XXX1*(1+((XXX2*(1-theta1nc[4]))/theta1nc[3])))))
  ETA20=log((TETA[1]*XXX1)/((TETA[2]*(1+(XXX2/TETA[3])))+(XXX1*(1+((XXX2*(1-lambdac))/TETA[3])))))
  sum(WW*((ETA10-ETA20)^2))
}


#delta_1 (lack of fit sum of squares of the optimum designs, when the competitive log model is fixed)
DELTA1=function(TETA,XX1,XX2,WW){   
  ETA11=log((theta1c[1]*XX1)/((theta1c[2]*(1+(XX2/theta1c[3])))+(XX1*(1+((XX2*(1-theta1c[4]))/theta1c[3])))))
  ETA21=log((TETA[1]*XX1)/((TETA[2]*(1+(XX2/TETA[3])))+(XX1*(1+((XX2*(1-lambdanc))/TETA[3])))))
  sum(WW*((ETA11-ETA21)^2))
}

#--------------------------------------------------------------------------------------------------------------
#T_efficiency formulation
Teff<-function(d1,d2){
  round((d1/d2)*100, 2)
}

del00<-DELTA0(MEANTH.HATC,dest1[,1],dest1[,2],wt1)
del20<-DELTA0(MEANTH.HATC,dest2[,1],dest2[,2],wt2)
del30<-DELTA0(MEANTH.HATC,dest3[,1],dest3[,2],wt3)
del40<-DELTA0(MEANTH.HATC,dest4[,1],dest4[,2],wt4)

Tef00<-Teff(del00,del00)
Tef20<-Teff(del20,del00)
Tef30<-Teff(del30,del00)
Tef40<-Teff(del40,del00)
TEF0<-c(Tef00,Tef20,Tef30,Tef40)

del21<-DELTA1(MEANTH.HATNC,dest1[,1],dest1[,2],wt1)
del31<-DELTA1(MEANTH.HATNC,dest2[,1],dest2[,2],wt2)
del41<-DELTA1(MEANTH.HATNC,dest3[,1],dest3[,2],wt3)
del11<-DELTA1(MEANTH.HATNC,dest4[,1],dest4[,2],wt4)

Tef21<-Teff(del21,del11)
Tef31<-Teff(del31,del11)
Tef41<-Teff(del41,del11)
Tef11<-Teff(del11,del11)
TEF1<-c(Tef21,Tef31,Tef41,Tef11)

TEF<-cbind(TEF0,TEF1)
TEF
