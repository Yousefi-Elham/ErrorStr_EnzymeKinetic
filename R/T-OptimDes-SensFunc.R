rm(list=ls())
library(numDeriv)
library(rgl)

par(mfrow=c(1,1))
theta1nc=theta1nc=c(12.0125, 8.5359, 5.6638,0)
theta2.hatc=c(5.538842, 2.193293, 1.624387,1) #mean of est. parameters of the fitted competitive log model


#mean of est. parameters of the fitted models in the compound T criteria
theta25.hatc=c(8.347629,  5.641800,  2.736343,1) 
theta25.hatnc=c(15.159141, 15.706634, 10.127313,0) 


theta1c=c(6.0645, 3.2799, 3.3153,1)
theta2.hatnc=c(11.910634, 10.457768,  9.561559,0) #mean of est. parameters of the fitted non competitive log model


#design space
epsilon <- 0.02
x1prim<- seq(epsilon,30, length=60)
x2prim<- seq(0,60,length=60)
gridprim <- expand.grid(x1prim,x2prim)

#--------------------------------------------------------------------------------------------------------------
# Note that the following designs have been resulted from separate programs!
# T optimum designs, Table 5, Log case 

des.psi0<-matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE)
weight0<-c(0.0095,0.1402,0.3600,0.4903)

des.psi5<-matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE)
weight5<-c(0.1688,0.1818,0.3002,0.3492)

des.psi1<-matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE)
weight1<-c(0.2500,0.2502,0.2500,0.2498)

#--------------------------------------------------------------------------------------------------------------
#delta formulation in T-optimality definition (to be used in sensitivity function formulation)
DELTA=function(TETA1,TETA2,DESIGN,WEIGHT){
  ETA1nc=log((TETA1[1]*DESIGN[,1])/((TETA1[2]*(1+(DESIGN[,2]/TETA1[3])))+(DESIGN[,1]*(1+((DESIGN[,2]*(1-TETA1[4]))/TETA1[3]))))) 
  ETA2c=log((TETA2[1]*DESIGN[,1])/((TETA2[2]*(1+(DESIGN[,2]/TETA2[3])))+DESIGN[,1]*(1+(((1-TETA2[4])*DESIGN[,2])/TETA2[3])))) 
  ETA3nc.c=sum(WEIGHT*((ETA1nc-ETA2c)^2))
  ETA3nc.c
}


DELTA0=DELTA(theta1nc,theta2.hatc,des.psi0,weight0)
DELTA05=DELTA(theta1nc,theta25.hatc,des.psi5,weight5)
DELTA15=DELTA(theta1c,theta25.hatnc,des.psi5,weight5)
DELTA1=DELTA(theta1c,theta2.hatnc,des.psi1,weight1)

#--------------------------------------------------------------------------------------------------------------
# Sensitivity functions in T-optimality (presented in the Appendix)
psi0=function(xest11,xest22){
  xest1=xest11
  xest2=xest22
  eta11=log((theta1nc[1]*xest1)/ ((theta1nc[2]*(1+(xest2/theta1nc[3]))) + xest1*(1+(((1-theta1nc[4])*xest2)/theta1nc[3]))))
  eta22=log((theta2.hatc[1]*xest1)/ ((theta2.hatc[2]*(1+(xest2/theta2.hatc[3]))) + xest1*(1+(((1-theta2.hatc[4])*xest2)/theta2.hatc[3]))))
  eta33=(eta11-eta22)^2
  eta33/DELTA0
}

psi0.5=function(xest11,xest22){
  xest1=xest11
  xest2=xest22
  eta11nc=log((theta1nc[1]*xest1)/ ((theta1nc[2]*(1+(xest2/theta1nc[3]))) + xest1*(1+(((1-theta1nc[4])*xest2)/theta1nc[3])))) #11
  eta22c=log((theta25.hatc[1]*xest1)/ ((theta25.hatc[2]*(1+(xest2/theta25.hatc[3]))) + xest1*(1+(((1-theta25.hatc[4])*xest2)/theta25.hatc[3])))) #10
  eta33nc.c=(eta11nc-eta22c)^2
  
  eta11c=log((theta1c[1]*xest1)/ ((theta1c[2]*(1+(xest2/theta1c[3]))) + xest1*(1+(((1-theta1c[4])*xest2)/theta1c[3])))) #10
  eta22nc=log((theta25.hatnc[1]*xest1)/ ((theta25.hatnc[2]*(1+(xest2/theta25.hatnc[3]))) + xest1*(1+(((1-theta25.hatnc[4])*xest2)/theta25.hatnc[3])))) #11
  eta33c.nc=(eta11c-eta22nc)^2
  ((1-0.5)*(eta33nc.c/DELTA05)+(0.5)*(eta33c.nc/DELTA15))
  
}

psi1=function(xest11,xest22){
  xest1=xest11
  xest2=xest22
  eta111=log((theta1c[1]*xest1)/ ((theta1c[2]*(1+(xest2/theta1c[3]))) + xest1*(1+(((1-theta1c[4])*xest2)/theta1c[3]))))
  eta221=log((theta2.hatnc[1]*xest1)/ ((theta2.hatnc[2]*(1+(xest2/theta2.hatnc[3]))) + xest1*(1+(((1-theta2.hatnc[4])*xest2)/theta2.hatnc[3]))))
  eta331=(eta111-eta221)^2
  eta331/DELTA1 
}
#--------------------------------------------------------------------------------------------------------------
#3D plot of A1
d.psi0<-psi0(gridprim[,1],gridprim[,2])
max(d.psi0)
psi0(des.psi0[,1],des.psi0[,2])

z1=matrix(d.psi0,60,60)
pmat1=persp(x1prim,x2prim,z1,xlab = "x1", ylab = "x2", zlab ="Psi",main = expression(A[1]),
            theta=120,phi=30, expand =0.4,
            shade=0.2,box = TRUE,col="gray",ticktype="detailed")
points(trans3d(des.psi0[,1],des.psi0[,2],max(d.psi0),pmat =pmat1), col = 2, pch = 16,cex=1.3)

#--------------------------------------------------------------------------------------------------------------
#3D plot of A2
d.psi5<-psi0.5(gridprim[,1],gridprim[,2])
max(d.psi5)
psi0.5(des.psi5[,1],des.psi5[,2])


z2=matrix(d.psi5,60,60)
pmat2=persp(x1prim,x2prim,z2,xlab = "x1", ylab = "x2", zlab ="Psi",main = expression(A[2]),
            theta=120,phi=30, expand =0.4,
            shade=0.2,box = TRUE,col="gray",ticktype="detailed")
points(trans3d(des.psi0[,1],des.psi0[,2],max(d.psi5),pmat =pmat2), col = 2, pch = 16,cex=1.3)

#--------------------------------------------------------------------------------------------------------------
#3D plot of A4
d.psi1<-psi1(gridprim[,1],gridprim[,2])
max(d.psi1)
psi1(des.psi1[,1],des.psi1[,2])

z3=matrix(d.psi1,60,60)
pmat3=persp(x1prim,x2prim,z3,xlab = "x1", ylab = "x2", zlab ="Psi",main = expression(A[4]),
            theta=120,phi=30, expand =0.4,
            shade=0.2,box = TRUE,col="gray",ticktype="detailed")
points(trans3d(des.psi0[,1],des.psi0[,2],max(d.psi1),pmat =pmat3), col = 2, pch = 16,cex=1.3)

