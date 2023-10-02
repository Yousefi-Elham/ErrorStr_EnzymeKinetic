# plot of the sensitivity functions related to D and DS optimal designs (in the log case)

rm(list=ls())
library(numDeriv)
library(rgl)


#Prior parameter values
priortheta.lognc=c(12.0125, 8.5359, 5.6638,0) #prior values for non competitive log model
priortheta.logc=c(6.0645, 3.2799, 3.3153,1)  #prior values for competitive log model
priortheta.ds.logcomb8=c(6.9897, 3.9799, 3.7380 , 0.8737) #prior for lambda=0.8


# the design space
epsilon <- 0.02
x1prim<- seq(epsilon,30, length=60)
x2prim<- seq(0,60,length=60)
gridprim <- expand.grid(x1prim,x2prim)

#--------------------------------------------------------------------------------------------------
# log model functions

#encompassing log model
combd.logmodel=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  theta4=par[4]
  
  eta=log((theta1*xval[,1])/ ((theta2*(1+(xval[,2]/theta3)))+xval[,1]*(1+(((1-theta4)*xval[,2])/theta3)))) 
}


#encompassing log model when the last parameter is fixed
combd.logmodel3=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  
  eta=log((theta1*xval[,1])/ ((theta2*(1+(xval[,2]/theta3)))+xval[,1]*(1+(((1-0.8737)*xval[,2])/theta3)))) 
}


#non competitive log model
noncomp.logmodel=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  
  eta=log((theta1*xval[,1]) / ((theta2+xval[,1])*(1+(xval[,2]/theta3)))) 
}


#competitive log model
comp.logmodel=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  
  eta=log((theta1*xval[,1]) / ((theta2*(1+(xval[,2]/theta3)))+xval[,1])) 
}
#--------------------------------------------------------------------------------------------------
# D and Ds optimal designs (resulting from D&DsOptimalDesigns R script)
# will be used here to plot sensitivity functions

des.D4<-matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE)
weight.4D<-rep(0.25,4)


des.D3<-matrix(c(0.02,0,30.00,0,0.02,60),nrow=3,byrow=TRUE)
weight.3D<-rep(1/3,3)


DsC<-matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE)
wsC<-c(0.017,0.173,0.327,0.483)


#Ds-optimal design for discrimination
dest3<-matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE)
wt3<-c(0.1633,0.2189,0.2811,0.3367)

#--------------------------------------------------------------------------------------------------
# D-optimality sensitivity function: used to check if a calculated design is actually D-optimal
D2.func<-function(model,prior,des.x,des.w,GRID2){
  
  D1.func<-function(GRID1){
    F.mat1<-jacobian(func=model, x=prior, method="Richardson",xval=des.x)
    w.mat<-diag(des.w,nrow=length(des.w))
    M<-t(F.mat1)%*%w.mat%*%F.mat1
    
    F.mat2<-jacobian(func=model, x=prior, method="Richardson",xval=GRID1)
    D=F.mat2%*%solve(M)%*%t(F.mat2)
    D
  }
  d.output<-c()
  for(i in 1:nrow(GRID2)){
    d.output[i]<-D1.func(GRID2[i,])
  }
  d.output
  
}


#--------------------------------------------------------------------------------------------------
# Ds-optimality sensitivity function: used to check if a calculated design is actually Ds-optimal
D2.func1 <- function(model,model1,prior,des.x,des.w,GRID2){
  
  D1.func<-function(GRID1){
    F.mat1 <- jacobian(func=model, x=prior, method="Richardson",xval=des.x)
    w.mat <- diag(des.w,nrow=length(des.w))
    M <- t(F.mat1)%*%w.mat%*%F.mat1
    
    F.mat2<-jacobian(func=model, x=prior, method="Richardson",xval=GRID1)
    D=F.mat2%*%solve(M)%*%t(F.mat2)
    D
  }
  
  D1star.func <- function(GRID1){
    F.mat1star <- jacobian(func=model1, x=prior[1:3], method="Richardson",xval=des.x)
    w.mat <- diag(des.w,nrow=length(des.w))
    M22 <- t(F.mat1star)%*%w.mat%*%F.mat1star
    
    F.mat2star <- jacobian(func=model1, x=prior[1:3], method="Richardson",xval=GRID1)
    Dstar=F.mat2star%*%solve(M22)%*%t(F.mat2star)
    Dstar
  }
  ds.output <- c()
  for(i in 1:nrow(GRID2)){
    ds.output[i]<-D1.func(GRID2[i,])-D1star.func(GRID2[i,])
  }
  ds.output
  
}
#--------------------------------------------------------------------------------------------------
#                           PLOT OF THE SENSITIVITY FUNCTIONS (PRESENTED IN THE APPENDIX)
#--------------------------------------------------------------------------------------------------
d.4NC <- D2.func(combd.logmodel,priortheta.lognc,des.D4,weight.4D,gridprim)
max(d.4NC)

z=matrix(d.4NC,60,60)
pmat1=persp(x1prim,x2prim,z,xlab = "x1", ylab = "x2", zlab =" Sensitivity function",main = "4DN",
            theta=120,phi=30, expand =0.4,
            shade=0.2,box = TRUE,col="gray",ticktype="detailed")
points(trans3d(des.D4[,1],des.D4[,2],max(d.4NC),pmat =pmat1), col = 2, pch = 16,cex=1.3)

#--------------------------------------------------------------------------------------------------
d.3NC<-D2.func(noncomp.logmodel,priortheta.lognc[1:3],des.D4,weight.4D,gridprim)
max(d.3NC)

z=matrix(d.3NC,60,60)
pmat2=persp(x1prim,x2prim,z,xlab = "x1", ylab = "x2", zlab =" Sensitivity function",main = "3DN",
            theta=120,phi=30, expand =0.4,
            shade=0.2,box = TRUE,col="gray",ticktype="detailed")
points(trans3d(des.D4[,1],des.D4[,2],max(d.3NC),pmat =pmat2), col = 2, pch = 16,cex=1.3)

#--------------------------------------------------------------------------------------------------
d.4C<-D2.func(combd.logmodel,priortheta.logc,des.D4,weight.4D,gridprim)
max(d.4C)

z=matrix(d.4C,60,60)
pmat3=persp(x1prim,x2prim,z,xlab = "x1", ylab = "x2", zlab =" Sensitivity function",main = "4DC",
            theta=120,phi=30, expand =0.4,
            shade=0.2,box = TRUE,col="gray",ticktype="detailed")
points(trans3d(des.D4[,1],des.D4[,2],max(d.4C),pmat =pmat3), col = 2, pch = 16,cex=1.3)

#--------------------------------------------------------------------------------------------------
d.3C<-D2.func(comp.logmodel,priortheta.logc[1:3],des.D3,weight.3D,gridprim)
max(d.3C)

z=matrix(d.3C,60,60)
pmat4=persp(x1prim,x2prim,z,xlab = "x1", ylab = "x2", zlab =" Sensitivity function",main = "3DC",
            theta=120,phi=30, expand =0.4,
            shade=0.2,box = TRUE,col="gray",ticktype="detailed")
points(trans3d(des.D3[,1],des.D3[,2],max(d.3C),pmat =pmat4), col = 2, pch = 16,cex=1.3)

#--------------------------------------------------------------------------------------------------
d.4Comb<-D2.func(combd.logmodel,priortheta.ds.logcomb8,des.D4,weight.4D,gridprim)
max(d.4Comb)

z=matrix(d.4Comb,60,60)
pmat5=persp(x1prim,x2prim,z,xlab = "x1", ylab = "x2", zlab =" Sensitivity function",main = "4DE",
            theta=120,phi=30, expand =0.4,
            shade=0.2,box = TRUE,col="gray",ticktype="detailed")
points(trans3d(des.D4[,1],des.D4[,2],max(d.4Comb),pmat =pmat5), col = 2, pch = 16,cex=1.3)

#--------------------------------------------------------------------------------------------------
dsNC<-D2.func1(combd.logmodel,noncomp.logmodel,priortheta.lognc,des.D4,weight.4D,gridprim)
max(dsNC)

z=matrix(dsNC,60,60)
pmat6=persp(x1prim,x2prim,z,xlab = "x1", ylab = "x2", zlab ="Sensitivity func.",main = "DsN",
            theta=120,phi=30, expand =0.3,
            shade=0.2,box = TRUE,col="gray",ticktype="detailed")
points(trans3d(des.D4[,1],des.D4[,2],1,pmat =pmat6), col = 2, pch = 16,cex=1.3)


#--------------------------------------------------------------------------------------------------
dsC<-D2.func1(combd.logmodel,comp.logmodel,priortheta.logc,DsC,wsC,gridprim)
max(dsC)

z2=matrix(dsC,60,60)
pmat7=persp(x1prim,x2prim,z2,xlab = "x1", ylab = "x2", zlab ="Sensitivity func.",,main = "DsC",
            theta=120,phi=30, expand =0.3,
            shade=0.2,box = TRUE,col="gray",ticktype="detailed")
points(trans3d(DsC[,1],DsC[,2],1,pmat =pmat7), col = 2, pch = 16,cex=1.3)

#--------------------------------------------------------------------------------------------------

ds.disc<-D2.func1(combd.logmodel,combd.logmodel3,priortheta.ds.logcomb8,dest3,wt3,gridprim)
max(ds.disc)

z3=matrix(ds.disc,60,60)
pmat8=persp(x1prim,x2prim,z3,xlab = "x1", ylab = "x2", zlab ="Sensitivity func.",,main = expression(A[3]),
            theta=120,phi=30, expand =0.4,
            shade=0.2,box = TRUE,col="gray",ticktype="detailed")
points(trans3d(dest3[,1],dest3[,2],1,pmat =pmat8), col = 2, pch = 16,cex=1.3)

