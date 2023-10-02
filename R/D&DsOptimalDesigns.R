# D and Ds optimal designs (Table 3)

rm(list=ls())
library(OptimalDesign)
library(numDeriv)
library(boot)
#-----------------------------------------------------------------------------------------------------------------------
#                                                        starting values
#-----------------------------------------------------------------------------------------------------------------------
priortheta.nc=c(8.6957,8.0664,12.0566,0) # lambda=0, prior par. values for non competitive model
priortheta.c=c(7.2976,4.3860,2.5821,1)   # lambda=1, prior par. values for competitive model
priortheta.ds.comb=c(7.4253,4.6808,3.0581,0.9636)  # prior par. values for encompassing model


# design space (standard case)
x1=seq(0,30,length=225)
x2=seq(0,60,length=270)
grid=expand.grid(x1,x2) 

#vector of weights
weight.4D=rep(0.25,4)
weight.3D=rep(1/3,3)

# c.v is a vector used for computation of c-optimal designs for a linear combination of parameters of interest
# in our case parameter of interest is the last parameter of the encompassing model
c.v=c(0,0,0,1)
#---------------------------------------------------------------------------------------------------------------------
priortheta.lognc=c(12.0125, 8.5359, 5.6638,0) # lambda=0, prior par. values for non competitive log model
priortheta.logc=c(6.0645, 3.2799, 3.3153,1)   # lambda=1, prior par. values for competitive log model
priortheta.ds.logcomb8=c(6.9897, 3.9799, 3.7380 , 0.8737) # prior par. values for encompassing log model

# design space (log case)
epsilon <- 0.02
x1prim<- c(epsilon,1:30)
x2prim<- 0:60
gridprim=expand.grid(x1prim,x2prim)
#-----------------------------------------------------------------------------------------------------------------------
#                                                         OPTIMAL DESIGNS FOR STANDARD MODELS
#-----------------------------------------------------------------------------------------------------------------------
#encompassing model
combd.model=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  theta4=par[4]
  
  eta=(theta1*xval[,1])/ ((theta2*(1+(xval[,2]/theta3)))+xval[,1]*(1+(((1-theta4)*xval[,2])/theta3))) 
}

#noncompetitive model
noncomp.model=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  
  eta=(theta1*xval[,1]) / ((theta2+xval[,1])*(1+(xval[,2]/theta3))) 
}

#competitive model
comp.model=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  
  eta=(theta1*xval[,1]) / ((theta2*(1+(xval[,2]/theta3)))+xval[,1]) 
}

#-----------------------------------------------------------------------------------------------------------------------
# D-optimal designs

# 4point D-optimal design for non competitive model (Denoted by 4D_N in Table 3)
F.4nc=jacobian(func=combd.model, x=priortheta.nc, method="Richardson",xval=grid)
D.4nc <- od_REX(as.matrix(F.4nc), crit = "D",alg.AA="REX", t.max = 50, eff = 1-1e-16)
supportD.4nc=grid[D.4nc$w.best!=0,]
Dopt.4nc=cbind(supportD.4nc,weight.4D)
Dopt.4nc
#-----------------------------------------

# 4point D-optimal design for competitive model (Denoted by 4D_C in Table 3)
F.4c=jacobian(func=combd.model, x=priortheta.c, method="Richardson",xval=grid)
D.4c <- od_REX(as.matrix(F.4c), crit = "D",alg.AA="REX", t.max = 50, eff = 1-1e-16)
supportD.4c=grid[D.4c$w.best!=0,]
Dopt.4c=cbind(supportD.4c,weight.4D)
Dopt.4c
#-----------------------------------------

# 3point D-optimal design for non-competitive model (Not presented in Table 3, but discussed in the Paper)
F.3nc=jacobian(func=noncomp.model, x=priortheta.nc[1:3], method="Richardson",xval=grid)
D.3nc <- od_REX(as.matrix(F.3nc), crit = "D",alg.AA="REX", t.max = 50, eff = 1-1e-16)
supportD.3nc=grid[D.3nc$w.best!=0,]
Dopt.3nc=cbind(supportD.3nc,weight.3D)
Dopt.3nc
#-----------------------------------------

# 3point D-optimal design for competitive model (Not presented in Table 3, but discussed in the Paper)
F.3c=jacobian(func=comp.model, x=priortheta.c[1:3], method="Richardson",xval=grid)
D.3c <- od_REX(as.matrix(F.3c), crit = "D",alg.AA="REX", t.max = 50, eff = 1-1e-16)
supportD.3c=grid[D.3c$w.best!=0,]
Dopt.3c=cbind(supportD.3c,weight.3D)
Dopt.3c
#-----------------------------------------

# 4point D-optimal design for encompassing model (Denoted by 4D_E in the Table 3)
F.4comb=jacobian(func=combd.model, x=priortheta.ds.comb, method="Richardson",xval=grid)
D.4comb <- od_REX(as.matrix(F.4comb), crit = "D",alg.AA="REX", t.max = 50, eff = 1-1e-16)
supportD.4comb=grid[D.4comb$w.best!=0,]
Dopt.4comb=cbind(supportD.4comb,weight.4D)
Dopt.4comb
#--------------------------------------------------------------------------------------------------------------
# D-s optimal designs

# optimal weights is computed through simplex method of LP (Harman R., Jurik T.(2008))

# function to calculate Ds-optimal designs
Ds.weightfun=function(F,cv)
{
  
  n=dim(F)[1]
  m=dim(F)[2]
  A3=cbind(t(F),-t(F),-cv)
  A3=rbind(A3,c(rep(1,2*n),0))
  b3=c(rep(0,m),1)
  a=c(rep(0,2*n),1)
  w.sim=simplex(a=a, A3 = A3, b3 = b3, maxi = TRUE)
  w.vec=w.sim$soln[1:((2*n))]   
  w=w.vec[1:n]+w.vec[(n+1):(2*n)]
  num.w=which(w!=0)
  list(weight = as.vector(w[num.w]), indicator = num.w, optvalue = w.sim$value)
  
}
#-----------------------------------------

# D-s optimal design for non competitive model (Denoted by Ds_N in Table 3)
Ds.nc=Ds.weightfun(F.4nc,c.v)
supportDs.nc=grid[Ds.nc$indicator,]
weightDs.nc=Ds.nc$weight
Dsopt.nc=cbind(supportDs.nc,weightDs.nc)
Dsopt.nc
#-----------------------------------------

# D-s optimal design for competitive model (Denoted by Ds_C in Table 3)
Ds.c=Ds.weightfun(F.4c,c.v)
supportDs.c=grid[Ds.c$indicator,]
weightDs.c=Ds.c$weight
Dsopt.c=cbind(supportDs.c,weightDs.c)
Dsopt.c
#-----------------------------------------

# D-s optimal design for model discrimination using the Encompassing model.
# It is presented in the Standard case of Table 5, denoted by A3.
F.Dscomb=jacobian(func=combd.model, x=priortheta.ds.comb, method="Richardson",xval=grid)
Ds.comb=Ds.weightfun(F.Dscomb,c.v)
supportDs.comb=grid[Ds.comb$indicator,]
weightDs.comb=Ds.comb$weight
Dsopt.comb=cbind(supportDs.comb,weightDs.comb)
Dsopt.comb
#-----------------------------------------------------------------------------------------------------------------------
#                                                       OPTIMAL DESIGNS FOR THE LOG-MODELS 
#-----------------------------------------------------------------------------------------------------------------------
#encompassing log model
combd.logmodel=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  theta4=par[4]
  
  eta=log((theta1*xval[,1])/ ((theta2*(1+(xval[,2]/theta3)))+xval[,1]*(1+(((1-theta4)*xval[,2])/theta3))))
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

#---------------------------------------------------------------------------------------------------
# D-optimal designs for log models


# 4point D-optimal design for non competitive log model (Denoted by 4D_N in the Log case of Table 3)
F.4lognc=jacobian(func=combd.logmodel, x=priortheta.lognc, method="Richardson",xval=gridprim)
D.4lognc <- od_REX(as.matrix(F.4lognc), crit = "D",alg.AA="REX", t.max = 50, eff = 1-1e-16)
supportD.4lognc=gridprim[D.4lognc$w.best!=0,]
Dopt.4lognc=cbind(supportD.4lognc,weight.4D)
Dopt.4lognc

#-----------------------------------------
# 4point D-optimal design for competitive log model (Denoted by 4D_C in the Log case of Table 3)
F.4logc=jacobian(func=combd.logmodel, x=priortheta.logc, method="Richardson",xval=gridprim)
D.4logc <- od_REX(as.matrix(F.4logc), crit = "D",alg.AA="REX", t.max = 50, eff = 1-1e-16)
supportD.4logc=gridprim[D.4logc$w.best!=0,]
Dopt.4logc=cbind(supportD.4logc,weight.4D)
Dopt.4logc
#-----------------------------------------
# D-optimal designs for non competitive log model (Denoted by 3D_N in the Log case of Table 3)
# linearized at first three respective prior values
# it has 4 points of support according to computation below!
F.3lognc=jacobian(func=noncomp.logmodel, x=priortheta.lognc[1:3], method="Richardson",xval=gridprim)
D.3lognc <- od_REX(as.matrix(F.3lognc), crit = "D",alg.AA="REX", t.max = 50, eff = 1-1e-16)
supportD.3lognc=gridprim[D.3lognc$w.best!=0,]
Dopt.3lognc=cbind(supportD.3lognc,weight.4D)
Dopt.3lognc

#-----------------------------------------
# 3point D-optimal designs for competitive log model (Denoted by 3D_C in the Log case of Table 3)
F.3logc=jacobian(func=comp.logmodel, x=priortheta.logc[1:3], method="Richardson",xval=gridprim)
D.3logc <- od_REX(as.matrix(F.3logc), crit = "D",alg.AA="REX", t.max = 50, eff = 1-1e-16)
supportD.3logc=gridprim[D.3logc$w.best!=0,]
Dopt.3logc=cbind(supportD.3logc,weight.3D)
Dopt.3logc

#-----------------------------------------
# 4point D-optimal designs for encompassing log model (Denoted by 4D_E in the Log case of Table 3)
F.4logcombd=jacobian(func=combd.logmodel, x=priortheta.ds.logcomb8, method="Richardson",xval=gridprim)
D.4logcombd <- od_REX(as.matrix(F.4logcombd), crit = "D",alg.AA="REX", t.max = 50, eff = 1-1e-16)
supportD.4logcombd=gridprim[D.4logcombd$w.best!=0,]
Dopt.4logcombd=cbind(supportD.4logcombd,weight.4D)
Dopt.4logcombd

#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#robustness check of D-optimal designs w.r.t par. prior values (for the encompassing) [NOT presented in the paper]
#all the designs are the same on the corners as before
# while the prior values could be any values from the uniform distribution

F.4logcombd=jacobian(func=combd.logmodel, x=c(runif(3, min = 0, max = 60),runif(1,min=0,max=1)), method="Richardson",xval=gridprim)
D.4logcombd <- od_REX(as.matrix(F.4logcombd), crit = "D",alg.AA="REX", t.max = 50, eff = 1-1e-16)
supportD.4logcombd=gridprim[D.4logcombd$w.best!=0,]
Dopt.4logcombd=cbind(supportD.4logcombd,weight.4D)
Dopt.4logcombd

##################################################################################################
# D-s optimal designs for log models

#D-s optimal designs for non competitive log model (Denoted by Ds_N in the Log case of Table 3)
Ds.lognc=Ds.weightfun(F.4lognc,c.v)
supportDs.lognc=gridprim[Ds.lognc$indicator,]
weightDs.lognc=Ds.lognc$weight
Dsopt.lognc=cbind(supportDs.lognc,weightDs.lognc)
Dsopt.lognc

#D-s optimal designs for competitive log model (Denoted by Ds_C in the Log case of Table 3)
Ds.logc=Ds.weightfun(F.4logc,c.v)
supportDs.logc=gridprim[Ds.logc$indicator,]
weightDs.logc=Ds.logc$weight
Dsopt.logc=cbind(supportDs.logc,weightDs.logc)
Dsopt.logc

##-------------------------------------------------------------------------------------------------
# D-s optimal design for model discrimination using the Encompassing log model.
# It is presented in the Log case of Table 5, denoted by A3.
Ds.logcomb.fun=function(prior){
  F.Ds.logcomb=jacobian(func=combd.logmodel, x=prior, method="Richardson",xval=gridprim)
  Ds.logcomb=Ds.weightfun(F.Ds.logcomb,c.v)
  supportDs.logcomb=gridprim[Ds.logcomb$indicator,]
  weightDs.logcomb=Ds.logcomb$weight
  Dsopt.logcomb=cbind(supportDs.logcomb,weightDs.logcomb)
  Dsopt.logcomb
}

Ds.logcomb.fun(priortheta.ds.logcomb8)
#-------------------------------------------------------------------------------------------------
# Robustness check of Ds optimum design w.r.t par. prior values (for the encompassing) [NOT presented in the paper] 
# shows that the support points remain the same
# while the optimum weights do not
Ds.logcomb.fun(c(runif(3, min = 0, max = 60),runif(1,min=0,max=1)))
