# computation of the T-optimum design A2, presented in the table 5 in the log case 

rm(list=ls())
library(mgcv)
par(mfrow=c(1,1))
#-----------------------------------------------------------------------------------------------------------------------
#                                                        starting values
#-----------------------------------------------------------------------------------------------------------------------
# the design space
epsilon <- 0.02
x1 <- c(epsilon,1:30)
x2 <- 0:60
grid <- expand.grid(x1,x2)


#starting design
xstart1 <- c(epsilon,30,epsilon,30,15)
xstart2 <- c(0,0,60,60,30)
xstart <- rbind(xstart1,xstart2)
n <- ncol(xstart) 


# starting weights
wstart <- c(rep((1/n),n))


# initial parameter values
theta1nc=c(12.0125, 8.5359, 5.6638,0) #nu=0, non-comp. log model is fixed
theta2.startc=c(6.0645, 3.2799, 3.3153)
lambdac=1   # lambdac=1, comp. log model is fitted


theta1c=c(6.0645, 3.2799, 3.3153,1) #nu=1, comp. log model is fixed
theta2.startnc=c(12.0125, 8.5359, 5.6638)
lambdanc=0   # lambdanc=0, non-comp. log model is fitted



#Maximum number of algorithm iterations
max.iter=5000


THETA.HAT=matrix(0,max.iter,6)
num.max=numeric(max.iter)
#-----------------------------------------------------------------------------------------------------------------------
#                                                        THETA.HAT
#-----------------------------------------------------------------------------------------------------------------------
iter=0

while(iter< max.iter){
  iter=iter+1
  plot(iter,xlim=c(iter-1,iter+1))
  text(x=iter+1,y=iter,iter,col="red")
  
  
  #first step optimization (minimization with respect to parameters of the fitted model)
  Tp.theta=function(par){
    theta21c=par[1]
    theta22c=par[2]
    theta23c=par[3]
    
    theta21nc=par[4]
    theta22nc=par[5]
    theta23nc=par[6]
    
    eta1nc=log((theta1nc[1]*xstart[1,])/((theta1nc[2]*(1+(xstart[2,]/theta1nc[3]))) + xstart[1,]*(1+(((1-theta1nc[4])*xstart[2,])/theta1nc[3])))) #11
    eta2c=log((theta21c*xstart[1,])/((theta22c*(1+(xstart[2,]/theta23c)))+xstart[1,]*(1+(((1-lambdac)*xstart[2,])/theta23c)))) #10
    eta3nc.c=wstart*(eta1nc-eta2c)^2
    
    eta1c=log((theta1c[1]*xstart[1,])/((theta1c[2]*(1+(xstart[2,]/theta1c[3]))) + xstart[1,]*(1+(((1-theta1c[4])*xstart[2,])/theta1c[3])))) #10
    eta2nc=log((theta21nc*xstart[1,])/((theta22nc*(1+(xstart[2,]/theta23nc)))+xstart[1,]*(1+(((1-lambdanc)*xstart[2,])/theta23nc)))) #11
    eta3c.nc=wstart*(eta1c-eta2nc)^2
    
    (1-0.5)*log(sum(eta3nc.c))+(0.5)*log(sum(eta3c.nc))
  }
  
  theta2.hat=optim(par=c(theta2.startc,theta2.startnc),Tp.theta,lower=c(epsilon,epsilon,epsilon,epsilon,epsilon,epsilon),upper=c(Inf,Inf,Inf,Inf,Inf,Inf),method="L-BFGS-B")$par
  
  
  theta21.hatc=theta2.hat[1]
  theta22.hatc=theta2.hat[2]
  theta23.hatc=theta2.hat[3]
  
  theta21.hatnc=theta2.hat[4]
  theta22.hatnc=theta2.hat[5]
  theta23.hatnc=theta2.hat[6]
  THETA.HAT[iter,]=theta2.hat
  
  theta2.startc=theta2.hat[1:3]
  theta2.startnc=theta2.hat[4:6]
  
  #-----------------------------------------------------------------------------------------------------------------------
  #                                                        X.HAT
  #-----------------------------------------------------------------------------------------------------------------------
  #delta formulation in T-optimality definition (lack of fit sum of squares related to the compound T criteria)
  DELTA=function(TETA1,TETA2){
    ETA1=log((TETA1[1]*xstart[1,])/((TETA1[2]*(1+(xstart[2,]/TETA1[3])))+(xstart[1,]*(1+((xstart[2,]*(1-TETA1[4]))/TETA1[3])))) )
    ETA2=log((TETA2[1]*xstart[1,])/((TETA2[2]*(1+(xstart[2,]/TETA2[3])))+xstart[1,]*(1+(((1-TETA2[4])*xstart[2,])/TETA2[3]))) )
    ETA3=sum(wstart*((ETA1-ETA2)^2))
    ETA3
  }
  
  DELTA1=DELTA(theta1nc,c(theta2.hat[1:3],lambdac))
  DELTA2=DELTA(theta1c,c(theta2.hat[4:6],lambdanc))
  #----------------------------------------------------
  
  #second step optimization (maximization with respect to design points)
  psi=function(xest11,xest22){
    xest1=xest11
    xest2=xest22
    eta11nc=log((theta1nc[1]*xest1)/ ((theta1nc[2]*(1+(xest2/theta1nc[3]))) + xest1*(1+(((1-theta1nc[4])*xest2)/theta1nc[3])))) #11
    eta22c=log((theta21.hatc*xest1)/ ((theta22.hatc*(1+(xest2/theta23.hatc))) + xest1*(1+(((1-lambdac)*xest2)/theta23.hatc)))) #10
    eta33nc.c=(eta11nc-eta22c)^2
    
    eta11c=log((theta1c[1]*xest1)/ ((theta1c[2]*(1+(xest2/theta1c[3]))) + xest1*(1+(((1-theta1c[4])*xest2)/theta1c[3])))) #10
    eta22nc=log((theta21.hatnc*xest1)/ ((theta22.hatnc*(1+(xest2/theta23.hatnc))) + xest1*(1+(((1-lambdanc)*xest2)/theta23.hatnc)))) #11
    eta33c.nc=(eta11c-eta22nc)^2
    (1-0.5)*((eta33nc.c)/DELTA1)+(0.5)*((eta33c.nc)/DELTA2)
    
  }
  psi.vec=psi(grid[,1],grid[,2])
  num.max[iter]=which.max(psi.vec)
  x.hat=c(grid[num.max[iter],1],grid[num.max[iter],2])
  
  xstart=cbind(xstart,x.hat)
  alpha=1/(iter+1)
  wstart=c(((1-alpha)*wstart),alpha)
  design=rbind(xstart,wstart)
  n=ncol(xstart)
  
}

#-----------------------------------------------------------------------------------------------
# some information fur further use

THETA.HAT[1:10,]

#mean of parameters of the fitted competitive log model in 3d-plots
TOTALMEAN=apply(THETA.HAT,2,mean) 
TOTALMEAN
#----------------------------------------------------
#optimum designs
design1.uniq=uniquecombs(t(design[-3,]))
count1=attr(design1.uniq,"index")
count2=unique(sort(count1))
weight=numeric(nrow(design1.uniq))
for(j in 1:length(count2)){
  m=which(count1==count2[j])
  weight[j]=sum(t(design)[m,3])
}

design.opt=cbind(design1.uniq,weight)
design.opt


# here, support points with considerable weights are chosen
mm=which(design.opt[,3]>0.001)
design.opt[mm,]


#> TOTALMEAN
#[1]  8.347629  5.641800  2.736343 15.159141 15.706634 10.127313

#> MEAN.THETA.HAT
#[,1]     [,2]     [,3]     [,4]     [,5]      [,6]
#[1,] 9.796334 6.801346 2.623006 12.19546 10.77184  9.888290
#[2,] 8.025490 5.521164 3.576009 14.34868 13.86346  8.715931
#[3,] 6.028741 4.125426 3.828016 18.05085 19.76623  8.998280
#[4,] 6.931412 4.579335 3.499738 16.57495 20.76304 11.360682
#[5,] 7.519792 5.029093 3.365096 15.65553 18.71617 11.016862
#[6,] 7.880233 5.303097 3.180910 15.25751 17.66152 11.020484
#[7,] 7.904181 5.323609 3.139802 15.29325 17.37939 10.814678
#[8,] 7.878915 5.306903 3.062082 15.46916 17.39685 10.738872
#[9,] 8.079909 5.442669 2.990207 15.13671 17.10506 10.964576
#[10,] 8.245622 5.569956 2.938399 14.90225 16.56420 10.920826
########################################################
#########################
#########################
#> design.opt[mm,]
#xstart1 xstart2    weight
#0.02       0 0.1688062 #0.1688
#30.00       0 0.1818036 #0.1818
#0.02      60 0.3001800 #0.3002
#30.00      60 0.3491702 #0.3492
#> sum(design.opt[mm,3])
#[1] 0.99996
#> DELTA.NAHAEI1
#[1] 0.2092339