# computation of the T-optimum design \nu = 0, presented in the table 5 in the standard case 
rm(list=ls())
#-----------------------------------------------------------------------------------------------------------------------
#                                                        starting values
#-----------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,1))

# the design space
x1=seq(0,30,length=225)
x2=seq(0,60,length=270)
grid=expand.grid(x1,x2)


#starting design
xstart1=c(0,30,0,30,15)
xstart2=c(0,0,60,60,30)
xstart=rbind(xstart1,xstart2)
n=ncol(xstart)


# starting weights
wstart=c(rep((1/n),n)) 


# initial parameter values
theta1=c(8.6957,8.0664,12.0566,0) #nu=0, non competitive model is fixed

theta2.start=c(7.2976,4.3860,2.5821) #competitive model is fitted
lambda=1   


#Maximum number of algorithm iterations
max.iter=20000
efficiency <- numeric(max.iter)


THETA.HAT=matrix(0,max.iter,3)
num.max=numeric(max.iter)
#-----------------------------------------------------------------------------------------------------------------------
#                                                        THETA.HAT
#-----------------------------------------------------------------------------------------------------------------------
iter=0
delta.0=0
psi.0=1

while(iter< max.iter){
  
  iter=iter+1
  plot(iter,xlim=c(iter-1,iter+1))
  text(x=iter+1,y=iter,iter,col="red")
  
  
  #first step optimization (minimization with respect to parameters of the fitted model)
  Tp.theta=function(par){
    theta21=par[1]
    theta22=par[2]
    theta23=par[3]
    
    eta1=(theta1[1]*xstart[1,])/ ((theta1[2]*(1+(xstart[2,]/theta1[3]))) + xstart[1,]*(1+(((1-theta1[4])*xstart[2,])/theta1[3]))) #11
    eta2=(theta21*xstart[1,])/ ((theta22*(1+(xstart[2,]/theta23)))+xstart[1,]*(1+(((1-lambda)*xstart[2,])/theta23))) #10
    eta3=wstart*(eta1-eta2)^2
    
    log(sum(eta3))
  }
  
  
  theta2.hat=optim(par=theta2.start,Tp.theta,lower=c(0.001,0.001,0.001),upper=c(Inf,Inf,Inf),method="L-BFGS-B")$par
  theta21.hat=theta2.hat[1]
  theta22.hat=theta2.hat[2]
  theta23.hat=theta2.hat[3]
  THETA.HAT[iter,]=theta2.hat
  
  theta2.start=theta2.hat
  
  #-----------------------------------------------------------------------------------------------------------------------
  #                                                        X.HAT
  #-----------------------------------------------------------------------------------------------------------------------
  
  #delta_0 in T-optimality definition (lack of fit sum of squares when the non competitive model is fixed)
  DELTA=function(TETA1,TETA2){
    ETA1=(TETA1[1]*xstart[1,])/((TETA1[2]*(1+(xstart[2,]/TETA1[3])))+(xstart[1,]*(1+((xstart[2,]*(1-TETA1[4]))/TETA1[3])))) 
    ETA2=(TETA2[1]*xstart[1,])/((TETA2[2]*(1+(xstart[2,]/TETA2[3])))+xstart[1,]*(1+(((1-TETA2[4])*xstart[2,])/TETA2[3]))) 
    ETA3=sum(wstart*((ETA1-ETA2)^2))
    ETA3
  }
  
  DELTA1=DELTA(theta1,c(theta2.hat,lambda))
  #----------------------------------------------------
  
  #second step optimization (maximization with respect to design points)
  psi=function(xest11,xest22){
    xest1=xest11
    xest2=xest22
    eta11=(theta1[1]*xest1)/ ((theta1[2]*(1+(xest2/theta1[3]))) + xest1*(1+(((1-theta1[4])*xest2)/theta1[3])))
    eta22=(theta21.hat*xest1)/ ((theta22.hat*(1+(xest2/theta23.hat))) + xest1*(1+(((1-lambda)*xest2)/theta23.hat)))
    eta33=(eta11-eta22)^2
    eta33/DELTA1 
  }
  
  
  psi.vec=psi(grid[,1],grid[,2])
  num.max[iter]=which.max(psi.vec)
  x.hat=c(grid[num.max[iter],1],grid[num.max[iter],2])
  
  
  xstart=cbind(xstart,x.hat)
  alpha=1/(iter+1)
  wstart=c(((1-alpha)*wstart),alpha)
  design=rbind(xstart,wstart)
  n=ncol(xstart)
  
  
  delta.0=Tp.theta(theta2.hat)
  psi.0=psi(x.hat[1],x.hat[2])
  efficiency[iter]=delta.0/psi.0
  if(efficiency[iter]>=0.99999)break
  
}

#-----------------------------------------------------------------------------------------------
# some information fur further use in T-efficiency

THETA.HAT[1:10,]
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


# If the optimum support point are more that four, according to T-optimal deigns reported in Atkinson (2012)
# while they may differ only a little,
# this means the algorithm may only need higher numbers of iterations.
# here, support points with considerable weights are chosen
mm=which(design.opt[,3]>0.01)
design.opt[mm,]
