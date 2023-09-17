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
theta1nc=c(8.6957,8.0664,12.0566,0) #nu=0, non-competitive model is fixed
theta2.startc=c(7.2976,4.3860,2.5821)
lambdac=1   # lambdac=1, competitive model is fitted

theta1c=c(7.2976,4.3860,2.5821,1) #nu=1, competitive model is fixed
theta2.startnc=c(8.6957,8.0664,12.0566)
lambdanc=0   # lambdanc=0, competitive model is fitted



#Maximum number of algorithm iterations
max.iter=20000
efficiency <- numeric(max.iter)

THETA.HAT=matrix(0,max.iter,6)
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
    theta21c=par[1]
    theta22c=par[2]
    theta23c=par[3]
    
    theta21nc=par[4]
    theta22nc=par[5]
    theta23nc=par[6]
    
    eta1nc=(theta1nc[1]*xstart[1,])/((theta1nc[2]*(1+(xstart[2,]/theta1nc[3]))) + xstart[1,]*(1+(((1-theta1nc[4])*xstart[2,])/theta1nc[3]))) #11
    eta2c=(theta21c*xstart[1,])/((theta22c*(1+(xstart[2,]/theta23c)))+xstart[1,]*(1+(((1-lambdac)*xstart[2,])/theta23c))) #10
    eta3nc.c=wstart*(eta1nc-eta2c)^2
    
    eta1c=(theta1c[1]*xstart[1,])/((theta1c[2]*(1+(xstart[2,]/theta1c[3]))) + xstart[1,]*(1+(((1-theta1c[4])*xstart[2,])/theta1c[3]))) #10
    eta2nc=(theta21nc*xstart[1,])/((theta22nc*(1+(xstart[2,]/theta23nc)))+xstart[1,]*(1+(((1-lambdanc)*xstart[2,])/theta23nc))) #11
    eta3c.nc=wstart*(eta1c-eta2nc)^2
    
    (1-0.5)*log(sum(eta3nc.c))+(0.5)*log(sum(eta3c.nc))
  }
  
  theta2.hat=optim(par=c(theta2.startc,theta2.startnc),Tp.theta,lower=c(0.001,0.001,0.001,0.001,0.001,0.001),upper=c(Inf,Inf,Inf,Inf,Inf,Inf),method="L-BFGS-B")$par
  
  
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
  #delta in T-optimality definition (lack of fit sum of squares as a compound combination of rival models)
  DELTA=function(TETA1,TETA2){
    ETA1=(TETA1[1]*xstart[1,])/((TETA1[2]*(1+(xstart[2,]/TETA1[3])))+(xstart[1,]*(1+((xstart[2,]*(1-TETA1[4]))/TETA1[3])))) 
    ETA2=(TETA2[1]*xstart[1,])/((TETA2[2]*(1+(xstart[2,]/TETA2[3])))+xstart[1,]*(1+(((1-TETA2[4])*xstart[2,])/TETA2[3]))) 
    ETA3=sum(wstart*((ETA1-ETA2)^2))
    ETA3
  }
  
  DELTA1=DELTA(theta1nc,c(theta2.hat[1:3],lambdac))
  DELTA2=DELTA(theta1c,c(theta2.hat[4:6],lambdanc))
  
  
  #second step optimization (maximization with respect to design points)
  psi=function(xest11,xest22){
    xest1=xest11
    xest2=xest22
    eta11nc=(theta1nc[1]*xest1)/ ((theta1nc[2]*(1+(xest2/theta1nc[3]))) + xest1*(1+(((1-theta1nc[4])*xest2)/theta1nc[3]))) #11
    eta22c=(theta21.hatc*xest1)/ ((theta22.hatc*(1+(xest2/theta23.hatc))) + xest1*(1+(((1-lambdac)*xest2)/theta23.hatc))) #10
    eta33nc.c=(eta11nc-eta22c)^2
    
    eta11c=(theta1c[1]*xest1)/ ((theta1c[2]*(1+(xest2/theta1c[3]))) + xest1*(1+(((1-theta1c[4])*xest2)/theta1c[3]))) #10
    eta22nc=(theta21.hatnc*xest1)/ ((theta22.hatnc*(1+(xest2/theta23.hatnc))) + xest1*(1+(((1-lambdanc)*xest2)/theta23.hatnc))) #11
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
  
  delta.0=Tp.theta(theta2.hat)
  psi.0=psi(x.hat[1],x.hat[2])
  efficiency[iter]=delta.0/psi.0
  if(efficiency[iter]>=0.99999) break
  
}
#-----------------------------------------------------------------------------------------------
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

