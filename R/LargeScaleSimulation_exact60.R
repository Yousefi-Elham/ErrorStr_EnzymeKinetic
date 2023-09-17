## Section 4.3, A simulation study of discriminating designs, subsection 4.3.2 Exact designs, N = 60

rm(list=ls())
library(BB)

th0E2 <- c(6.0645, 3.2799, 3.3153) #estimated values from nls on original obs. for competitive log model
s0E2 <- c(0.9260 , 0.7288,  0.6041) #std.dev values from nls on original obs. for competitive log model
th1E2 <- c(12.0125, 8.5359, 5.6638) #estimated values from nls on original obs. for non competitive log model
s1E2 <- c(2.0553, 1.5721, 0.8879) #std.dev values from nls on original obs. for non competitive log model


# the standard deviation sigma (computed using nls function for combined log model in 120 original data) hat justified in the paper 
sigh <- 0.5128   



#the design space
epsilon <- 0.02
x1 <- c(epsilon,1:30)
x2 <- 0:60
XE2 <- expand.grid(x1,x2)



# number of Monte Carlo simulations and the number of experiments in exact design of size 60
N.Big<-1000
Nexp.Big <- 100


#coefficient of "sigh" (residual standard error) in big design (exact design of size 60)
csig.big <- 4
#-------------------------------------------------------------------------------------------------
generate_yLN_Example2 <- function(xi, model, th, sig) {
  
  # Generates data from one of the two competitive and non competitive models (considering log transformation)
  # in the design points given by the design xi. 
  
  
  # Arguments:
  # xi ... a matrix nx2, where n is the number of observations
  # model ... the model (either 0 or 1) used for the simulations
  # th ... the 3D parameter theta used for the simulations
  # sig ... the standard deviation of errors used for the simulation
  
  # Model 0
  f0 <- function(S, I, th){
    V <- th[1]
    Km <- th[2]
    Kic <- th[3]
    V * S / (Km * (1 + I/Kic) + S)
  }
  
  # Model 1
  f1 <- function(S, I, th){
    V <- th[1]
    Km <- th[2]
    Kin <- th[3]
    V * S / ((Km + S) * (1 + I/Kin))
  }
  
  if (model == 0) y <- log(f0(xi[, 1], xi[, 2], th)) + rnorm(nrow(xi),0, sd =( sig))
  if (model == 1) y <- log(f1(xi[, 1], xi[, 2], th)) + rnorm(nrow(xi),0, sd =( sig))
  
  return(y)
}

#******************************************************************************************************************
decide_Example2 <- function(xi, y, th0i, th1i) 
{
  # Makes a decision about the true model based on the minimum sum of squares
  # Competing models are the competitive and non competitive models
  #
  # Arguments:
  # xi ... the nx2 matrix of design points (with possible replications)
  # y ... the nx1 vector of observations in the design points
  # th0i, th1i ... initial parameters for the optimization (both are 3D vectors)
  
  # The sum of squares for Model 0
  S0 <- function(x, xi, y){
    sum((y - log((x[1] * xi[, 1] / (x[2] * (1 + xi[, 2]/x[3]) + xi[, 1]))))^2) 
  }
  
  # The sum of squares for Model 1
  S1 <- function(x, xi, y){
    sum((y - log((x[1] * xi[, 1] / ((x[2] + xi[, 1]) * (1 + xi[, 2]/x[3])) )))^2) 
  }
  
  res0<-BBoptim(th0i, S0,xi = xi, y = y, gr=NULL,control=list(ftol=0.001, trace=F), method=1, 
                lower=c(0.001, 0.001, 0.001), upper=c(Inf, 60, 30), quiet=T)
  par0<- res0$par
  
  
  res1<-BBoptim(th1i, S1,xi = xi, y = y, gr=NULL,control=list(ftol=0.001, trace=F), method=1,
                lower=c(0.001, 0.001, 0.001), upper=c(Inf, 60, 30), quiet=T)
  par1<- res1$par
  
  
  if (abs(res0$value- res1$value) < 1e-6) {
    result <- sample(c(0, 1), 1)
  } else {
    result <- as.numeric(res0$value > res1$value)
  }
  
  # The result is either 0 or 1 (=the prediction which model is the true one)
  return(list(result=result, par0=par0, par1=par1))
  
}
#-------------------------------------------------------------------------------------------------
simulate_LN_Example2 <- function(n.experiments, xi, model, th.sim, sig.sim)
{
  
  # Simulate n.experiments from the design xi using model 0 or 1 
  #
  # Arguments:
  # n.experiments ... the number of simulated experiments
  # xi ... an n times 2 matrix of design points for the n observations
  # model ... the model to be used for the simulations (0 or 1)
  # th.sim ... the mean value parameter for the simulations
  # sig.sim ... the standard deviation of errors for the simulations
  #
  # Notes:
  # After each virtual experiment make a decision about the correct model.
  
  # The following thetas are (only) used as initial solutions in the optimization
  th0 <- c(6.0645, 3.2799, 3.3153)
  th1 <- c(12.0125, 8.5359, 5.6638)

  dec <- c(0, 0)
  count2=count3=count4=0
  count1=1
  
  #We repeat until we have Nexp.Big=100 admissible experiments
  repeat{
    
    y <- generate_yLN_Example2(xi, model, th.sim, sig.sim)
    dc.c<- decide_Example2(xi, y, th0, th1)
    dc<- dc.c$result
    PAR0=dc.c$par0
    PAR1=dc.c$par1
    
    
    if(length(PAR0[PAR0<=0])==0&length(PAR1[PAR1<=0])==0&
       PAR0[1]>0.001&PAR0[1]<Inf&PAR1[1]>0.001&PAR1[1]<Inf&
       PAR0[2]>0.001&PAR0[2]<60& PAR1[2]>0.001&PAR1[2]<60&
       PAR0[3]>0.001&PAR0[3]<30& PAR1[3]>0.001&PAR1[3]<30){
      dec[dc + 1] <- dec[dc + 1] + 1
      count1<-count1+1
    }else{
      dec[dc + 1] <- dec[dc + 1] + 0
      count4=count4+1
      
      if(length(PAR0[PAR0<=0])!=0 |length(PAR1[PAR1<=0])!=0 ) count2<-count2+1
      
    }
    if(count1>n.experiments)break
  }
  
  percent=dec/n.experiments
  neg<-count2
  worse=count4
  # percent: The vector with two elements representing
  # the proportions of predictions that the true model is Model 0 and Model 1 
  # neg and worse count the non suitable parameter estimates
  return(list(percent=percent,neg=neg,worse=worse))
}

##################################################################################################################
des.A1.Big <- matrix(c(rep(c(0.02, 0),1),rep(c(30.00, 0),8),rep(c(0.02, 60),22),
                       rep(c(30.00, 60),29)), ncol = 2, byrow = TRUE)

des.A2.Big <- matrix(c(rep(c(0.02, 0),10),rep(c(30.00, 0),11),
                       rep(c(0.02, 60),18),rep(c(30.00, 60),21)), ncol = 2, byrow = TRUE)

des.A3.Big <- matrix(c(rep(c(0.02, 0),10),rep(c(30.00, 0),13),
                       rep(c(0.02, 60),17),rep(c(30.00, 60),20)), ncol = 2, byrow = TRUE)

des.A4.Big <- matrix(c(rep(c(0.02, 0),15),rep(c(30.00, 0),15),
                       rep(c(0.02, 60),15),rep(c(30.00, 60),15)), ncol = 2, byrow = TRUE)
#----------------------------------------------------------------------------------------

des.delta1.Big <- matrix(c(rep(c(30.00, 0),11),rep(c(0.02, 60),19),
                           rep(c(30.00, 60),30)), ncol = 2, byrow = TRUE)

des.delta2.Big <- matrix(c(rep(c(0.02, 0),1),rep(c(30.00, 0),10),
                           rep(c(0.02, 60),20),rep(c(30.00, 60),29)), ncol = 2, byrow = TRUE)

#des.delta3.Big <- matrix(c(rep(c(0.02, 0),1),rep(c(30.00, 0),10),
#                           rep(c(0.02, 60),20),rep(c(30.00, 60),29)), ncol = 2, byrow = TRUE)

#des.delta4.Big <- matrix(c(rep(c(0.02, 0),1),rep(c(30.00, 0),10),
#                           rep(c(0.02, 60),20),rep(c(30.00, 60),29)), ncol = 2, byrow = TRUE)
#---------------------------------------------------------------------------------------------------
des.delta5.Big <- matrix(c(rep(c(1, 0),2),rep(c(30.00, 0),11),
                           rep(c(0.02, 60),19),rep(c(30.00, 60),28)), ncol = 2, byrow = TRUE)

#des.delta51.Big <- matrix(c(rep(c(0.02, 0),1),rep(c(30.00, 0),10),
#                            rep(c(0.02, 60),20),rep(c(30.00, 60),29)), ncol = 2, byrow = TRUE)

#des.delta52.Big <- matrix(c(rep(c(0.02, 0),1),rep(c(30.00, 0),10),
#                           rep(c(0.02, 60),20),rep(c(30.00, 60),29)), ncol = 2, byrow = TRUE)

#des.delta2.Big=des.delta3.Big=des.delta4.Big=des.delta51.Big=des.delta52.Big (therefore we run the simulations for one of them)
#----------------------------------------------------------------------------------------------------
des.delta6newa.Big <- matrix(c(rep(c(1, 0),12),rep(c(30.00, 1),18),
                               rep(c(1, 60),12),rep(c(30.00, 60),18)), ncol = 2, byrow = TRUE)

des.delta6newb.Big <- matrix(c(rep(c(2, 0),12),rep(c(30.00, 1),5),rep(c(30.00, 2),11),
                               rep(c(3, 60),14),rep(c(30.00, 60),18)), ncol = 2, byrow = TRUE)

des.delta6newc.Big <- matrix(c(rep(c(0.02, 0),7),rep(c(4, 0),18),rep(c(30, 1),3),rep(c(30, 2),12),
                               rep(c(5, 60),6), rep(c(6, 60),2), rep(c(30, 60),12) ), ncol = 2, byrow = TRUE)

#---------------------------------------------------------------------------------------------------
des.delta6.Big <- matrix(c(rep(c(0.02, 0),7),rep(c(4, 0),18), rep(c(30.00, 1),3),rep(c(30.00, 2),12),
                           rep(c(5.00, 60),8),rep(c(30.00, 60),12)), ncol = 2, byrow = TRUE)

des.delta61.Big <- matrix(c(rep(c(0.02, 0),7),rep(c(4, 0),19), rep(c(30.00, 1),4),rep(c(30.00, 2),12),
                            rep(c(4.00, 60),1),rep(c(5.00, 60),6),rep(c(30.00, 60),11)), ncol = 2, byrow = TRUE)

des.delta62.Big <- matrix(c(rep(c(0.02, 0),5),rep(c(4, 0),15), rep(c(30.00, 3),15),
                            rep(c(4.00, 60),10),rep(c(30.00, 60),15)), ncol = 2, byrow = TRUE)
#---------------------------------------------------------------------------------------------------
#des.delta7.Big=des.delta71.Big=des.delta72.Big=des.delta61.Big (therefore we run the simulations for one of them)

#des.delta7.Big <- matrix(c(rep(c(0.02, 0),7),rep(c(4.00, 0),19),
#                           rep(c(30.00, 1),4), rep(c(30.00, 2),12),rep(c(4.00, 60),1),rep(c(5.00, 60),6),
#                           rep(c(30.00, 60),11)), ncol = 2, byrow = TRUE)

#des.delta71.Big <- matrix(c(rep(c(0.02, 0),7),rep(c(4.00, 0),19),
#                           rep(c(30.00, 1),4),rep(c(30.00, 2),12),rep(c(4.00, 60),1),rep(c(5.00, 60),6),
#                           rep(c(30.00, 60),11)), ncol = 2, byrow = TRUE)

#des.delta72.Big <- matrix(c(rep(c(0.02, 0),7),rep(c(4.00, 0),19),
#                            rep(c(30.00, 1),4),rep(c(30.00, 2),12),rep(c(4.00, 60),1),rep(c(5.00, 60),6),
#                            rep(c(30.00, 60),11)), ncol = 2, byrow = TRUE)


----------------------------------------------------------------------------------------------------
# Monte-Carlo study of the classification performance of des.A1.Big
set.seed(123456789)
res.A1.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.A1.Big=C1.A1.Big <- numeric(N.Big)
Neg.C0.A1.Big=Neg.C1.A1.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.A1.Big, 0, th0E2, csig.big*sigh)
      res.A1.Big[i, 1] <- sim0$percent[1]
      C0.A1.Big[i] <- sim0$worse
      Neg.C0.A1.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.A1.Big, 1, th1E2, csig.big*sigh)
      res.A1.Big[i, 2] <- sim1$percent[2]
      C1.A1.Big[i] <- sim1$worse
      Neg.C1.A1.Big[i] <- sim1$neg
    }
    
  }
}
C0.A1.Big
Neg.C0.A1.Big


sum(C0.A1.Big)
sum(C1.A1.Big)


C1.A1.Big
Neg.C1.A1.Big

sum(Neg.C0.A1.Big)
sum(Neg.C1.A1.Big)
#------------------------------------
# Monte-Carlo study of the classification performance of des.A2.Big
set.seed(123456789)
res.A2.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.A2.Big=C1.A2.Big <- numeric(N.Big)
Neg.C0.A2.Big=Neg.C1.A2.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.A2.Big, 0, th0E2, csig.big*sigh)
      res.A2.Big[i, 1] <- sim0$percent[1]
      C0.A2.Big[i] <- sim0$worse
      Neg.C0.A2.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.A2.Big, 1, th1E2, csig.big*sigh)
      res.A2.Big[i, 2] <- sim1$percent[2]
      C1.A2.Big[i] <- sim1$worse
      Neg.C1.A2.Big[i] <- sim1$neg
    }
    
  }
}
C0.A2.Big
Neg.C0.A2.Big


sum(C0.A2.Big)
sum(C1.A2.Big)


C1.A2.Big
Neg.C1.A2.Big

sum(Neg.C0.A2.Big)
sum(Neg.C1.A2.Big)
#---------------------------------------
# Monte-Carlo study of the classification performance of des.A3.Big
set.seed(123456789)
res.A3.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.A3.Big=C1.A3.Big <- numeric(N.Big)
Neg.C0.A3.Big=Neg.C1.A3.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.A3.Big, 0, th0E2, csig.big*sigh)
      res.A3.Big[i, 1] <- sim0$percent[1]
      C0.A3.Big[i] <- sim0$worse
      Neg.C0.A3.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.A3.Big, 1, th1E2, csig.big*sigh)
      res.A3.Big[i, 2] <- sim1$percent[2]
      C1.A3.Big[i] <- sim1$worse
      Neg.C1.A3.Big[i] <- sim1$neg
    }
    
  }
}
C0.A3.Big
Neg.C0.A3.Big


sum(C0.A3.Big)
sum(C1.A3.Big)


C1.A3.Big
Neg.C1.A3.Big

sum(Neg.C0.A3.Big)
sum(Neg.C1.A3.Big)
#---------------------------------------
# Monte-Carlo study of the classification performance of des.A4.Big
set.seed(123456789)
res.A4.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.A4.Big=C1.A4.Big <- numeric(N.Big)
Neg.C0.A4.Big=Neg.C1.A4.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.A4.Big, 0, th0E2, csig.big*sigh)
      res.A4.Big[i, 1] <- sim0$percent[1]
      C0.A4.Big[i] <- sim0$worse
      Neg.C0.A4.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.A4.Big, 1, th1E2, csig.big*sigh)
      res.A4.Big[i, 2] <- sim1$percent[2]
      C1.A4.Big[i] <- sim1$worse
      Neg.C1.A4.Big[i] <- sim1$neg
    }
    
  }
}
C0.A4.Big
Neg.C0.A4.Big


sum(C0.A4.Big)
sum(C1.A4.Big)


C1.A4.Big
Neg.C1.A4.Big

sum(Neg.C0.A4.Big)
sum(Neg.C1.A4.Big)
#---------------------------------------
# Monte-Carlo study of the classification performance of des.delta1.Big
set.seed(123456789)
res.delta1.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.delta1.Big=C1.delta1.Big <- numeric(N.Big)
Neg.C0.delta1.Big=Neg.C1.delta1.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.delta1.Big, 0, th0E2, csig.big*sigh)
      res.delta1.Big[i, 1] <- sim0$percent[1]
      C0.delta1.Big[i] <- sim0$worse
      Neg.C0.delta1.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.delta1.Big, 1, th1E2, csig.big*sigh)
      res.delta1.Big[i, 2] <- sim1$percent[2]
      C1.delta1.Big[i] <- sim1$worse
      Neg.C1.delta1.Big[i] <- sim1$neg
    }
    
  }
}
C0.delta1.Big
Neg.C0.delta1.Big


sum(C0.delta1.Big)
sum(C1.delta1.Big)


C1.delta1.Big
Neg.C1.delta1.Big

sum(Neg.C0.delta1.Big)
sum(Neg.C1.delta1.Big)
#---------------------------------------
# Monte-Carlo study of the classification performance of des.delta2.Big
set.seed(123456789)
res.delta2.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.delta2.Big=C1.delta2.Big <- numeric(N.Big)
Neg.C0.delta2.Big=Neg.C1.delta2.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.delta2.Big, 0, th0E2, csig.big*sigh)
      res.delta2.Big[i, 1] <- sim0$percent[1]
      C0.delta2.Big[i] <- sim0$worse
      Neg.C0.delta2.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.delta2.Big, 1, th1E2, csig.big*sigh)
      res.delta2.Big[i, 2] <- sim1$percent[2]
      C1.delta2.Big[i] <- sim1$worse
      Neg.C1.delta2.Big[i] <- sim1$neg
    }
    
  }
}
C0.delta2.Big
Neg.C0.delta2.Big


sum(C0.delta2.Big)
sum(C1.delta2.Big)


C1.delta2.Big
Neg.C1.delta2.Big

sum(Neg.C0.delta2.Big)
sum(Neg.C1.delta2.Big)
#---------------------------------------
# Monte-Carlo study of the classification performance of des.delta5.Big
set.seed(123456789)
res.delta5.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.delta5.Big=C1.delta5.Big <- numeric(N.Big)
Neg.C0.delta5.Big=Neg.C1.delta5.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.delta5.Big, 0, th0E2, csig.big*sigh)
      res.delta5.Big[i, 1] <- sim0$percent[1]
      C0.delta5.Big[i] <- sim0$worse
      Neg.C0.delta5.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.delta5.Big, 1, th1E2, csig.big*sigh)
      res.delta5.Big[i, 2] <- sim1$percent[2]
      C1.delta5.Big[i] <- sim1$worse
      Neg.C1.delta5.Big[i] <- sim1$neg
    }
    
  }
}
C0.delta5.Big
Neg.C0.delta5.Big


sum(C0.delta5.Big)
sum(C1.delta5.Big)


C1.delta5.Big
Neg.C1.delta5.Big

sum(Neg.C0.delta5.Big)
sum(Neg.C1.delta5.Big)
#---------------------------------------
# Monte-Carlo study of the classification performance of des.delta6newa.Big
set.seed(123456789)
res.delta6newa.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.delta6newa.Big=C1.delta6newa.Big <- numeric(N.Big)
Neg.C0.delta6newa.Big=Neg.C1.delta6newa.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.delta6newa.Big, 0, th0E2, csig.big*sigh)
      res.delta6newa.Big[i, 1] <- sim0$percent[1]
      C0.delta6newa.Big[i] <- sim0$worse
      Neg.C0.delta6newa.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.delta6newa.Big, 1, th1E2, csig.big*sigh)
      res.delta6newa.Big[i, 2] <- sim1$percent[2]
      C1.delta6newa.Big[i] <- sim1$worse
      Neg.C1.delta6newa.Big[i] <- sim1$neg
    }
    
  }
}
C0.delta6newa.Big
Neg.C0.delta6newa.Big


sum(C0.delta6newa.Big)
sum(C1.delta6newa.Big)


C1.delta6newa.Big
Neg.C1.delta6newa.Big

sum(Neg.C0.delta6newa.Big)
sum(Neg.C1.delta6newa.Big)
#------------------------------------------------------------
# Monte-Carlo study of the classification performance of des.delta6newb.Big
set.seed(123456789)
res.delta6newb.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.delta6newb.Big=C1.delta6newb.Big <- numeric(N.Big)
Neg.C0.delta6newb.Big=Neg.C1.delta6newb.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.delta6newb.Big, 0, th0E2, csig.big*sigh)
      res.delta6newb.Big[i, 1] <- sim0$percent[1]
      C0.delta6newb.Big[i] <- sim0$worse
      Neg.C0.delta6newb.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.delta6newb.Big, 1, th1E2, csig.big*sigh)
      res.delta6newb.Big[i, 2] <- sim1$percent[2]
      C1.delta6newb.Big[i] <- sim1$worse
      Neg.C1.delta6newb.Big[i] <- sim1$neg
    }
    
  }
}
C0.delta6newb.Big
Neg.C0.delta6newb.Big


sum(C0.delta6newb.Big)
sum(C1.delta6newb.Big)


C1.delta6newb.Big
Neg.C1.delta6newb.Big

sum(Neg.C0.delta6newb.Big)
sum(Neg.C1.delta6newb.Big)
#------------------------------------------------------------
# Monte-Carlo study of the classification performance of des.delta6newc.Big
set.seed(123456789)
res.delta6newc.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.delta6newc.Big=C1.delta6newc.Big <- numeric(N.Big)
Neg.C0.delta6newc.Big=Neg.C1.delta6newc.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.delta6newc.Big, 0, th0E2, csig.big*sigh)
      res.delta6newc.Big[i, 1] <- sim0$percent[1]
      C0.delta6newc.Big[i] <- sim0$worse
      Neg.C0.delta6newc.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.delta6newc.Big, 1, th1E2, csig.big*sigh)
      res.delta6newc.Big[i, 2] <- sim1$percent[2]
      C1.delta6newc.Big[i] <- sim1$worse
      Neg.C1.delta6newc.Big[i] <- sim1$neg
    }
    
  }
}
C0.delta6newc.Big
Neg.C0.delta6newc.Big


sum(C0.delta6newc.Big)
sum(C1.delta6newc.Big)


C1.delta6newc.Big
Neg.C1.delta6newc.Big

sum(Neg.C0.delta6newc.Big)
sum(Neg.C1.delta6newc.Big)
#------------------------------------------------------------
# Monte-Carlo study of the classification performance of des.delta6.Big
set.seed(123456789)
res.delta6.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.delta6.Big=C1.delta6.Big <- numeric(N.Big)
Neg.C0.delta6.Big=Neg.C1.delta6.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.delta6.Big, 0, th0E2, csig.big*sigh)
      res.delta6.Big[i, 1] <- sim0$percent[1]
      C0.delta6.Big[i] <- sim0$worse
      Neg.C0.delta6.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.delta6.Big, 1, th1E2, csig.big*sigh)
      res.delta6.Big[i, 2] <- sim1$percent[2]
      C1.delta6.Big[i] <- sim1$worse
      Neg.C1.delta6.Big[i] <- sim1$neg
    }
    
  }
}
C0.delta6.Big
Neg.C0.delta6.Big


sum(C0.delta6.Big)
sum(C1.delta6.Big)


C1.delta6.Big
Neg.C1.delta6.Big

sum(Neg.C0.delta6.Big)
sum(Neg.C1.delta6.Big)
#------------------------------------------------------------
# Monte-Carlo study of the classification performance of des.delta61.Big
set.seed(123456789)
res.delta61.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.delta61.Big=C1.delta61.Big <- numeric(N.Big)
Neg.C0.delta61.Big=Neg.C1.delta61.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.delta61.Big, 0, th0E2, csig.big*sigh)
      res.delta61.Big[i, 1] <- sim0$percent[1]
      C0.delta61.Big[i] <- sim0$worse
      Neg.C0.delta61.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.delta61.Big, 1, th1E2, csig.big*sigh)
      res.delta61.Big[i, 2] <- sim1$percent[2]
      C1.delta61.Big[i] <- sim1$worse
      Neg.C1.delta61.Big[i] <- sim1$neg
    }
    
  }
}
C0.delta61.Big
Neg.C0.delta61.Big


sum(C0.delta61.Big)
sum(C1.delta61.Big)


C1.delta61.Big
Neg.C1.delta61.Big

sum(Neg.C0.delta61.Big)
sum(Neg.C1.delta61.Big)
#---------------------------------------
# Monte-Carlo study of the classification performance of des.delta62.Big
set.seed(123456789)
res.delta62.Big <- matrix(0, nrow = N.Big, ncol = 2)
C0.delta62.Big=C1.delta62.Big <- numeric(N.Big)
Neg.C0.delta62.Big=Neg.C1.delta62.Big <- numeric(N.Big)
for (i in 1:N.Big) {
  
  for (im in 1:2) {
    
    par(mfrow=c(1,2))
    plot(i,i,xlim=c(i-1,i+1))
    text(x=i+1,y=i,i,col="red")
    
    plot(im,im,xlim=c(im-1,im+1))
    text(x=im+1,y=im,im,col="darkgreen")
    
    
    if (im == 1){
      
      sim0<- simulate_LN_Example2(Nexp.Big, des.delta62.Big, 0, th0E2, csig.big*sigh)
      res.delta62.Big[i, 1] <- sim0$percent[1]
      C0.delta62.Big[i] <- sim0$worse
      Neg.C0.delta62.Big[i] <- sim0$neg
      
    }
    
    if (im == 2){
      
      sim1<-simulate_LN_Example2(Nexp.Big, des.delta62.Big, 1, th1E2, csig.big*sigh)
      res.delta62.Big[i, 2] <- sim1$percent[2]
      C1.delta62.Big[i] <- sim1$worse
      Neg.C1.delta62.Big[i] <- sim1$neg
    }
    
  }
}
C0.delta62.Big
Neg.C0.delta62.Big


sum(C0.delta62.Big)
sum(C1.delta62.Big)


C1.delta62.Big
Neg.C1.delta62.Big

sum(Neg.C0.delta62.Big)
sum(Neg.C1.delta62.Big)


#---------------------------------------------------------------------------------------------
## Note that the following parts regarding the information on unsuitable parameter estimates are not reported in the paper
p.Big<- c( sum(C0.A1.Big),sum(C1.A1.Big),sum(C0.A2.Big),sum(C1.A2.Big),sum(C0.A3.Big),sum(C1.A3.Big),sum(C0.A4.Big),sum(C1.A4.Big),
           sum(C0.delta1.Big),sum(C1.delta1.Big),sum(C0.delta2.Big), sum(C1.delta2.Big), sum(C0.delta5.Big), sum(C1.delta5.Big),
           sum(C0.delta6newa.Big), sum(C1.delta6newa.Big),sum(C0.delta6newb.Big), sum(C1.delta6newb.Big),sum(C0.delta6newc.Big),
           sum(C1.delta6newc.Big),sum(C0.delta6.Big), sum(C1.delta6.Big),sum(C0.delta61.Big), sum(C1.delta61.Big),
           sum(C0.delta62.Big), sum(C1.delta62.Big) )


p.Big

pmean<-c(mean(c(sum(C0.A1.Big),sum(C1.A1.Big))),mean(c(sum(C0.A2.Big),sum(C1.A2.Big))),
         mean(c(sum(C0.A3.Big),sum(C1.A3.Big))),mean(c(sum(C0.A4.Big),sum(C1.A4.Big))),
         mean(c(sum(C0.delta1.Big),sum(C1.delta1.Big))),mean(c(sum(C0.delta2.Big),sum(C1.delta2.Big))),
         mean(c(sum(C0.delta5.Big),sum(C1.delta5.Big))),mean(c(sum(C0.delta6newa.Big),sum(C1.delta6newa.Big))),
         mean(c(sum(C0.delta6newb.Big),sum(C1.delta6newb.Big))),mean(c(sum(C0.delta6newc.Big),sum(C1.delta6newc.Big))),
         mean(c(sum(C0.delta6.Big),sum(C1.delta6.Big))),mean(c(sum(C0.delta61.Big),sum(C1.delta61.Big))),
         mean(c(sum(C0.delta62.Big),sum(C1.delta62.Big))) )


pmean

NN.Big<- c( sum(Neg.C0.A1.Big),sum(Neg.C1.A1.Big),sum(Neg.C0.A2.Big),sum(Neg.C1.A2.Big),sum(Neg.C0.A3.Big),sum(Neg.C1.A3.Big),sum(Neg.C0.A4.Big),sum(Neg.C1.A4.Big),
            sum(Neg.C0.delta1.Big),sum(Neg.C1.delta1.Big),sum(Neg.C0.delta2.Big), sum(Neg.C1.delta2.Big), sum(Neg.C0.delta5.Big), sum(Neg.C1.delta5.Big),
            sum(Neg.C0.delta6newa.Big),sum(Neg.C1.delta6newa.Big),sum(Neg.C0.delta6newb.Big), sum(Neg.C1.delta6newb.Big),sum(Neg.C0.delta6newc.Big),
            sum(Neg.C1.delta6newc.Big),sum(Neg.C0.delta6.Big), sum(Neg.C1.delta6.Big),sum(Neg.C0.delta61.Big), sum(Neg.C1.delta61.Big),
            sum(Neg.C0.delta62.Big), sum(Neg.C1.delta62.Big) )
NN.Big
#############################################################################################
#Boxplots for the total correct classification rates of all designs (sec. 4.3)
par(mfrow=c(1,1))
RES.BIG=cbind(res.A1.Big, res.A2.Big, res.A3.Big,res.A4.Big,
              res.delta1.Big, res.delta2.Big, res.delta5.Big,
              res.delta6newa.Big,res.delta6newb.Big,res.delta6newc.Big,res.delta6.Big,
              res.delta61.Big,res.delta62.Big)

NAMES=c(expression(A[1]),expression(A[2]),expression(A[3]),expression(A[4]),
        expression(delta[1]),expression(delta[2]),expression(delta[5*a]),
        expression(delta[6*a]),expression(delta[6*b]),expression(delta[6*c]),
        expression(delta[10*a]),expression(delta[10*b]),expression(delta[10*c]) )


boxplot(RES.BIG,col = rep(c("white", "gray"),13),xaxt="n")
axis(1, at = seq(1.5 , 25.5 , 2), labels = NAMES , tick=FALSE,cex=0.2)
for(i in seq(0.5 , 27.5 , 2)){ 
  abline(v=i,lty=3, col="grey",lwd=1)
}
xtick=seq(1 , 26 , 1)
axis(1,at=xtick,labels=F)
#############################################################################################
#Boxplots for averages of classification rates of all designs (sec. 4.3)

newres.A1.Big=apply(res.A1.Big,1,mean)
newres.A2.Big=apply(res.A2.Big,1,mean)
newres.A3.Big=apply(res.A3.Big,1,mean)
newres.A4.Big=apply(res.A4.Big,1,mean)
newres.delta1.Big=apply(res.delta1.Big,1,mean)
newres.delta2.Big=apply(res.delta2.Big,1,mean)
newres.delta5.Big=apply(res.delta5.Big,1,mean)
newres.delta6newa.Big=apply(res.delta6newa.Big,1,mean)
newres.delta6newb.Big=apply(res.delta6newb.Big,1,mean)
newres.delta6newc.Big=apply(res.delta6newc.Big,1,mean)
newres.delta6.Big=apply(res.delta6.Big,1,mean)
newres.delta61.Big=apply(res.delta61.Big,1,mean)
newres.delta62.Big=apply(res.delta62.Big,1,mean)

par(mfrow=c(1,1))
NEWRES.BIG=cbind(newres.A1.Big,newres.A2.Big,newres.A3.Big,newres.A4.Big,
                 newres.delta1.Big,newres.delta2.Big,newres.delta5.Big,newres.delta6newa.Big,newres.delta6newb.Big,
                 newres.delta6newc.Big,newres.delta6.Big,newres.delta61.Big,newres.delta62.Big)

NAMES=c(expression(A[1]),expression(A[2]),expression(A[3]),expression(A[4]),
        expression(delta[1]),expression(delta[2]),expression(delta[5*a]),
        expression(delta[6*a]),expression(delta[6*b]),expression(delta[6*c]),
        expression(delta[10*a]),expression(delta[10*b]),expression(delta[10*c]))

boxplot(NEWRES.BIG,col = rep("grey",13),xaxt="n")
axis(1, at = seq(1 , 13.5 , 1), labels = NAMES , tick=FALSE,cex=0.2)
xtick=seq(1 , 13 , 1)
axis(1,at=xtick,labels=F)