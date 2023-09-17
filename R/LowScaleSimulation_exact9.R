## Section 4.3, A simulation study of discriminating designs, subsection 4.3.1 Exact designs, N = 9

rm(list=ls())
library(BB)

th0E2 <- c(6.0645, 3.2799, 3.3153) #estimated values from nls on original obs. for competitive log model
s0E2 <- c(0.9260 , 0.7288,  0.6041) #std.dev values from nls on original obs. for competitive log model
th1E2 <- c(12.0125, 8.5359, 5.6638) #estimated values from nls on original obs. for non competitive log model
s1E2 <- c(2.0553, 1.5721, 0.8879) #std.dev values from nls on original obs. for non competitive log model


# the standard deviation sigma (computed using nls function for combined log model in 120 original data) hat justified in the paper 
sigh <- 0.5128


# the values used for perturbing the parameters in first simulation part 
# previously various values was used but then we observed that perturbation does not effect the general performance of designs
# therefore we did NO perturbation (c.vec= 0)
c.vec <- 0

#the design space
epsilon <- 0.02
x1 <- c(epsilon,1:30)
x2 <- 0:60
XE2 <- expand.grid(x1,x2)


# number of Monte Carlo simulations and the number of experiments for exact designs of size 6-9
N1<-100
Nexp <- 100


#coefficient of sigh (residual standard error) in exact design of size 6-9
csig <- 1
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

#-------------------------------------------------------------------------------------------------
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
  
  #We repeat until we have Nexp=100 admissible experiments
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

###################################################################################################

#N=9 number of exact design


des.A1 <- matrix(c(rep(c(30.00, 0),1),rep(c(0.02, 60),3),
                   rep(c(30.00, 60),5)), ncol = 2, byrow = TRUE)


des.A2 <- matrix(c(rep(c(0.02, 0),1),rep(c(30.00, 0),2),
                   rep(c(0.02, 60),3),rep(c(30.00, 60),3)), ncol = 2, byrow = TRUE)


#des.A3=des.A2 <- matrix(c(rep(c(0.02, 0),1),rep(c(30.00, 0),2),
#                   rep(c(0.02, 60),3),rep(c(30.00, 60),3)), ncol = 2, byrow = TRUE)


des.A4 <- matrix(c(rep(c(0.02, 0),2),rep(c(30.00, 0),3),
                   rep(c(0.02, 60),2),rep(c(30.00, 60),2)), ncol = 2, byrow = TRUE)

#---------------------------------------------------------------------------------------------------
des.delta1  <- matrix(c(30.00,    0,
                        30.00,    0,
                        0.02,   60,
                        0.02,   60,
                        0.02,   60,
                        30.00,   60,
                        30.00,   60,
                        30.00,   60,
                        30.00,   60) , ncol = 2, byrow = TRUE)


des.delta4  <- matrix(c(0.02,    0,
                        30.00,    0,
                        30.00,    0,
                        0.02,   60,
                        0.02,   60,
                        30.00,   60,
                        30.00,   60,
                        30.00,   60,
                        30.00,   60) , ncol = 2, byrow = TRUE)


des.delta5  <- matrix(c(1.00,    0,
                        30.00,    0,
                        30.00,    0,
                        0.02,   60,
                        0.02,   60,
                        30.00,   60,
                        30.00,   60,
                        30.00,   60,
                        30.00,   60) , ncol = 2, byrow = TRUE)


des.delta52  <- matrix(c(1.00,    0,
                         30.00,    0,
                         0.02,   60,
                         0.02,   60,
                         0.02,   60,
                         30.00,   60,
                         30.00,   60,
                         30.00,   60,
                         30.00,   60) , ncol = 2, byrow = TRUE)


des.delta6  <- matrix(c(0.02,    0,
                        3.00,    0,
                        4.00,    0,
                        4.00,    0,
                        30.00,    2,
                        30.00,    2,
                        5.00,   60,
                        30.00,   60,
                        30.00,   60) , ncol = 2, byrow = TRUE)


des.delta61  <- matrix(c(0.02,    0,
                         4.00,    0,
                         4.00,    0,
                         4.00,    0,
                         30.00,    1,
                         30.00,    3,
                         5.00,   60,
                         30.00,   60,
                         30.00,   60) , ncol = 2, byrow = TRUE)


des.delta62  <- matrix(c(0.02,    0,
                         4.00,    0,
                         4.00,    0,
                         30.00,    3,
                         30.00,    3,
                         4.00,   60,
                         4.00,   60,
                         30.00,   60,
                         30.00,   60) , ncol = 2, byrow = TRUE)



####################################################################################################
# Note that we don't report the matrices C0, C1, Neg.C0 and Neg.C1 for all designs, in the paper.


# Classification performance of des.A1
set.seed(123456789)
res.A1 <- matrix(0, nrow = length(c.vec), ncol = 2)
C0.A1 <- matrix(0, nrow = N1, ncol = length(c.vec))
C1.A1 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C0.A1 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C1.A1 <- matrix(0, nrow = N1, ncol = length(c.vec))

for (ic in 1:length(c.vec)) {
  for (im in 1:2) {
    res.A1[ic, im] <- 0
    for (i in 1:N1) {
      
      
      par(mfrow=c(1,3))
      
      plot(ic,ic,xlim=c(ic-1,ic+1))
      text(x=ic+1,y=ic,ic,col="blue")
      
      plot(im,im,xlim=c(im-1,im+1))
      text(x=im+1,y=im,im,col="darkgreen")
      
      plot(i,i,xlim=c(i-1,i+1))
      text(x=i+1,y=i,i,col="red")
      
      
      if (im == 1) {
        th0.pert <- runif(3, min = th0E2 - c.vec[ic]*s0E2,
                          max = th0E2 + c.vec[ic]*s0E2)
        
        sim0<-simulate_LN_Example2(Nexp, des.A1, 0, th0.pert, csig*sigh)
        
        res.A1[ic, im] <- res.A1[ic, im] + sim0$percent[1]
        C0.A1[i, ic] <-sim0$worse
        Neg.C0.A1[i, ic] <- sim0$neg
        
      }
      if (im == 2) {
        th1.pert <- runif(3, min = th1E2 - c.vec[ic]*s1E2,
                          max = th1E2 + c.vec[ic]*s1E2)
        
        sim1<-simulate_LN_Example2(Nexp, des.A1, 1, th1.pert, csig*sigh)
        
        res.A1[ic, im] <- res.A1[ic, im] + sim1$percent[2]
        C1.A1[i, ic] <-sim1$worse
        Neg.C1.A1[i, ic] <- sim1$neg
      }
    }
  }
}
res.A1 <- res.A1 / N1

#---------------------------------------------------------------------------------------------------
# Classification performance of des.A2
set.seed(123456789)
res.A2 <- matrix(0, nrow = length(c.vec), ncol = 2)
C0.A2 <- matrix(0, nrow = N1, ncol = length(c.vec))
C1.A2 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C0.A2 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C1.A2 <- matrix(0, nrow = N1, ncol = length(c.vec))

for (ic in 1:length(c.vec)) {
  for (im in 1:2) {
    res.A2[ic, im] <- 0
    for (i in 1:N1) {
      
      
      par(mfrow=c(1,3))
      
      plot(ic,ic,xlim=c(ic-1,ic+1))
      text(x=ic+1,y=ic,ic,col="blue")
      
      plot(im,im,xlim=c(im-1,im+1))
      text(x=im+1,y=im,im,col="darkgreen")
      
      plot(i,i,xlim=c(i-1,i+1))
      text(x=i+1,y=i,i,col="red")
      
      
      if (im == 1) {
        th0.pert <- runif(3, min = th0E2 - c.vec[ic]*s0E2,
                          max = th0E2 + c.vec[ic]*s0E2)
        
        sim0<-simulate_LN_Example2(Nexp, des.A2, 0, th0.pert, csig*sigh)
        
        res.A2[ic, im] <- res.A2[ic, im] + sim0$percent[1]
        C0.A2[i, ic] <-sim0$worse
        Neg.C0.A2[i, ic] <- sim0$neg
        
      }
      if (im == 2) {
        th1.pert <- runif(3, min = th1E2 - c.vec[ic]*s1E2,
                          max = th1E2 + c.vec[ic]*s1E2)
        
        sim1<-simulate_LN_Example2(Nexp, des.A2, 1, th1.pert, csig*sigh)
        
        res.A2[ic, im] <- res.A2[ic, im] + sim1$percent[2]
        C1.A2[i, ic] <-sim1$worse
        Neg.C1.A2[i, ic] <- sim1$neg
      }
    }
  }
}
res.A2 <- res.A2 / N1

#---------------------------------------------------------------------------------------------------
# Classification performance of des.A4 
set.seed(123456789)
res.A4 <- matrix(0, nrow = length(c.vec), ncol = 2)
C0.A4 <- matrix(0, nrow = N1, ncol = length(c.vec))
C1.A4 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C0.A4 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C1.A4 <- matrix(0, nrow = N1, ncol = length(c.vec))

for (ic in 1:length(c.vec)) {
  for (im in 1:2) {
    res.A4[ic, im] <- 0
    for (i in 1:N1) {
      
      
      par(mfrow=c(1,3))
      
      plot(ic,ic,xlim=c(ic-1,ic+1))
      text(x=ic+1,y=ic,ic,col="blue")
      
      plot(im,im,xlim=c(im-1,im+1))
      text(x=im+1,y=im,im,col="darkgreen")
      
      plot(i,i,xlim=c(i-1,i+1))
      text(x=i+1,y=i,i,col="red")
      
      
      if (im == 1) {
        th0.pert <- runif(3, min = th0E2 - c.vec[ic]*s0E2,
                          max = th0E2 + c.vec[ic]*s0E2)
        
        sim0<-simulate_LN_Example2(Nexp, des.A4, 0, th0.pert, csig*sigh)
        
        res.A4[ic, im] <- res.A4[ic, im] + sim0$percent[1]
        C0.A4[i, ic] <-sim0$worse
        Neg.C0.A4[i, ic] <- sim0$neg
        
      }
      if (im == 2) {
        th1.pert <- runif(3, min = th1E2 - c.vec[ic]*s1E2,
                          max = th1E2 + c.vec[ic]*s1E2)
        
        sim1<-simulate_LN_Example2(Nexp, des.A4, 1, th1.pert, csig*sigh)
        
        res.A4[ic, im] <- res.A4[ic, im] + sim1$percent[2]
        C1.A4[i, ic] <-sim1$worse
        Neg.C1.A4[i, ic] <- sim1$neg
      }
    }
  }
}
res.A4 <- res.A4 / N1

#---------------------------------------------------------------------------------------------------
# Classification performance of des.delta1
set.seed(123456789)
res.delta1 <- matrix(0, nrow = length(c.vec), ncol = 2)
C0.delta1 <- matrix(0, nrow = N1, ncol = length(c.vec))
C1.delta1 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C0.delta1 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C1.delta1 <- matrix(0, nrow = N1, ncol = length(c.vec))

for (ic in 1:length(c.vec)) {
  for (im in 1:2) {
    res.delta1[ic, im] <- 0
    for (i in 1:N1) {
      
      
      par(mfrow=c(1,3))
      
      plot(ic,ic,xlim=c(ic-1,ic+1))
      text(x=ic+1,y=ic,ic,col="blue")
      
      plot(im,im,xlim=c(im-1,im+1))
      text(x=im+1,y=im,im,col="darkgreen")
      
      plot(i,i,xlim=c(i-1,i+1))
      text(x=i+1,y=i,i,col="red")
      
      
      if (im == 1) {
        th0.pert <- runif(3, min = th0E2 - c.vec[ic]*s0E2,
                          max = th0E2 + c.vec[ic]*s0E2)
        
        sim0<-simulate_LN_Example2(Nexp, des.delta1, 0, th0.pert, csig*sigh)
        
        res.delta1[ic, im] <- res.delta1[ic, im] + sim0$percent[1]
        C0.delta1[i, ic] <-sim0$worse
        Neg.C0.delta1[i, ic] <- sim0$neg
        
      }
      if (im == 2) {
        th1.pert <- runif(3, min = th1E2 - c.vec[ic]*s1E2,
                          max = th1E2 + c.vec[ic]*s1E2)
        
        sim1<-simulate_LN_Example2(Nexp, des.delta1, 1, th1.pert, csig*sigh)
        
        res.delta1[ic, im] <- res.delta1[ic, im] + sim1$percent[2]
        C1.delta1[i, ic] <-sim1$worse
        Neg.C1.delta1[i, ic] <- sim1$neg
      }
    }
  }
}
res.delta1 <- res.delta1 / N1

#---------------------------------------------------------------------------------------------------
# Classification performance of des.delta4
set.seed(123456789)
res.delta4 <- matrix(0, nrow = length(c.vec), ncol = 2)
C0.delta4 <- matrix(0, nrow = N1, ncol = length(c.vec))
C1.delta4 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C0.delta4 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C1.delta4 <- matrix(0, nrow = N1, ncol = length(c.vec))

for (ic in 1:length(c.vec)) {
  for (im in 1:2) {
    res.delta4[ic, im] <- 0
    for (i in 1:N1) {
      
      
      par(mfrow=c(1,3))
      
      plot(ic,ic,xlim=c(ic-1,ic+1))
      text(x=ic+1,y=ic,ic,col="blue")
      
      plot(im,im,xlim=c(im-1,im+1))
      text(x=im+1,y=im,im,col="darkgreen")
      
      plot(i,i,xlim=c(i-1,i+1))
      text(x=i+1,y=i,i,col="red")
      
      
      if (im == 1) {
        th0.pert <- runif(3, min = th0E2 - c.vec[ic]*s0E2,
                          max = th0E2 + c.vec[ic]*s0E2)
        
        sim0<-simulate_LN_Example2(Nexp, des.delta4, 0, th0.pert, csig*sigh)
        
        res.delta4[ic, im] <- res.delta4[ic, im] + sim0$percent[1]
        C0.delta4[i, ic] <-sim0$worse
        Neg.C0.delta4[i, ic] <- sim0$neg
        
      }
      if (im == 2) {
        th1.pert <- runif(3, min = th1E2 - c.vec[ic]*s1E2,
                          max = th1E2 + c.vec[ic]*s1E2)
        
        sim1<-simulate_LN_Example2(Nexp, des.delta4, 1, th1.pert, csig*sigh)
        
        res.delta4[ic, im] <- res.delta4[ic, im] + sim1$percent[2]
        C1.delta4[i, ic] <-sim1$worse
        Neg.C1.delta4[i, ic] <- sim1$neg
      }
    }
  }
}
res.delta4 <- res.delta4 / N1

#---------------------------------------------------------------------------------------------------
# Classification performance of des.delta5
set.seed(123456789)
res.delta5 <- matrix(0, nrow = length(c.vec), ncol = 2)
C0.delta5 <- matrix(0, nrow = N1, ncol = length(c.vec))
C1.delta5 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C0.delta5 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C1.delta5 <- matrix(0, nrow = N1, ncol = length(c.vec))

for (ic in 1:length(c.vec)) {
  for (im in 1:2) {
    res.delta5[ic, im] <- 0
    for (i in 1:N1) {
      
      
      par(mfrow=c(1,3))
      
      plot(ic,ic,xlim=c(ic-1,ic+1))
      text(x=ic+1,y=ic,ic,col="blue")
      
      plot(im,im,xlim=c(im-1,im+1))
      text(x=im+1,y=im,im,col="darkgreen")
      
      plot(i,i,xlim=c(i-1,i+1))
      text(x=i+1,y=i,i,col="red")
      
      
      if (im == 1) {
        th0.pert <- runif(3, min = th0E2 - c.vec[ic]*s0E2,
                          max = th0E2 + c.vec[ic]*s0E2)
        
        sim0<-simulate_LN_Example2(Nexp, des.delta5, 0, th0.pert, csig*sigh)
        
        res.delta5[ic, im] <- res.delta5[ic, im] + sim0$percent[1]
        C0.delta5[i, ic] <-sim0$worse
        Neg.C0.delta5[i, ic] <- sim0$neg
        
      }
      if (im == 2) {
        th1.pert <- runif(3, min = th1E2 - c.vec[ic]*s1E2,
                          max = th1E2 + c.vec[ic]*s1E2)
        
        sim1<-simulate_LN_Example2(Nexp, des.delta5, 1, th1.pert, csig*sigh)
        
        res.delta5[ic, im] <- res.delta5[ic, im] + sim1$percent[2]
        C1.delta5[i, ic] <-sim1$worse
        Neg.C1.delta5[i, ic] <- sim1$neg
      }
    }
  }
}
res.delta5 <- res.delta5 / N1

#---------------------------------------------------------------------------------------------------
# Classification performance of des.delta52
set.seed(123456789)
res.delta52 <- matrix(0, nrow = length(c.vec), ncol = 2)
C0.delta52 <- matrix(0, nrow = N1, ncol = length(c.vec))
C1.delta52 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C0.delta52 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C1.delta52 <- matrix(0, nrow = N1, ncol = length(c.vec))

for (ic in 1:length(c.vec)) {
  for (im in 1:2) {
    res.delta52[ic, im] <- 0
    for (i in 1:N1) {
      
      
      par(mfrow=c(1,3))
      
      plot(ic,ic,xlim=c(ic-1,ic+1))
      text(x=ic+1,y=ic,ic,col="blue")
      
      plot(im,im,xlim=c(im-1,im+1))
      text(x=im+1,y=im,im,col="darkgreen")
      
      plot(i,i,xlim=c(i-1,i+1))
      text(x=i+1,y=i,i,col="red")
      
      
      if (im == 1) {
        th0.pert <- runif(3, min = th0E2 - c.vec[ic]*s0E2,
                          max = th0E2 + c.vec[ic]*s0E2)
        
        sim0<-simulate_LN_Example2(Nexp, des.delta52, 0, th0.pert, csig*sigh)
        
        res.delta52[ic, im] <- res.delta52[ic, im] + sim0$percent[1]
        C0.delta52[i, ic] <-sim0$worse
        Neg.C0.delta52[i, ic] <- sim0$neg
        
      }
      if (im == 2) {
        th1.pert <- runif(3, min = th1E2 - c.vec[ic]*s1E2,
                          max = th1E2 + c.vec[ic]*s1E2)
        
        sim1<-simulate_LN_Example2(Nexp, des.delta52, 1, th1.pert, csig*sigh)
        
        res.delta52[ic, im] <- res.delta52[ic, im] + sim1$percent[2]
        C1.delta52[i, ic] <-sim1$worse
        Neg.C1.delta52[i, ic] <- sim1$neg
      }
    }
  }
}
res.delta52 <- res.delta52 / N1

#---------------------------------------------------------------------------------------------------
# Classification performance of des.delta6
set.seed(123456789)
res.delta6 <- matrix(0, nrow = length(c.vec), ncol = 2)
C0.delta6 <- matrix(0, nrow = N1, ncol = length(c.vec))
C1.delta6 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C0.delta6 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C1.delta6 <- matrix(0, nrow = N1, ncol = length(c.vec))

for (ic in 1:length(c.vec)) {
  for (im in 1:2) {
    res.delta6[ic, im] <- 0
    for (i in 1:N1) {
      
      
      par(mfrow=c(1,3))
      
      plot(ic,ic,xlim=c(ic-1,ic+1))
      text(x=ic+1,y=ic,ic,col="blue")
      
      plot(im,im,xlim=c(im-1,im+1))
      text(x=im+1,y=im,im,col="darkgreen")
      
      plot(i,i,xlim=c(i-1,i+1))
      text(x=i+1,y=i,i,col="red")
      
      
      if (im == 1) {
        th0.pert <- runif(3, min = th0E2 - c.vec[ic]*s0E2,
                          max = th0E2 + c.vec[ic]*s0E2)
        
        sim0<-simulate_LN_Example2(Nexp, des.delta6, 0, th0.pert, csig*sigh)
        
        res.delta6[ic, im] <- res.delta6[ic, im] + sim0$percent[1]
        C0.delta6[i, ic] <-sim0$worse
        Neg.C0.delta6[i, ic] <- sim0$neg
        
      }
      if (im == 2) {
        th1.pert <- runif(3, min = th1E2 - c.vec[ic]*s1E2,
                          max = th1E2 + c.vec[ic]*s1E2)
        
        sim1<-simulate_LN_Example2(Nexp, des.delta6, 1, th1.pert, csig*sigh)
        
        res.delta6[ic, im] <- res.delta6[ic, im] + sim1$percent[2]
        C1.delta6[i, ic] <-sim1$worse
        Neg.C1.delta6[i, ic] <- sim1$neg
      }
    }
  }
}
res.delta6 <- res.delta6 / N1

#---------------------------------------------------------------------------------------------------
# Classification performance of des.delta61
set.seed(123456789)
res.delta61 <- matrix(0, nrow = length(c.vec), ncol = 2)
C0.delta61 <- matrix(0, nrow = N1, ncol = length(c.vec))
C1.delta61 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C0.delta61 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C1.delta61 <- matrix(0, nrow = N1, ncol = length(c.vec))

for (ic in 1:length(c.vec)) {
  for (im in 1:2) {
    res.delta61[ic, im] <- 0
    for (i in 1:N1) {
      
      
      par(mfrow=c(1,3))
      
      plot(ic,ic,xlim=c(ic-1,ic+1))
      text(x=ic+1,y=ic,ic,col="blue")
      
      plot(im,im,xlim=c(im-1,im+1))
      text(x=im+1,y=im,im,col="darkgreen")
      
      plot(i,i,xlim=c(i-1,i+1))
      text(x=i+1,y=i,i,col="red")
      
      
      if (im == 1) {
        th0.pert <- runif(3, min = th0E2 - c.vec[ic]*s0E2,
                          max = th0E2 + c.vec[ic]*s0E2)
        
        sim0<-simulate_LN_Example2(Nexp, des.delta61, 0, th0.pert, csig*sigh)
        
        res.delta61[ic, im] <- res.delta61[ic, im] + sim0$percent[1]
        C0.delta61[i, ic] <-sim0$worse
        Neg.C0.delta61[i, ic] <- sim0$neg
        
      }
      if (im == 2) {
        th1.pert <- runif(3, min = th1E2 - c.vec[ic]*s1E2,
                          max = th1E2 + c.vec[ic]*s1E2)
        
        sim1<-simulate_LN_Example2(Nexp, des.delta61, 1, th1.pert, csig*sigh)
        
        res.delta61[ic, im] <- res.delta61[ic, im] + sim1$percent[2]
        C1.delta61[i, ic] <-sim1$worse
        Neg.C1.delta61[i, ic] <- sim1$neg
      }
    }
  }
}
res.delta61 <- res.delta61 / N1

#---------------------------------------------------------------------------------------------------
# Classification performance of des.delta62
set.seed(123456789)
res.delta62 <- matrix(0, nrow = length(c.vec), ncol = 2)
C0.delta62 <- matrix(0, nrow = N1, ncol = length(c.vec))
C1.delta62 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C0.delta62 <- matrix(0, nrow = N1, ncol = length(c.vec))
Neg.C1.delta62 <- matrix(0, nrow = N1, ncol = length(c.vec))

for (ic in 1:length(c.vec)) {
  for (im in 1:2) {
    res.delta62[ic, im] <- 0
    for (i in 1:N1) {
      
      
      par(mfrow=c(1,3))
      
      plot(ic,ic,xlim=c(ic-1,ic+1))
      text(x=ic+1,y=ic,ic,col="blue")
      
      plot(im,im,xlim=c(im-1,im+1))
      text(x=im+1,y=im,im,col="darkgreen")
      
      plot(i,i,xlim=c(i-1,i+1))
      text(x=i+1,y=i,i,col="red")
      
      
      if (im == 1) {
        th0.pert <- runif(3, min = th0E2 - c.vec[ic]*s0E2,
                          max = th0E2 + c.vec[ic]*s0E2)
        
        sim0<-simulate_LN_Example2(Nexp, des.delta62, 0, th0.pert, csig*sigh)
        
        res.delta62[ic, im] <- res.delta62[ic, im] + sim0$percent[1]
        C0.delta62[i, ic] <-sim0$worse
        Neg.C0.delta62[i, ic] <- sim0$neg
        
      }
      if (im == 2) {
        th1.pert <- runif(3, min = th1E2 - c.vec[ic]*s1E2,
                          max = th1E2 + c.vec[ic]*s1E2)
        
        sim1<-simulate_LN_Example2(Nexp, des.delta62, 1, th1.pert, csig*sigh)
        
        res.delta62[ic, im] <- res.delta62[ic, im] + sim1$percent[2]
        C1.delta62[i, ic] <-sim1$worse
        Neg.C1.delta62[i, ic] <- sim1$neg
      }
    }
  }
}
res.delta62 <- res.delta62 / N1

####################################################################################################
# some of the designs are equal to each other and therefore their simulations would result into same results

Table3 <- rbind( as.vector(t(res.A1)), as.vector(t(res.A2)),
                 as.vector(t(res.A4)),as.vector(t(res.delta1)),as.vector(t(res.delta4)),
                 as.vector(t(res.delta5)),as.vector(t(res.delta52)),
                 as.vector(t(res.delta6)),as.vector(t(res.delta61)),as.vector(t(res.delta62)) )



print(Table3)

#--------------------------------------------------------------------------------------
# we have finally reported average values of hit rates when both models contribute equally in simulations
# related to Table 7, part N = 9
round(rowMeans(Table3)*100,digits=3)
