
# Likelihood ration test for models selection
# The following part in NOT presented in the final version of the manuscript, but exists in the arXiv version.
library(BB)

# initial parameter values for competitive (th0E2) and non competitive (th1E2) models (standard case)
th0E2 <- c(7.298, 4.386, 2.582)
th1E2 <- c(8.696, 8.066, 12.057)

#the standard deviation sigma (computed using 120 original data) hat justified in the paper for the rival models
sigh0<-0.1553
sigh1<-0.2272

# simulation number
N <- 10000

basia = matrix(c(6.640,	30.00,	0.0,
                 6.040,	20.00,	0.0,
                 4.755,	10.00,	0.0,
                 4.105,	7.00,	0.0,
                 3.090,	3.00,	0.0,
                 2.505,	2.00,	0.0,
                 1.455,	1.00,	0.0,
                 1.265,	0.70,	0.0,
                 0.646,	0.30,	0.0,
                 0.525,	0.20,	0.0,
                 0.367,	0.10,	0.0,
                 0.206,	0.07,	0.0,
                 0.088,	0.03,	0.0,
                 0.069,	0.02,	0.0,
                 0.009,	0.01,	0.0,
                 6.085,	30.00,	0.9,
                 5.495,	20.00,	0.9,
                 3.955,	10.00,	0.9,
                 3.540,	7.00,	0.9,
                 2.380,	3.00,	0.9,
                 1.905,	2.00,	0.9,
                 1.185,	1.00,	0.9,
                 0.812,	0.70,	0.9,
                 0.413,	0.30,	0.9,
                 0.341,	0.20,	0.9,
                 0.169,	0.10,	0.9,
                 0.098,	0.07,	0.9,
                 0.062,	0.03,	0.9,
                 0.045,	0.02,	0.9,
                 0.006,	0.01,	0.9,
                 6.030,	30.00,	2.8,
                 5.380,	20.00,	2.8,
                 3.930,	10.00,	2.8,
                 3.320,	7.00,	2.8,
                 2.090,	3.00,	2.8,
                 1.550,	2.00,	2.8,
                 0.865,	1.00,	2.8,
                 0.692,	0.70,	2.8,
                 0.350,	0.30,	2.8,
                 0.245,	0.20,	2.8,
                 0.135,	0.10,	2.8,
                 0.104,	0.07,	2.8,
                 0.043,	0.03,	2.8,
                 0.026,	0.02,	2.8,
                 0.003,	0.01,	2.8,
                 5.570,	30.00,	3.7,
                 4.935,	20.00,	3.7,
                 3.210,	10.00,	3.7,
                 2.695,	7.00,	3.7,
                 2.065,	3.00,	3.7,
                 1.210,	2.00,	3.7,
                 0.659,	1.00,	3.7,
                 0.512,	0.70,	3.7,
                 0.238,	0.30,	3.7,
                 0.188,	0.20,	3.7,
                 0.087,	0.10,	3.7,
                 0.090,	0.07,	3.7,
                 0.026,	0.03,	3.7,
                 0.020,	0.02,	3.7,
                 0.004,	0.01,	3.7,
                 4.480,	30.00,	7.5,
                 4.020,	20.00,	7.5,
                 2.605,	10.00,	7.5,
                 1.935,	7.00,	7.5,
                 1.094,	3.00,	7.5,
                 0.836,	2.00,	7.5,
                 0.467,	1.00,	7.5,
                 0.355,	0.70,	7.5,
                 0.157,	0.30,	7.5,
                 0.123,	0.20,	7.5,
                 0.054,	0.10,	7.5,
                 0.040,	0.07,	7.5,
                 0.021,	0.03,	7.5,
                 0.016,	0.02,	7.5,
                 0.001,	0.01,	7.5,
                 3.530,	30.00,	15.0,
                 2.825,	20.00,	15.0,
                 1.765,	10.00, 15.0,
                 1.460,	7.00,	15.0,
                 0.675,	3.00,	15.0,
                 0.531,	2.00,	15.0,
                 0.275,	1.00,	15.0,
                 0.209,	0.70,	15.0,
                 0.093,	0.30,	15.0,
                 0.079,	0.20,	15.0,
                 0.034,	0.10,	15.0,
                 0.026,	0.07,	15.0,
                 0.013,	0.03,	15.0,
                 0.011,	0.02,	15.0,
                 0.000,	0.01,	15.0,
                 2.335,	30.00,	30.0,
                 2.010,	20.00,	30.0,
                 1.044,	10.00,	30.0,
                 0.970,	7.00,	30.0,
                 0.385,	3.00,	30.0,
                 0.368,	2.00,	30.0,
                 0.151,	1.00,	30.0,
                 0.151,	0.70,	30.0,
                 0.061,	0.30,	30.0,
                 0.048,	0.20,	30.0,
                 0.026,	0.10,	30.0,
                 0.021,	0.07,	30.0,
                 0.007,	0.03,	30.0,
                 0.007,	0.02,	30.0,
                 0.001,	0.01,	30.0,
                 1.175,	30.00,	60.0,
                 1.370,	20.00,	60.0,
                 0.566,	10.00, 60.0,
                 0.713,	7.00,	60.0,
                 0.245,	3.00,	60.0,
                 0.250,	2.00,	60.0,
                 0.104,	1.00,	60.0,
                 0.118,	0.70,	60.0,
                 0.046,	0.30,	60.0,
                 0.038,	0.20,	60.0,
                 0.019,	0.10,	60.0,
                 0.016,	0.07,	60.0,
                 0.006,	0.03,	60.0,
                 0.005,	0.02,	60.0,
                 0.000,	0.01,	60.0),nrow=120,ncol=3,byrow=TRUE)

#-------------------------------------------------------------------------------------------------
xmat<-basia[,2-3]
n<-nrow(basia)


generate_y_Example2 <- function(xi, model, th, sig) {
  
  # Generates data from one of the two models
  # in the design points given by the design xi. 
  
  # Arguments:
  # xi ... a matrix nx2, where n is the number of observations
  # model ... the model (either 0 or 1) used for the simulations
  # th ... the 3D parameter theta used for the simulations
  # sig ... the standard deviation of errors used for the simulation
  
  # Model 0 (competitive inhibition)
  f0 <- function(S, I, th){
    V <- th[1]
    Km <- th[2]
    Kic <- th[3]
    V * S / (Km * (1 + I/Kic) + S)
  }
  
  # Model 1 (non-competitive inhibition)
  f1 <- function(S, I, th){
    V <- th[1]
    Km <- th[2]
    Kin <- th[3]
    V * S / ((Km + S) * (1 + I/Kin))
  }
  
  if (model == 0) y <- f0(xi[, 1], xi[, 2], th) + rnorm(nrow(xi), sd = sig)
  if (model == 1) y <- f1(xi[, 1], xi[, 2], th) + rnorm(nrow(xi), sd = sig)
  
  # y is the simulated vector of observations
  return(y)
}

#-------------------------------------------------------------------------------------------------

decide_Example2 <- function(xi, y, th0i, th1i, model) 
{
  # Makes a decision about the true model based on the minimum sum of squares
  # Competing models are the two competitive and non competitive models
  #
  # Arguments:
  # xi ... the nx2 matrix of design points (with possible replications)
  # y ... the nx1 vector of observations in the design points
  # th0i, th1i ... initial parameters for the optimization (both are 3D vectors)
  
  
  
  # Model 0 (competitive inhibition)
  f0 <- function(S, I, th){
    V <- th[1]
    Km <- th[2]
    Kic <- th[3]
    V * S / (Km * (1 + I/Kic) + S)
  }
  
  # Model 1 (non-competitive inhibition)
  f1 <- function(S, I, th){
    V <- th[1]
    Km <- th[2]
    Kin <- th[3]
    V * S / ((Km + S) * (1 + I/Kin))
  }
  
  # log likelihood function for model 0
  l0 <- function(x, xi, y){
    ((-n/2)*log(2*pi*(sigh0^2)))-(1/2*(sigh0^2))*sum((y - (x[1] * xi[, 1] / (x[2] * (1 + xi[, 2]/x[3]) + xi[, 1])))^2) 
  }
  
  # log likelihood function for model 1
  l1 <- function(x, xi, y){
    ((-n/2)*log(2*pi*(sigh1^2)))-(1/2*(sigh1^2))*sum((y - (x[1] * xi[, 1] / ((x[2] + xi[, 1]) * (1 + xi[, 2]/x[3]))))^2)
  }
  
  # Calculate MLE for Model 0
  lik.func0 <- BBoptim(th0i, l0,xi = xi, y = y, gr=NULL,control=list(ftol=0.001, trace=F, maximize = TRUE), method=1, 
                       lower=c(0.001, 0.001, 0.001), upper=c(Inf, 60, 30), quiet=T)
  par0<-lik.func0$par
  
  
  # Calculate MLE for Model 1 
  lik.func1 <- BBoptim(th1i, l1,xi = xi, y = y, gr=NULL,control=list(ftol=0.001, trace=F, maximize = TRUE), method=1, 
                       lower=c(0.001, 0.001, 0.001), upper=c(Inf, 60, 30), quiet=T)
  
  par1<-lik.func1$par
  
  
  if(model==0)  T.stat <-lik.func0$value - lik.func1$value 
  if(model==1)  T.stat <-lik.func1$value - lik.func0$value
  
  
  return(list(T.stat=T.stat, par0=par0, par1=par1))
}

t0<-decide_Example2(xmat,basia[,1],th0E2,th1E2,0)$T.stat
t1<-decide_Example2(xmat,basia[,1],th0E2,th1E2,1)$T.stat
#-------------------------------------------------------------------------------------------------
set.seed(123456789)
T.AB<- C <- Neg.C <- matrix(0, nrow = N, ncol = 2)


for(i in 1:N){
  for (im in 1:2){
    if (im == 1){
      count01=0 #all undesired cases until a good one
      count02=0 #all negative pars
      repeat{
        y0<- generate_y_Example2(xmat, 0, th0E2,sigh0)
        decide0<-decide_Example2(xmat, y0, th0E2, th1E2, 0) 
        PAR00=decide0$par0
        PAR01=decide0$par1
        if(length(PAR00[PAR00<=0])==0&length(PAR01[PAR01<=0])==0&
           PAR00[1]>0.001&PAR00[1]<Inf&PAR01[1]>0.001&PAR01[1]<Inf&
           PAR00[2]>0.001&PAR00[2]<60& PAR01[2]>0.001&PAR01[2]<60&
           PAR00[3]>0.001&PAR00[3]<30& PAR01[3]>0.001&PAR01[3]<30){
          T.AB[i,im]<-decide0$T.stat
          break
        }else{
          count01=count01+1 
          if(length(PAR00[PAR00<=0])!=0 |length(PAR01[PAR01<=0])!=0 ) count02=count02+1
        }
      }
      C[i,im]<-count01
      Neg.C[i,im]<-count02
    }
    if (im == 2){
      count11=0 #all undesired cases until a good one
      count12=0 #all negative pars
      repeat{
        y1<- generate_y_Example2(xmat, 1, th1E2,sigh1)
        decide1<-decide_Example2(xmat, y1, th0E2, th1E2, 1)
        PAR10=decide1$par0
        PAR11=decide1$par1
        if(length(PAR10[PAR10<=0])==0&length(PAR11[PAR11<=0])==0&
           PAR10[1]>0.001&PAR10[1]<Inf&PAR11[1]>0.001&PAR11[1]<Inf&
           PAR10[2]>0.001&PAR10[2]<60& PAR11[2]>0.001&PAR11[2]<60&
           PAR10[3]>0.001&PAR10[3]<30& PAR11[3]>0.001&PAR11[3]<30){
          T.AB[i,im]<-decide1$T.stat
          break
        }else{
          count11=count11+1
          if(length(PAR10[PAR10<=0])!=0 |length(PAR11[PAR11<=0])!=0 ) count12=count12+1
        }
      }
      C[i,im]<-count11
      Neg.C[i,im]<-count12
    }
  } 
  
}

p.value0<-length(T.AB[,1][ T.AB[,1] <= t0 ])/N
p.value1<-length(T.AB[,2][ T.AB[,2] <= t1 ])/N

p.value0
p.value1
