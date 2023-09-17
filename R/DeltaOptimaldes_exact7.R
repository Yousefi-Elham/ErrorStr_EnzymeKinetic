# Section 4.3, A simulation study of discriminating designs, subsection 4.3.1 Exact designs, N = 7
rm(list=ls())
par(mfrow=c(1,1))

th0E2 <- c(6.0645, 3.2799, 3.3153) #estimated values from nls on original obs. for competitive log model
s0E2 <- c(0.9260 , 0.7288,  0.6041) #std.dev values from nls on original obs. for competitive log model
th1E2 <- c(12.0125, 8.5359, 5.6638) #estimated values from nls on original obs. for non competitive log model
s1E2 <- c(2.0553, 1.5721, 0.8879) #std.dev values from nls on original obs. for non competitive log model


# discretization of the 2D design space
epsilon <- 0.02
x1 <- c(epsilon,1:30)
x2 <- 0:60
XE2 <- expand.grid(x1,x2)


##################################################################################################################
od_delta_KL_Example2 <- function(th0, th0l, th0u, th1, th1l, th1u, X, n, t.max=60) {
  
  # The KL-exchange algorithm for computing delta-optimal designs
  #
  # Arguments:
  # th0 ... the nominal value of theta_0 (a 3D vector)
  # th0l, th0u ... 3D vectors representing the lower and upper bounds for theta_0
  # th1 ... the nominal value of theta_1 (a 3D vector)
  # th1l, th1u ... 3D vectors representing the lower and upper bounds for theta_1
  # X ... the N x 2 matrix, the N-point discretization of the 2D design space
  # n ... the required size of the experiment (the number of observations)
  # t.max ... the computation time
  
  # The size of the design space   
  N <- nrow(X)
  
  # Calculate F0, F1, a0, a1 for all design points
  F0.full <- matrix(0, nrow = N, ncol = 3)
  a0.full <- rep(0, N)
  F1.full <- matrix(0, nrow = N, ncol = 3)
  a1.full <- rep(0, N)
  
  for (i in 1:N) {
    x1 <- X[i, 1]
    x2 <- X[i, 2]
    
    a <- log(th0[1] * x1) 
    b <- 1 + x2 / th0[3]
    d <- log(b * th0[2] + x1)
    F0.full[i, 1] <- 1 / th0[1]
    F0.full[i, 2] <- -b /(b * th0[2] + x1) 
    F0.full[i, 3] <- (th0[2]*x2)/((th0[3]^2) *(b * th0[2] + x1))
    a0.full[i] <- a - d - sum(F0.full[i,] * th0)
    
    a <- log(th1[1] * x1) 
    b <- 1 + x2 / th1[3]
    d <- log(b) + log(th1[2] + x1)
    F1.full[i, 1] <- 1 / th1[1]
    F1.full[i, 2] <- -1 / (th1[2] + x1)
    F1.full[i, 3] <- x2 / (th1[3]^2 * b)
    a1.full[i] <- a - d - sum(F1.full[i,] * th1)
  }
  
  # The square of the delta-criterion for the design given by weights w
  crit_delta_sq <- function(w) 
  {
    sel <- rep(1:N, w)
    F0 <- F0.full[sel,]; F1 <- F1.full[sel,]
    a0 <- a0.full[sel]; a1 <- a1.full[sel]
    r <- bvls::bvls(cbind(F0, -F1), a1 - a0, c(th0l, th1l), c(th0u, th1u))$deviance
    
  }
  
  start <- as.numeric(proc.time()[3])
  info <- paste("Running od_delta_KL_Example2 for cca", t.max, "seconds")
  info <- paste(info, " starting at ", Sys.time(), ".", sep = "")
  print(info, quote = FALSE)
  
  next.sec <- 0
  n.ex <- 0
  n.rest <- 0
  
  finish.all <- FALSE
  crit.best <- -1
  
  while (!finish.all) {
    n.rest <- n.rest + 1
    w <- rep(0, N)
    for (k in 1:n) {
      i <- sample(1:N, 1)
      w[i] <- w[i] + 1
    }
    crit.w <- crit_delta_sq(w)
    
    finish.all <- finish <- as.numeric(proc.time()[3]) > start + t.max
    while (!finish) {
      tm <- as.numeric(proc.time()[3]) - start
      if (tm > next.sec) {
        info <- paste("Time:", round(tm, 1), "secs, Best value:", crit.best)
        print(info, quote = FALSE)
        next.sec <- ceiling(tm)
      }
      
      non.supp <- w < 1e-09
      Kd.fun <- runif(N)
      Kd.fun[non.supp] <- Inf
      Kord <- order(Kd.fun)
      Kact <- N - sum(non.supp)
      Kind <- Kord[1:Kact]
      Lind <- sample(1:N)
      
      imp <- FALSE
      for (iL in Lind) {
        for (iK in Kind) {
          w.temp <- w
          w.temp[iL] <- w[iL] + 1
          w.temp[iK] <- w[iK] - 1
          crit.temp <- crit_delta_sq(w.temp)
          
          if (crit.temp > crit.w + 1e-12) {
            w <- w.temp
            crit.w <- crit.temp
            n.ex <- n.ex + 1
            imp <- TRUE
            break
          }
        }
        if (imp) break
      }
      
      if (as.numeric(proc.time()[3]) > start + t.max) finish.all <- TRUE
      if (finish.all || !imp) finish <- TRUE
    }
    
    if (crit.w > crit.best) {
      w.best <- w
      crit.best <- crit.w
    }
    
    plot(X[, 1], X[, 2], type = "n", main = paste("Best design after restart no.", n.rest)); grid()
    points(X[as.logical(w.best), 1], X[as.logical(w.best), 2], pch = 19,
           cex = sqrt(w.best[as.logical(w.best)]))
  }
  
  t.act <- round(as.numeric(proc.time()[3]) - start, 2)
  info <- paste("od_delta_KL_Example2 finished after", t.act, "seconds at", Sys.time())
  print(info, quote = FALSE)
  
  info <- paste("with", n.rest, "restarts and", n.ex, "exchanges.")
  print(info, quote = FALSE)
  
  plot(X[, 1], X[, 2], type = "n", main = "Best design found"); grid()
  points(X[as.logical(w.best),1], X[as.logical(w.best),2], pch = 19, 
         cex = sqrt(w.best[as.logical(w.best)]))
  
  
  X.design <- cbind(X,w.best)
  w.pos=which(X.design[,3]>0)
  best.design <- X.design[w.pos,]
  plot(best.design[, 1], best.design[, 2], type = "n",xlab="X1",ylab="X2", main = "Best design found"); grid()
  for (j in 1:nrow(best.design)) {
    text(best.design[j, 1], best.design[j, 2], best.design[j, 3])
  }
  
  # Output values: 
  # w.best ... the weights of the best design found
  # xi.best ... the best design found
  # delta.sq.best ... the criterion value of xi.best
  # t.act ... the actual computation time
  list(w.best = w.best, xi.best = X[rep(1:N, w.best),],
       delta.sq.best = crit.best, t.act = t.act)
}

#-------------------------------------------------------------------------------------------------
# Computation of delta optimal designs, which are not directly presented in the paper but used for simulations.


#delta1, r=1
set.seed(123456789) 
delta1 <- od_delta_KL_Example2(th0E2, th0E2 - s0E2, th0E2 + s0E2, 
                               th1E2, th1E2 - s1E2, th1E2 + s1E2, XE2, n = 7, t.max = 120)

#delta2, r=2
set.seed(123456789) 
delta2 <- od_delta_KL_Example2(th0E2, th0E2 - 2*s0E2, th0E2 + 2*s0E2, 
                               th1E2, th1E2 - 2*s1E2, th1E2 + 2*s1E2, XE2, n = 7, t.max = 120)

#delta3, r=3
set.seed(123456789)
delta3 <- od_delta_KL_Example2(th0E2, th0E2 - 3*s0E2, th0E2 + 3*s0E2, 
                               th1E2, th1E2 - 3*s1E2, th1E2 + 3*s1E2, XE2, n = 7, t.max = 120)

#delta4, r=4
set.seed(123456789)
delta4 <- od_delta_KL_Example2(th0E2, th0E2 - 4*s0E2, th0E2 + 4*s0E2, 
                               th1E2, th1E2 - 4*s1E2, th1E2 + 4*s1E2, XE2, n = 7, t.max = 120)

#-------------------------------------------------------------------------------------------------
# the bounds for parameters for r=5 contain negative values so we used 3 alternatives
th0E2 - 5*s0E2
th0E2 + 5*s0E2

th1E2 - 5*s1E2
th1E2 + 5*s1E2

# alternative a) to avoid negative lower bands, when the tuning parameter r=5
set.seed(123456789)
delta5 <- od_delta_KL_Example2(th0E2, c(1.4345, 0.00,  0.2948), th0E2 + 5*s0E2, 
                               th1E2, th1E2 - 5*s1E2, th1E2 + 5*s1E2, XE2, n = 7, t.max = 120)

# alternative b) to avoid negative lower bands, when the tuning parameter r=5
set.seed(123456789) 
delta51 <- od_delta_KL_Example2(th0E2, c(1.4345, 0.00,  0.2948), c(10.6945,  7.288,  6.3358), 
                                th1E2, th1E2 - 5*s1E2, th1E2 + 5*s1E2, XE2, n = 7, t.max = 120)

# alternative c) to avoid negative lower bands, when the tuning parameter r=5
set.seed(123456789)
delta52 <- od_delta_KL_Example2(th0E2, c(1.4345, 1.079830,  0.2948), c(10.6945,  9.962444,  6.3358), 
                                th1E2, th1E2 - 5*s1E2, th1E2 + 5*s1E2, XE2, n = 7, t.max = 120)

#-------------------------------------------------------------------------------------------------
# the bounds for parameters for r=10 contain negative values so we used 3 alternatives
th0E2 - 10*s0E2
th0E2 + 10*s0E2

th1E2 - 10*s1E2
th1E2 + 10*s1E2

#alternative a) to avoid negative lower bands, when the tuning parameter r=10
set.seed(123456789)
delta6 <- od_delta_KL_Example2(th0E2, c(0.00, 0.00, 0.00), th0E2 + 10*s0E2, 
                               th1E2, c(0.00, 0.00, 0.00), th1E2 + 10*s1E2, XE2, n = 7, t.max = 120)

#alternative b) to avoid negative lower bands, when the tuning parameter r=10
set.seed(123456789)
delta61 <- od_delta_KL_Example2(th0E2, c(0.00, 0.00, 0.00), c(18.52,  14.576,  12.082), 
                                th1E2, c(0.00, 0.00, 0.00), c(41.106, 31.442,  17.758), XE2, n = 7, t.max = 120)

#alternative c) to avoid negative lower bands, when the tuning parameter r=10
set.seed(123456789)
delta62 <- od_delta_KL_Example2(th0E2, c(1.3172328, 0.3555085, 0.5360061), c(27.92078, 30.26016, 20.50576), 
                                th1E2, c(5.557143,  3.634513,  1.949351), c(25.96661,  20.04714,  16.45605),
                                XE2, n = 7, t.max = 120)

#-------------------------------------------------------------------------------------------------
# the bounds for parameters for r=15 contain negative values so we used 3 alternatives
th0E2 - 15*s0E2
th0E2 + 15*s0E2

th1E2 - 15*s1E2
th1E2 + 15*s1E2

#alternative a) to avoid negative lower bands, when the tuning parameter r=15
set.seed(123456789)
delta7 <- od_delta_KL_Example2(th0E2, c(0.00, 0.00, 0.00), th0E2 + 15*s0E2, 
                               th1E2, c(0.00, 0.00, 0.00), th1E2 + 15*s1E2, XE2, n = 7, t.max = 120)

#alternative b) to avoid negative lower bands, when the tuning parameter r=15
set.seed(123456789)
delta71 <- od_delta_KL_Example2(th0E2, c(0.00, 0.00, 0.00), c(27.780,  21.864,  18.123), 
                                th1E2, c(0.00, 0.00, 0.00), c(61.659,  47.163,  26.637), XE2, n = 7, t.max = 120)

#alternative c) to avoid negative lower bands, when the tuning parameter r=15
set.seed(123456789)
delta72 <- od_delta_KL_Example2(th0E2, c(0.6138981, 0.1170428, 0.2155228), c(59.90923, 91.91290, 50.99792), 
                                th1E2, c(3.779729, 2.371618,  1.143619), c(38.17738, 30.72232, 28.05011),
                                XE2, n = 7, t.max = 120)


##################################################################################################
des.delta1<-delta1$xi.best
des.delta2<-delta2$xi.best
des.delta3<-delta3$xi.best
des.delta4<-delta4$xi.best

des.delta5<-delta5$xi.best
des.delta51<-delta51$xi.best
des.delta52<-delta52$xi.best


des.delta6<-delta6$xi.best
des.delta61<-delta61$xi.best
des.delta62<-delta62$xi.best

des.delta7<-delta7$xi.best
des.delta71<-delta71$xi.best
des.delta72<-delta72$xi.best

des.delta1
des.delta2
des.delta3
des.delta4

des.delta5
des.delta51
des.delta52


des.delta6
des.delta61
des.delta62

des.delta7
des.delta71
des.delta72
