# Computation of the optimal design based on quadratic design criterion (Hamilton & Watts (1985)) in the log case
# considered model is the encompassing model

rm(list=ls())
library(numDeriv)
library(base)
library(nloptr)
library(MASS)


# initial parameter values for the encompassing model (Log case)
priortheta8=c(6.9897, 3.9799, 3.7380 , 0.8737)

#n=p number of parameters
p <- 4
n <- 4

# The quadratic design criterion (Hamilton & Watts (1985))
outputfunc <- function(xfull){
  
  # xfull: the design, the support points (xmat)  and the weights w1-w4
  x11=xfull[1]
  x21=xfull[2]
  x31=xfull[3]
  x41=xfull[4]
  x12=xfull[5]
  x22=xfull[6]
  x32=xfull[7]
  x42=xfull[8]
  
  w1=xfull[9]
  w2=xfull[10]
  w3=xfull[11]
  w4=xfull[12]
  
  
  xmat=matrix(c(x11,x12,x21,x22,x31,x32,x41,x42),nrow=4,byrow=TRUE)
  
  
  #encompassing log model
  combd.logmodel=function(par){
    theta1=par[1]
    theta2=par[2]
    theta3=par[3]
    theta4=par[4]
    
    eta=log((theta1*xmat[,1])/ ((theta2*(1+(xmat[,2]/theta3)))+xmat[,1]*(1+(((1-theta4)*xmat[,2])/theta3))))
  }
  
  combd.logmodelvec=function(par,xvec){
    theta1=par[1]
    theta2=par[2]
    theta3=par[3]
    theta4=par[4]
    
    eta=log((theta1*xvec[1])/ ((theta2*(1+(xvec[2]/theta3)))+xvec[1]*(1+(((1-theta4)*xvec[2])/theta3)))) #combined log model
  }
  
  
  
  
  V=jacobian(func=combd.logmodel, x=priortheta8, method="Richardson")
  
  
  Vd1=hessian(func=combd.logmodelvec, x=priortheta8, method="Richardson",xvec=xmat[1,])
  Vd2=hessian(func=combd.logmodelvec, x=priortheta8, method="Richardson",xvec=xmat[2,])
  Vd3=hessian(func=combd.logmodelvec, x=priortheta8, method="Richardson",xvec=xmat[3,])
  Vd4=hessian(func=combd.logmodelvec, x=priortheta8, method="Richardson",xvec=xmat[4,])
  
  
  QR <- qr(V)
  U <- qr.Q(QR)
  L <- qr.R(QR)
  L <- ginv(L)
  
  
  A1  <-A2  <-A3  <-A4  <- matrix(0, nrow=4,ncol=4,byrow=T)
  M  <- matrix(0, nrow=4,ncol=4,byrow=T)
  
  ########################################################################################
  # computation of A matrices
  ########################################################################################
  for(i in 1:p){
    if(i==1){
      for(j in 1:p){
        for(k in 1:p){
          
          
          for(r in 1:p){
            for(s in 1:p){
              for(u in 1:p){
                for(t in 1:n){
                  if(t==1){
                    A1[j,k]<- A1[j,k]+(V[t,u]*Vd1[r,s]*L[u,i]*L[r,j]*L[s,k])
                  }
                  if(t==2){
                    A1[j,k]<- A1[j,k]+(V[t,u]*Vd2[r,s]*L[u,i]*L[r,j]*L[s,k])
                  }
                  if(t==3){
                    A1[j,k]<- A1[j,k]+(V[t,u]*Vd3[r,s]*L[u,i]*L[r,j]*L[s,k])
                  } 
                  if(t==4){
                    A1[j,k]<- A1[j,k]+(V[t,u]*Vd4[r,s]*L[u,i]*L[r,j]*L[s,k])
                  } 
                  
                }
              }
              
            }
          }
          
          
          
        }
      }  
    }
    if(i==2){
      for(j in 1:p){
        for(k in 1:p){
          
          
          for(r in 1:p){
            for(s in 1:p){
              for(u in 1:p){
                for(t in 1:n){
                  if(t==1){
                    A2[j,k]<- A2[j,k]+(V[t,u]*Vd1[r,s]*L[u,i]*L[r,j]*L[s,k])
                  }
                  if(t==2){
                    A2[j,k]<- A2[j,k]+(V[t,u]*Vd2[r,s]*L[u,i]*L[r,j]*L[s,k])
                  }
                  if(t==3){
                    A2[j,k]<- A2[j,k]+(V[t,u]*Vd3[r,s]*L[u,i]*L[r,j]*L[s,k])
                  } 
                  if(t==4){
                    A2[j,k]<- A2[j,k]+(V[t,u]*Vd4[r,s]*L[u,i]*L[r,j]*L[s,k])
                  } 
                  
                }
              }
              
            }
          }
          
          
          
        }
      }   
    }
    if(i==3){
      for(j in 1:p){
        for(k in 1:p){
          
          
          for(r in 1:p){
            for(s in 1:p){
              for(u in 1:p){
                for(t in 1:n){
                  if(t==1){
                    A3[j,k]<- A3[j,k]+(V[t,u]*Vd1[r,s]*L[u,i]*L[r,j]*L[s,k])
                  }
                  if(t==2){
                    A3[j,k]<- A3[j,k]+(V[t,u]*Vd2[r,s]*L[u,i]*L[r,j]*L[s,k])
                  }
                  if(t==3){
                    A3[j,k]<- A3[j,k]+(V[t,u]*Vd3[r,s]*L[u,i]*L[r,j]*L[s,k])
                  } 
                  if(t==4){
                    A3[j,k]<- A3[j,k]+(V[t,u]*Vd4[r,s]*L[u,i]*L[r,j]*L[s,k])
                  } 
                  
                }
              }
              
            }
          }
          
          
          
        }
      }   
    }
    
    if(i==4){
      for(j in 1:p){
        for(k in 1:p){
          
          
          for(r in 1:p){
            for(s in 1:p){
              for(u in 1:p){
                for(t in 1:n){
                  if(t==1){
                    A4[j,k]<- A4[j,k]+(V[t,u]*Vd1[r,s]*L[u,i]*L[r,j]*L[s,k])
                  }
                  if(t==2){
                    A4[j,k]<- A4[j,k]+(V[t,u]*Vd2[r,s]*L[u,i]*L[r,j]*L[s,k])
                  }
                  if(t==3){
                    A4[j,k]<- A4[j,k]+(V[t,u]*Vd3[r,s]*L[u,i]*L[r,j]*L[s,k])
                  } 
                  if(t==4){
                    A4[j,k]<- A4[j,k]+(V[t,u]*Vd4[r,s]*L[u,i]*L[r,j]*L[s,k])
                  } 
                  
                }
              }
              
            }
          }
          
          
          
        }
      }   
    }
  }
  
  
  #-------------------------------------------
  for(i in 1:p){
    for(j in 1:p){
      M[i,j]<- 
        A1[1,i]*A1[1,j]+A1[1,j]*A1[i,1]+A1[1,1]*A1[i,j]+
        A1[1,i]*A2[2,j]+A2[1,j]*A1[i,2]+A2[2,1]*A1[i,j]+
        A1[1,i]*A3[3,j]+A3[1,j]*A1[i,3]+A3[3,1]*A1[i,j]+
        A1[1,i]*A4[4,j]+A4[1,j]*A1[i,4]+A4[4,1]*A1[i,j]+ 
        
        A2[2,i]*A1[1,j]+A1[2,j]*A2[i,1]+A1[1,2]*A2[i,j]+
        A2[2,i]*A2[2,j]+A2[2,j]*A2[i,2]+A2[2,2]*A2[i,j]+
        A2[2,i]*A3[3,j]+A3[2,j]*A2[i,3]+A3[3,2]*A2[i,j]+
        A2[2,i]*A4[4,j]+A4[2,j]*A2[i,4]+A4[4,2]*A2[i,j]+
        
        A3[3,i]*A1[1,j]+A1[3,j]*A3[i,1]+A1[1,3]*A3[i,j]+
        A3[3,i]*A2[2,j]+A2[3,j]*A3[i,2]+A2[2,3]*A3[i,j]+
        A3[3,i]*A3[3,j]+A3[3,j]*A3[i,3]+A3[3,3]*A3[i,j]+
        A3[3,i]*A4[4,j]+A4[3,j]*A3[i,4]+A4[4,3]*A3[i,j]+
        
        A4[4,i]*A1[1,j]+A1[4,j]*A4[i,1]+A1[1,4]*A4[i,j]+
        A4[4,i]*A2[2,j]+A2[4,j]*A4[i,2]+A2[2,4]*A4[i,j]+
        A4[4,i]*A3[3,j]+A3[4,j]*A4[i,3]+A3[3,4]*A4[i,j]+
        A4[4,i]*A4[4,j]+A4[4,j]*A4[i,4]+A4[4,4]*A4[i,j]
      
    }
  }
  
  #M
  
  R <- 120*diag(c(w1,w2,w3,w4))
  TRC <- sum(diag(ginv(R)%*%M))
  k02 <- ((0.5128)*qchisq(.95, df=4)/12)
  
  crit=det(R)^(-1/2)*(det(t(V)%*%V))^(-1/2)*(1+k02*TRC)
  return(crit)  
}
#-------------------------------------------------
#-------------------------------------------------


eval_grad_f0 <- function( xfull1 ) {
  jacob=jacobian(func=outputfunc, x=xfull1, method="Richardson")
  jacob
}

#-------------------------------------------------

eval_g_eq <- function( xfull ) {
  constr <- c( xfull[9]+xfull[10]+xfull[11]+xfull[12]-1 )
  grad <- c( 0,0,0,0,0,0,0,0,1,1,1,1 )
  return( list( "constraints"=constr, "jacobian"=grad ) )
}
#-------------------------------------------------
local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel" = 1.0e-7 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-7,
              "maxeval" = 1000,
              "local_opts" = local_opts )
#-------------------------------------------------
# initial values
zz <- c(0.02,30,0.02,30,0,0,60,60,0.25,0.25,0.25,0.25)

# 
res <- nloptr( x0=zz,eval_f=outputfunc ,eval_grad_f=eval_grad_f0,lb = c(0.02,0.02,0.02,0.02,0,0,0,0,0,0,0,0),
               ub=c(30,30,30,30,60,60,60,60,1,1,1,1),eval_g_eq =eval_g_eq ,opts=opts)


print( res )
cbind(round(res$solution[1:4],3),round(res$solution[5:8],3),round(res$solution[9:12],3))
