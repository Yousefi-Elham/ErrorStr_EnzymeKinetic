# Comparison of D and Ds optimal designs using D and Ds efficiencies (Table 4)

rm(list=ls())
library(numDeriv)


# prior parameter values resulting from parameter estimation 
# A dataset with 120 observations has been used for this purpose (Tables 1 & 2)

#prior values in the standard case
priortheta.nc=c(8.6957,8.0664,12.0566,0) 
priortheta.c=c(7.2976,4.3860,2.5821,1)  
priortheta.ds.comb=c(7.4253,4.6808,3.0581,0.9636)  

#prior values in the log case
priortheta.lognc=c(12.0125, 8.5359, 5.6638,0) 
priortheta.logc=c(6.0645, 3.2799, 3.3153,1)  
priortheta.ds.logcomb8=c(6.9897, 3.9799, 3.7380 , 0.8737) 

#--------------------------------------------------------------------------------------------------------------
# encompassing model as a function of parameters and design points
combd.model=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  theta4=par[4]
  
  eta=(theta1*xval[,1])/ ((theta2*(1+(xval[,2]/theta3)))+xval[,1]*(1+(((1-theta4)*xval[,2])/theta3))) #combined model
}


# non competitive model as a function of parameters and design points
noncomp.model=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  
  eta=(theta1*xval[,1]) / ((theta2+xval[,1])*(1+(xval[,2]/theta3))) 
}


# competitive model as a function of parameters and design points
comp.model=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  
  eta=(theta1*xval[,1]) / ((theta2*(1+(xval[,2]/theta3)))+xval[,1]) 
}
#--------------------------------------------------------------------------------------------------------------
# encompassing log model as a function of parameters and design points
combd.logmodel=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  theta4=par[4]
  
  eta=log((theta1*xval[,1])/ ((theta2*(1+(xval[,2]/theta3)))+xval[,1]*(1+(((1-theta4)*xval[,2])/theta3)))) 
}


# non competitive log model as a function of parameters and design points
noncomp.logmodel=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  
  eta=log((theta1*xval[,1]) / ((theta2+xval[,1])*(1+(xval[,2]/theta3)))) 
}

# competitive log model as a function of parameters and design points
comp.logmodel=function(par,xval){
  theta1=par[1]
  theta2=par[2]
  theta3=par[3]
  
  eta=log((theta1*xval[,1]) / ((theta2*(1+(xval[,2]/theta3)))+xval[,1]))
}
#--------------------------------------------------------------------------------------------------------------

# Optimal designs of Table 3, D and Ds-optimal designs (Standard case)

weight.6 <- c(0.086,0.208,0.206,0.500)
weight.7 <- c(0.027,0.088,0.371,0.514)

des1.4DNC<-data.frame(matrix(c(30,0,5.223,0,30,12.045,5.223,12.045),nrow=4,byrow=TRUE),rep(0.25,4),priortheta.nc)
des1.3DNC<-data.frame(matrix(c(30,0,5.223,0,30,12.045),nrow=3,byrow=TRUE),rep(1/3,3),priortheta.nc[1:3])
des1.4DC<-data.frame(matrix(c(30,0,3.348,0,30,20.297,7.902,7.137),nrow=4,byrow=TRUE),rep(0.25,4),priortheta.c)
des1.3DC<-data.frame(matrix(c(30,0,3.348,0,30,20.297),nrow=3,byrow=TRUE),rep(1/3,3),priortheta.c[1:3])
des1.4DComb<-data.frame(matrix(c(30,0,3.616,0,30,18.290,7.500,7.584),nrow=4,byrow=TRUE),rep(0.25,4),priortheta.ds.comb)
des1.DsN<-data.frame(matrix(c(30,0,3.884,0,30,16.952,3.884,16.952),nrow=4,byrow=TRUE),weight.6, priortheta.nc)
des1.DsC<-data.frame(matrix(c(30,0,2.545,0,30,28.550,7.098,8.253),nrow=4,byrow=TRUE),weight.7, priortheta.c)

des2.4DNC<-data.frame(matrix(c(30,0,5.223,0,30,12.045,5.223,12.045),nrow=4,byrow=TRUE),rep(0.25,4))
des2.3DNC<-data.frame(matrix(c(30,0,5.223,0,30,12.045),nrow=3,byrow=TRUE),rep(1/3,3))
des2.4DC<-data.frame(matrix(c(30,0,3.348,0,30,20.297,7.902,7.137),nrow=4,byrow=TRUE),rep(0.25,4))
des2.3DC<-data.frame(matrix(c(30,0,3.348,0,30,20.297),nrow=3,byrow=TRUE),rep(1/3,3))
des2.4DComb<-data.frame(matrix(c(30,0,3.616,0,30,18.290,7.500,7.584),nrow=4,byrow=TRUE),rep(0.25,4))
des2.DsN<-data.frame(matrix(c(30,0,3.884,0,30,16.952,3.884,16.952),nrow=4,byrow=TRUE),weight.6)
des2.DsC<-data.frame(matrix(c(30,0,2.545,0,30,28.550,7.098,8.253),nrow=4,byrow=TRUE),weight.7)


#--------------------------------------------------------------------------------------------------------------
# Optimal designs of Table 3, D and Ds-optimal designs (Log case)

des1.4DNlog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),rep(0.25,4), priortheta.lognc)
des1.3DNlog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),rep(0.25,4), c(priortheta.lognc[1:3],-1) )
des1.4DClog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),rep(0.25,4), priortheta.logc)
des1.3DClog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60),nrow=3,byrow=TRUE),rep(1/3,3), priortheta.logc[1:3])
des1.4DComblog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),rep(0.25,4), priortheta.ds.logcomb8)
des1.DsNlog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),rep(0.25,4), priortheta.lognc)
weight.10<-c(0.017,0.173,0.327,0.483)
des1.DsClog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),weight.10, priortheta.logc)

des2.4DNlog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),rep(0.25,4))
des2.3DNlog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),rep(0.25,4))
des2.4DClog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),rep(0.25,4))
des2.3DClog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60),nrow=3,byrow=TRUE),rep(1/3,3))
des2.4DComblog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),rep(0.25,4))
des2.DsNlog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),rep(0.25,4))
weight.10<-c(0.017,0.173,0.327,0.483)
des2.DsClog<-data.frame(matrix(c(0.02,0,30.00,0,0.02,60,30.00,60),nrow=4,byrow=TRUE),weight.10)

#--------------------------------------------------------------------------------------------------------------
# create a list of all optimal designs, above
design1=list(des1.4DNC=des1.4DNC,des1.3DNC=des1.3DNC,des1.4DC=des1.4DC,des1.3DC=des1.3DC,des1.4DComb=des1.4DComb,des1.DsN=des1.DsN,des1.DsC=des1.DsC,
             des1.4DNlog=des1.4DNlog,des1.3DNlog=des1.3DNlog,des1.4DClog=des1.4DClog,des1.3DClog=des1.3DClog,des1.4DComblog=des1.4DComblog,
             des1.DsNlog=des1.DsNlog,des1.DsClog=des1.DsClog)

design2=list(des2.4DNC=des2.4DNC,des2.3DNC=des2.3DNC,des2.4DC=des2.4DC,des2.3DC=des2.3DC,des2.4DComb=des2.4DComb,des2.DsN=des2.DsN,des2.DsC=des2.DsC,
             des2.4DNlog=des2.4DNlog,des2.3DNlog=des2.3DNlog,des2.4DClog=des2.4DClog,des2.3DClog=des2.3DClog,des2.4DComblog=des2.4DComblog,
             des2.DsNlog=des2.DsNlog,des2.DsClog=des2.DsClog)

NN=length(design1)


# The setting used for comparison in the D and Ds efficiencies
CASE=c(rep("stan.case",7),rep("log.case",7))
method=c(rep("D.OPT",5),rep("Ds.OPT",2),rep("D.OPT",5),rep("Ds.OPT",2))
MODEL=c("combd","noncomp","combd","comp","combd","combd","combd","combd","noncomp","combd","comp","combd","combd","combd")

#--------------------------------------------------------------------------------------------------------------
# Construction of the Information Matrices as a function of design, the respective model and the parameter estimates 

# Full information matrix used for D-optimality
M.func<-function(des.x,model,parm){
  F.mat<-jacobian(func=model, x=parm, method="Richardson",xval=des.x[,-3])
  w.mat<-diag(c(des.x[,3]),nrow=length(c(des.x[,3])))
  MM=det(t(F.mat)%*%w.mat%*%F.mat)
  return(MM)
}


# Information Matrix of the nuisance parameters, used in Ds-optimality procedure
Ms.func<-function(model1,model2,des.x,parm){
  F.mat<-jacobian(func=model1, x=parm, method="Richardson",xval=des.x[,-3])
  w.mat<-diag(c(des.x[,3]),nrow=length(c(des.x[,3])))
  m.mat<-det(t(F.mat)%*%w.mat%*%F.mat)
  
  F22.mat<-jacobian(func=model2, x=parm[1:3], method="Richardson",xval=des.x[,-3])
  m22.mat<-det(t(F22.mat)%*%w.mat%*%F22.mat)
  return(m.mat/m22.mat)
}


#--------------------------------------------------------------------------------------------------------------
# Efficiency formulation

effs.func<-function(ms1,ms2,s){
  return(((ms1/ms2)^(1/s))*100)
}

#--------------------------------------------------------------------------------------------------------------
RESULT=matrix(0,NN,NN)

for(j in 1:NN){
  
  
  #--------------------------------------------------------------------------------------------------------------
  #                                                            standard.case
  #--------------------------------------------------------------------------------------------------------------
  if(CASE[j]=="stan.case"){
    #**********************************************************************
    #                            D.OPT
    #**********************************************************************
    if(method[j]=="D.OPT"){
      
      
      if(MODEL[j]=="combd"){
        design.des=data.frame(design1[j]) 
        if(j==9){
          parmmm1=design.des[,4]
          parmmm=parmmm1[1:3]
        }else{
          parmmm=design.des[,4]
        }
        design.des44=design.des[,-4]
        WW1=M.func(design.des44,combd.model,parmmm)
        for(i in 1:NN){
          design.up=data.frame(design2[i])
          RESULT[i,j]=effs.func(M.func(design.up,combd.model,parmmm),WW1,4)
        }
      }
      
      
      if(MODEL[j]=="noncomp"){
        design.des=data.frame(design1[j])
        if(j==9){
          parmmm1=design.des[,4]
          parmmm=parmmm1[1:3]
        }else{
          parmmm=design.des[,4]
        }
        design.des44=design.des[,-4]
        WW1=M.func(design.des44,noncomp.model,parmmm)
        for(i in 1:NN){
          design.up=data.frame(design2[i])
          RESULT[i,j]=effs.func(M.func(design.up,noncomp.model,parmmm),WW1,3)
        }
      }
      
      
      
      if(MODEL[j]=="comp"){
        design.des=data.frame(design1[j])
        if(j==9){
          parmmm1=design.des[,4]
          parmmm=parmmm1[1:3]
        }else{
          parmmm=design.des[,4]
        }
        design.des44=design.des[,-4]
        WW1=M.func(design.des44,comp.model,parmmm)
        for(i in 1:NN){
          design.up=data.frame(design2[i])
          RESULT[i,j]=effs.func(M.func(design.up,comp.model,parmmm),WW1,3)
        }
      }
    }
    #**********************************************************************
    #                            Ds.OPT
    #**********************************************************************
    
    if(method[j]=="Ds.OPT"){
      
      
      if(j==6){
        design.des=data.frame(design1[j])
        if(j==9){
          parmmm1=design.des[,4]
          parmmm=parmmm1[1:3]
        }else{
          parmmm=design.des[,4]
        }
        design.des44=design.des[,-4]
        WW1=Ms.func(combd.model,noncomp.model,design.des44,parmmm)
        for(i in 1:NN){
          design.up=data.frame(design2[i])
          RESULT[i,j]=effs.func(Ms.func(combd.model,noncomp.model,design.up,parmmm),WW1,1)
        }
      }
      
      if(j==7){
        design.des=data.frame(design1[j])
        if(j==9){
          parmmm1=design.des[,4]
          parmmm=parmmm1[1:3]
        }else{
          parmmm=design.des[,4]
        }
        design.des44=design.des[,-4]
        WW1=Ms.func(combd.model,comp.model,design.des44,parmmm)
        for(i in 1:NN){
          design.up=data.frame(design2[i])
          RESULT[i,j]=effs.func(Ms.func(combd.model,comp.model,design.up,parmmm),WW1,1)
        }
      }
    }
  }
  
  #--------------------------------------------------------------------------------------------------------------
  #                                                            log.case
  #--------------------------------------------------------------------------------------------------------------
  
  if(CASE[j]=="log.case"){
    
    #**********************************************************************
    #                            D.OPT
    #**********************************************************************
    if(method[j]=="D.OPT"){
      
      
      if(MODEL[j]=="combd"){
        design.des=data.frame(design1[j])
        if(j==9){
          parmmm1=design.des[,4]
          parmmm=parmmm1[1:3]
        }else{
          parmmm=design.des[,4]
        }
        design.des44=design.des[,-4]
        WW1=M.func(design.des44,combd.logmodel,parmmm)
        for(i in 1:NN){
          design.up=data.frame(design2[i])
          RESULT[i,j]=effs.func(M.func(design.up,combd.logmodel,parmmm),WW1,4)
        }
      }
      
      
      if(MODEL[j]=="noncomp"){
        design.des=data.frame(design1[j])
        if(j==9){
          parmmm1=design.des[,4]
          parmmm=parmmm1[1:3]
        }else{
          parmmm=design.des[,4]
        }
        design.des44=design.des[,-4]
        WW1=M.func(design.des44,noncomp.logmodel,parmmm)
        for(i in 1:NN){
          design.up=data.frame(design2[i])
          RESULT[i,j]=effs.func(M.func(design.up,noncomp.logmodel,parmmm),WW1,3)
        }
      }
      
      
      
      if(MODEL[j]=="comp"){
        design.des=data.frame(design1[j])
        if(j==9){
          parmmm1=design.des[,4]
          parmmm=parmmm1[1:3]
        }else{
          parmmm=design.des[,4]
        }
        design.des44=design.des[,-4]
        WW1=M.func(design.des44,comp.logmodel,parmmm)
        for(i in 1:NN){
          design.up=data.frame(design2[i])
          RESULT[i,j]=effs.func(M.func(design.up,comp.logmodel,parmmm),WW1,3)
        }
      }
    }
    #**********************************************************************
    #                            Ds.OPT
    #**********************************************************************
    
    if(method[j]=="Ds.OPT"){
      
      
      if(j==13){
        design.des=data.frame(design1[j])
        if(j==9){
          parmmm1=design.des[,4]
          parmmm=parmmm1[1:3]
        }else{
          parmmm=design.des[,4]
        }
        design.des44=design.des[,-4]
        WW1=Ms.func(combd.logmodel,noncomp.logmodel,design.des44,parmmm)
        for(i in 1:NN){
          design.up=data.frame(design2[i])
          RESULT[i,j]=effs.func(Ms.func(combd.logmodel,noncomp.logmodel,design.up,parmmm),WW1,1)
        }
      }
      
      if(j==14){
        design.des=data.frame(design1[j])
        if(j==9){
          parmmm1=design.des[,4]
          parmmm=parmmm1[1:3]
        }else{
          parmmm=design.des[,4]
        }
        design.des44=design.des[,-4]
        WW1=Ms.func(combd.logmodel,comp.logmodel,design.des44,parmmm)
        for(i in 1:NN){
          design.up=data.frame(design2[i])
          RESULT[i,j]=effs.func(Ms.func(combd.logmodel,comp.logmodel,design.up,parmmm),WW1,1)
        }
      }
    }
  }
  
}
names1=c("4DN.stan","3DN.stan","4DC.stan","3DC.stan","4DCOM.stan","DSN.stan","DSC.stan")
names2=c("4DN.log","3DN.log","4DC.log","3DC.log","4DCOM.log","DSN.log","DSC.log")

NAME=c(names1,names2)
RESULT.final=round(RESULT,2)

colnames(RESULT.final)=NAME
rownames(RESULT.final)=NAME

View(RESULT.final)






