
# 120 observations from Bogacka et. al (2011)
basia = read.table(header=TRUE, text="
                   y x1 x2
                   6.640	30.00	0.0
                   6.040	20.00	0.0
                   4.755	10.00	0.0
                   4.105	7.00	0.0
                   3.090	3.00	0.0
                   2.505	2.00	0.0
                   1.455	1.00	0.0
                   1.265	0.70	0.0
                   0.646	0.30	0.0
                   0.525	0.20	0.0
                   0.367	0.10	0.0
                   0.206	0.07	0.0
                   0.088	0.03	0.0
                   0.069	0.02	0.0
                   0.009	0.01	0.0
                   6.085	30.00	0.9
                   5.495	20.00	0.9
                   3.955	10.00	0.9
                   3.540	7.00	0.9
                   2.380	3.00	0.9
                   1.905	2.00	0.9
                   1.185	1.00	0.9
                   0.812	0.70	0.9
                   0.413	0.30	0.9
                   0.341	0.20	0.9
                   0.169	0.10	0.9
                   0.098	0.07	0.9
                   0.062	0.03	0.9
                   0.045	0.02	0.9
                   0.006	0.01	0.9
                   6.030	30.00	2.8
                   5.380	20.00	2.8
                   3.930	10.00	2.8
                   3.320	7.00	2.8
                   2.090	3.00	2.8
                   1.550	2.00	2.8
                   0.865	1.00	2.8
                   0.692	0.70	2.8
                   0.350	0.30	2.8
                   0.245	0.20	2.8
                   0.135	0.10	2.8
                   0.104	0.07	2.8
                   0.043	0.03	2.8
                   0.026	0.02	2.8
                   0.003	0.01	2.8
                   5.570	30.00	3.7
                   4.935	20.00	3.7
                   3.210	10.00	3.7
                   2.695	7.00	3.7
                   2.065	3.00	3.7
                   1.210	2.00	3.7
                   0.659	1.00	3.7
                   0.512	0.70	3.7
                   0.238	0.30	3.7
                   0.188	0.20	3.7
                   0.087	0.10	3.7
                   0.090	0.07	3.7
                   0.026	0.03	3.7
                   0.020	0.02	3.7
                   0.004	0.01	3.7
                   4.480	30.00	7.5
                   4.020	20.00	7.5
                   2.605	10.00	7.5
                   1.935	7.00	7.5
                   1.094	3.00	7.5
                   0.836	2.00	7.5
                   0.467	1.00	7.5
                   0.355	0.70	7.5
                   0.157	0.30	7.5
                   0.123	0.20	7.5
                   0.054	0.10	7.5
                   0.040	0.07	7.5
                   0.021	0.03	7.5
                   0.016	0.02	7.5
                   0.001	0.01	7.5
                   3.530	30.00	15.0
                   2.825	20.00	15.0
                   1.765	10.00	15.0
                   1.460	7.00	15.0
                   0.675	3.00	15.0
                   0.531	2.00	15.0
                   0.275	1.00	15.0
                   0.209	0.70	15.0
                   0.093	0.30	15.0
                   0.079	0.20	15.0
                   0.034	0.10	15.0
                   0.026	0.07	15.0
                   0.013	0.03	15.0
                   0.011	0.02	15.0
                   0.000	0.01	15.0
                   2.335	30.00	30.0
                   2.010	20.00	30.0
                   1.044	10.00	30.0
                   0.970	7.00	30.0
                   0.385	3.00	30.0
                   0.368	2.00	30.0
                   0.151	1.00	30.0
                   0.151	0.70	30.0
                   0.061	0.30	30.0
                   0.048	0.20	30.0
                   0.026	0.10	30.0
                   0.021	0.07	30.0
                   0.007	0.03	30.0
                   0.007	0.02	30.0
                   0.001	0.01	30.0
                   1.175	30.00	60.0
                   1.370	20.00	60.0
                   0.566	10.00	60.0
                   0.713	7.00	60.0
                   0.245	3.00	60.0
                   0.250	2.00	60.0
                   0.104	1.00	60.0
                   0.118	0.70	60.0
                   0.046	0.30	60.0
                   0.038	0.20	60.0
                   0.019	0.10	60.0
                   0.016	0.07	60.0
                   0.006	0.03	60.0
                   0.005	0.02	60.0
                   0.000	0.01	60.0
                   ")

###############################################################################################################
with(basia, scatter.smooth(y ~ x1))
with(basia, scatter.smooth(y ~ x2))

# initial parameter estimation in the standard case 
# for competitive (m1), noncompetitive (m2)  and encompassing (m3) models

m1<-nls(y~b1*x1/(b2*(1+x2/b3)+x1), start= list(b1=10,b2=4.36,b3=2.58), data=basia)
summary(m1)

m2<-nls(y~b1*x1/((b2+x1)*(1+x2/b3)), start= list(b1=10,b2=4.36,b3=2*2.58), data=basia)
summary(m2)

m3<-nls(y~b1*x1/((b2*(1+x2/b3)+x1*(1+(1-b4)*x2/b3))), start= list(b1=10,b2=4.36,b3=1.2*2.58,b4=0.8), data=basia)
summary(m3)

#--------------------------------------------------------------------------------------------------------------
#Scatter plot of residuals versus fitted values, standard case

fmin1<-min(fitted(m3)) #0.0007691968
fmax1<-max(fitted(m3)) #6.423166
rmin1<-min(residuals(m3)) #-0.63836
rmax1<-max(residuals(m3)) #0.4562213


plot(fitted(m3),residuals(m3),xlim=c(0,7),ylim=c(-0.8,0.8),xaxt="n", yaxt="n",xlab="fitted",ylab="residuals")
xtick=seq(0,7,by=1)
ytick=seq(-0.8,0.8,by=.4)
axis(1,at=xtick,labels=T)
axis(2,at=ytick,labels=T)
abline(h=0,lty=2,col="blue")
mtext("(a)",side=3,line=0.4,outer=F,cex=1.3)

#--------------------------------------------------------------------------------------------------------------
# Here we checked different epsilons (small values) to be replaced by the zeros in the data. We tried different
# values and finally the epsilon which makes the least Residual standard error was chosen. 
# This is done for all three non-competitive, competitive and also combined log models. 
# Fortunately we observed that all three log models resulted in the same value which is 0.02.

for (i in seq_len(ncol(basia))) basia[[i]][basia[[i]] == 0] <- 0.001
for (i in seq_len(ncol(basia))) basia[[i]][basia[[i]] == 0.01] <- 0.02

# initial parameter estimation in the log case 
# for competitive (m4), noncompetitive (m5)  and encompassing (m6) models

with(basia, scatter.smooth(log(y) ~ x1))
with(basia, scatter.smooth(log(y) ~ x2))

m4<-nls(log(y)~log(b1*x1/(b2*(1+x2/b3)+x1)), start= list(b1=7.2976,b2=4.3860,b3=2.5821), data=basia)  
summary(m4)

m5<-nls(log(y)~log(b1*x1/((b2+x1)*(1+x2/b3))), start= list(b1=8.6957,b2=8.0664,b3=12.0566), data=basia) 
summary(m5)

m6<-nls(log(y)~log(b1*x1/((b2*(1+x2/b3)+x1*(1+(1-b4)*x2/b3)))), start= list(b1=7.42535,b2=4.68079,b3=3.05814,b4=0.96362), data=basia)
summary(m6)

#--------------------------------------------------------------------------------------------------------------
#Scatter plot of residuals versus fitted values, log case

fmin2<-min(fitted(m6)) #-6.185951
fmax2<-max(fitted(m6)) #1.819805
rmin2<-min(residuals(m6)) #-2.456082
rmax2<-max(residuals(m6)) #0.887634

plot(fitted(m6),residuals(m6),xlim=c(-7,2),ylim=c(-2.5,1),xaxt="n", yaxt="n",xlab="fitted",ylab="residuals")
xtick=seq(-7,2,by=1)
ytick=seq(-2.5,1,by=0.5)
axis(1,at=xtick,labels=T)
axis(2,at=ytick,labels=T)
abline(h=0,lty=2,col="blue")
mtext("(b)",side=3,line=0.4,outer=F,cex=1.3)

#--------------------------------------------------------------------------------------------------------------
#Scatter plot of residuals versus fitted values, back-transformed case

fmin3<-min(exp(fitted(m6))) #0.002058143
fmax3<-max(exp(fitted(m6))) # 6.170654
rmin3<-min(basia$y-exp(fitted(m6))) #-0.6307701
rmax3<-max(basia$y-exp(fitted(m6))) #0.7612291

plot(exp(fitted(m6)),basia$y-exp(fitted(m6)),xlim=c(0,7),ylim=c(-0.8,0.8),xaxt="n", yaxt="n",xlab="fitted",ylab="residuals")
xtick=seq(0,7,by=1)
ytick=seq(-0.8,0.8,by=.4)
axis(1,at=xtick,labels=T)
axis(2,at=ytick,labels=T)
abline(h=0,lty=2,col="blue")
mtext("(c)",side=3,line=0.4,outer=F,cex=1.3)


sqrt(sum((residuals(m3))^2)/116)
sqrt(sum((basia$y-exp(fitted(m6)))^2)/116)
sqrt(sum((residuals(m6))^2)/116)
