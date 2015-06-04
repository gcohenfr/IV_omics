#Illustrate IV with a simulated example
#Plots for talk at Overall's group
#April 20, 2015

library(ggplot2)
library(mvtnorm)

#set.seed(120291)

set.seed(132)
#Generate data with endogeneity
n <-100
beta<-2
Xtrue<-rnorm(n)
X<-Xtrue+rnorm(n,0,0.5)
y<-beta*Xtrue+rnorm(n)
Z<-0.5*X+rnorm(n)

#sample points to highlight
some.points<-sample(n,5)
col.ind<-rep(0,n)
col.ind[some.points]<-1
size.ind<-rep(1,n)
size.ind[some.points]<-1.1

#Arrage data
dat.example<-data.frame(y=y,p=X,p.true=Xtrue,g=Z,size.ind=size.ind,col.ind=factor(col.ind))

subset.example<-dat.example[some.points,]
x.arrow = subset.example[5,"p.true"]
y.arrow = subset.example[5,"y"]
xend.arrow = subset.example[5,"p"]
yend.arrow = subset.example[5,"y"]


#Plot 1: true data, true line
ggplot(dat.example,aes(p.true,y))+geom_point(colour="grey")+
  geom_abline(intercept = 0, slope = beta,colour="blue",size=1.5)+xlab("protein")+ylab("stenosis")

#Plot 2: observed data, true line, estimated line
ggplot(dat.example,aes(p.true,y))+
  geom_point(colour="grey")+
  geom_point(aes(p,y))+
  geom_abline(intercept = 0, slope = beta,colour="blue",size=1.5)+
  stat_smooth(aes(p,y),method="lm", se=FALSE,colour="red",size=1.5)+
  geom_text(x=1.8,y=4.8,label="true",colour="blue")+
  geom_text(x=2.5,y=2.8,label="LS",colour="red")+
  xlab("protein")+ylab("stenosis")

#Plot 3: highlight some points
ggplot(dat.example,aes(p.true,y))+
  geom_point(colour="grey")+
  geom_point(aes(p,y))+
  geom_abline(intercept = 0, slope = beta,colour="blue",size=1.5)+
  stat_smooth(aes(p,y),method="lm", se=FALSE,colour="red",size=1.5)+
  geom_point(data=subset.example,colour="blue",size=3)+
  geom_point(data=subset.example,aes(p,y),colour="red",size=3)+
  geom_segment(aes(x=x.arrow,y=y.arrow,xend=xend.arrow,yend=yend.arrow),colour="green")+
  geom_text(x=1.8,y=4.8,label="true",colour="blue")+
  geom_text(x=2.5,y=2.8,label="LS",colour="red")+
  xlab("protein")+ylab("stenosis")

# First stage prediction

p.1st<-predict(lm(p~g,dat.example))
dat.example<-transform(dat.example,p1st=p.1st)

#Plot 4: first stage
ggplot(dat.example,aes(g,p))+
  geom_point()+
  stat_smooth(method="lm", se=FALSE,colour="orange",size=1.5)+
  xlab("instrumental variable")+ylab("protein")

#Plot 5: protein prediction
ggplot(dat.example,aes(g,p))+
  geom_point()+
  geom_point(aes(g,p1st),colour="orange",size=3)+
  stat_smooth(method="lm", se=FALSE,colour="orange")+
  xlab("instrumental variable")+ylab("protein")  

#Plot 6: highlight some points
subset.example<-dat.example[some.points,]
x.arrow = subset.example[5,"g"]
y.arrow = subset.example[5,"p"]
xend.arrow = subset.example[5,"g"]
yend.arrow = subset.example[5,"p1st"]

ggplot(dat.example,aes(g,p))+
  geom_point()+
  geom_point(aes(g,p1st),colour="orange",size=3)+
  geom_point(data=subset.example,aes(g,p),colour="red",size=3)+
  stat_smooth(method="lm", se=FALSE,colour="orange")+
  geom_segment(aes(x=x.arrow,y=y.arrow,xend=xend.arrow,yend=yend.arrow),colour="green")+
  xlab("instrumental variable")+ylab("protein")  

#Plot 7: second stage
coef.2st=coef(lm(y~p1st,data=dat.example))

ggplot(dat.example,aes(p,y))+geom_point()+
  geom_point(aes(p1st,y),colour="blue")+
  geom_abline(intercept = 0, slope = beta,colour="blue",size=1.3)+
  geom_abline(intercept = coef.2st[1], slope = coef.2st[2],colour="green",size=1.3)+
  stat_smooth(aes(p,y),method="lm", se=FALSE,colour="red",size=1.5)+
  geom_text(x=1.8,y=4.8,label="true",colour="blue")+
  geom_text(x=2.5,y=2.8,label="LS",colour="red")+
  geom_text(x=3.1,y=4.8,label="IV",colour="green")+
  xlab("protein")+ylab("stenosis")

#Plot 8: more data
#sum.iv<-c()
#for(i in 1:1000){
n <-1000
beta<-2
Xtrue<-rnorm(n)
X<-Xtrue+rnorm(n,0,0.5)
eps<-rnorm(n)
y<-beta*Xtrue+eps
Z<-0.5*Xtrue+rnorm(n)

#Arrage data
dat.big<-data.frame(y=y,p=X,p.true=Xtrue,g=Z)
#iv<-ivreg(y~p|g,data=dat.big)
#lm.n<-lm(y~p,data=dat.big)

#sum.iv<-rbind(sum.iv,c(coef(iv)[2],coef(lm.n)[2]))}

#2SLS
p.1st<-predict(lm(p~g,dat.big))
dat.big<-transform(dat.big,p1st=p.1st)

coef.2st=coef(lm(y~p1st,data=dat.big))

ggplot(dat.big,aes(p,y))+geom_point()+
  geom_abline(intercept = 0, slope = beta,colour="blue",size=1.3)+
  geom_abline(intercept = coef.2st[1], slope = coef.2st[2],colour="green",size=1.3)+
  stat_smooth(aes(p,y),method="lm", se=FALSE,colour="red",size=1.5)+ 
  geom_text(x=2.8,y=6.8,label="true",colour="blue")+
  geom_text(x=4.5,y=6.8,label="LS",colour="red")+
  geom_text(x=4,y=7,label="IV",colour="green")+
  xlab("protein")+ylab("stenosis")



