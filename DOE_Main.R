####ANOVA table
score=c(10,15,17,19,10,17,19,21,15,14,21,14,10,15,12,21,27,17,18,19,19,3,21,19,
        0,0,0,0,0,7,3,7,7,0,0,8,0,12,3,0,17,5,13,3,0,5,0,0,
        21,19,17,8,19,21,19,15,17,11,21,10,14,15,21,21,21,19,17,14,12,19,21,17,
        10,11,13,17,17,17,17,21,17,17,12,3,15,15,15,19,12,14,19,10,8,19,19,8,
        3,3,17,0,10,7,7,17,7,3,14,0,10,12,3,10,5,5,7,6,17,6,0,7)
distance=c(10,8,8,8,12,12,8,12,12,10,8,8,12,10,10,8,8,10,12,10,10,12,12,10,
           12,12,10,10,8,8,10,8,8,10,12,8,12,8,12,12,8,10,8,10,10,12,12,10,
           8,10,10,12,10,10,8,12,12,12,8,12,12,10,8,8,8,12,8,10,12,10,10,8,
           12,8,10,12,8,8,10,8,10,8,12,10,12,10,12,10,12,10,8,12,12,10,8,8,
           10,12,8,12,12,12,10,8,8,12,12,12,10,10,8,8,8,8,10,10,8,10,12,10)
posture=c('K','S','K','K','S','S','S','K','S','S','S','K','S','K','K','S','K',
          'K','K','S','S','K','K','S','K','S','S','S','K','K','S','K','S','K',
          'S','S','K','K','K','S','S','K','S','S','K','S','K','K','S','S','K',
          'S','S','K','K','K','K','K','S','K','S','S','S','K','S','S','K','K',
          'S','S','K','K','S','S','K','S',"S",'K','K','K','K','K','K','K','K',
          'S','K','S','K','S','S','S','S','S','K','S','S','S','S','K','S','S',
          'K','K','S','S','K','K','K','S','K','K','S','K','S','S','S','K','K','K')
horizontal=c('H','R','H','R','R','H','H','R','H','R','R','H','R','R','R','H','R',
             'H','H','H','H','R','H','R','R','H','H','R','R','H','H','R','R','R',
             'R','H','R','H','H','H','H','H','R','R','R','R','H','H','R','R','H',
             'H','H','H','R','H','H','R','H','R','R','H','R','H','H','H','R','R',
             'R','R','R','H','H','R','H','H','H','H','R','R','R','R','R','H','H',
             'H','H','R','R','H','H','R','R','R','H','R','H','R','H','H','H','R',
             'R','H','R','H','H','R','H','R','H','R','R','R','H','R','H','H','R','R')
hand=c('ND','ND','D','D','D','ND','ND','D','D','ND','D','ND','ND','D','ND','D','ND',
       'D','D','ND','D','ND','ND','D','D','D','D','D','D','D','ND','ND','D','ND','D',
       'D','ND','ND','D','ND','ND','D','ND','ND','D','ND','ND','ND','D','D','ND','ND',
       'D','D','D','D','ND','ND','D','D','D','ND','ND','D','ND','D','ND','ND','ND',
       'ND','D','ND','ND','ND','ND','D','D','ND','ND','D','D','ND','D','D','D','D',
       'ND','ND','ND','ND','ND','ND','D','D','D','D','D','ND','ND','ND','D','D','ND',
       'ND','ND','ND','D','ND','ND','ND','D','ND','D','D','ND','D','D','D','D','D')
person=c(rep(1,24),rep(2,24),rep(3,24),rep(4,24),rep(5,24))
dart=data.frame(score=score,A=factor(distance),B=factor(horizontal),C=factor(posture),D=factor(hand),blocks=factor(person))
attach(dart)
interaction.plot(A,B,score)
interaction.plot(A,C,score)
interaction.plot(A,D,score)
interaction.plot(Hard,Press,score)
interaction.plot(Hard,Cook,Strength)
interaction.plot(Cook,Press,Strength)
outblock=lm(score~blocks+A*B*C*D,dart)
summary.aov(outblock)



#####Assumptions checking
par(mfrow=c(1,2))
qqnorm(outblock$residuals,main="NP Plot of Residuals")
qqline(outblock$residuals)
plot(outblock$fitted.values,outblock$residuals,main="Residuals vs. Fits")
plot(outblock$residuals,ylab="Residuals",type="l",main="Residuals vs. Index (Presumably Time)")
par(mfrow=c(1,1))


#####Box-Cox transformations
###user must specify values to calculate dfe in the lines below for the model considered ###

dfe=92
inc=.1	#use lambdas in increments of this value
lamlist=seq(-1,1,inc)
sse=rep(NA,length(lamlist))
dart$score2<-dart$score+0.01
ydot=prod(dart$score2)^(1/length(dart$score2))

for(i in 1:length(lamlist)) {
  lam=lamlist[i]
  ylam=(dart$score2^lam-1 )/(lam*(ydot^(lam-1)))
  if(lam==0) {ylam=ydot*log(dart$score2)}
  dart2=data.frame(ylam,A=factor(distance),B=factor(posture),C=factor(horizontal),D=factor(hand),blocks=factor(person))
  attach(dart2)
  sse[i]=sum(lm(ylam~blocks+A*B*C*D,dart2)$residuals^2)
  detach(dart2)
}

ssstar=min(sse)*(1+(qt(.025,dfe)^2)/dfe)
ssstar
plot(lamlist,sse,ylim=c(min(sse),max(ssstar,max(sse))))
title("Sse vs. lambda for Box-Cox: horizontal line = ss*")
lines(c(min(lamlist),max(lamlist)),c(ssstar,ssstar))
opt_lam=lamlist[sse==min(sse)]
sse
ssstar
min(sse)
opt_lam


library(MASS)
aovout=aov(score2~blocks+A*B*C*D,dart)
boxcox(aovout)

##### New model
dart$score2<-dart$score+0.01
dart$score3<-dart$score2^0.7
dart=data.frame(score3,A=factor(distance),B=factor(posture),C=factor(horizontal),D=factor(hand),blocks=factor(person))
attach(dart)
outblock3=lm(score3~blocks+A*B*C*D,dart)
summary.aov(outblock3)


##### New Assumptions checking
par(mfrow=c(1,2))
qqnorm(outblock3$residuals,main="NP Plot of Residuals")
qqline(outblock3$residuals)
plot(outblock3$fitted.values,outblock3$residuals,main="Residuals vs. Fits")
plot(outblock3$residuals,ylab="Residuals",type="l",main="Residuals vs. Index (Presumably Time)")
par(mfrow=c(1,1))

install.packages('car')
library('car')
ncvTest(outblock3)




