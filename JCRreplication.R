########## Replication file ########
## "Enemies Within"
## Casey Crisman-Cox
## Texas A&M
## c.crisman-cox@tamu.edu
## Journal of Conflict Resolution
## UPDATED: 15 May 2018
##     - Changed contact info


#######Libraries, Functions, and Data#####
library(lmtest) 
library(MASS)
library(VGAM)
library(multiwayvcov)
library(brglm)
library(ggplot2) 
library(gridExtra)
library(lme4)
library(maxLik)
library(brms) 
library(matrixStats)
library(stringr)



rm(list=ls())

load("Terror_JCR_Replication.Rdata") 
source("functions/extraFunctions_JCR.R")



####First reg table ####


dyadMod1 <- glm.nb(attacks~left
                   +center
                   +tpop
                   +milper
                   +rgdpl
                   +polity2
                   +PR
                   +govfrac
                   +lag.attacks
                   +press.use
                   +growth_gdp
                   +ethfrac
                   +relfrac
                   +coldwar,
                   data=data.full,
                   maxit=1000,
                   x=TRUE, 
                   y=TRUE)

dyadMod2 <- glm.nb(attacks~left
                   +center
                   +tpop
                   +milper
                   +rgdpl
                   +polity2
                   +PR
                   +govfrac
                   +lag.attacks
                   +press.use
                   +growth_gdp
                   +ethfrac
                   +relfrac
                   +coldwar
                   +I(Size>2)
                   +Incomp
                   +intraD
                   +duration
                   +I(duration^2)
                   , data=data.full,
                   maxit=1000,
                   x=TRUE,
                   y=TRUE)
##Convergence issues if we use Size <=2 , so try it in maxLik
dyadMod2a <- maxLik(NBreg, NBgrad,
                    start=c(dyadMod2$coef, 1/dyadMod2$theta), 
                    X=dyadMod2$x, Y=dyadMod2$y,
                    method="NR", control=list(tol=1e-20, gradtol=1e-20))
#MaxLik converged so let's use it
dyadMod2$coefficients <- dyadMod2a$est[-length(dyadMod2a$est)]
dyadMod2$theta <- 1/exp(dyadMod2a$est[length(dyadMod2a$est)])




dyadMod3 <- glm.nb(attacks~right
                   +tpop
                   +milper
                   +rgdpl
                   +polity2
                   +PR
                   +govfrac
                   +lag.attacks
                   +press.use
                   +growth_gdp
                   +ethfrac
                   +relfrac
                   +coldwar
                   +I(Size>2)
                   +Incomp
                   +intraD
                   +duration
                   +I(duration^2)
                   , data=data.full,
                   maxit=1000,
                   x=TRUE,
                   y=TRUE)
##Convergence issues if we use Size <=2 , so try it in maxLik
dyadMod3a <- maxLik(NBreg, NBgrad,
                    start=c(dyadMod3$coef, 1/dyadMod3$theta), 
                    X=dyadMod3$x, Y=dyadMod3$y,
                    method="NR", control=list(tol=1e-20, gradtol=1e-20))

#MaxLik converged so let's use it
dyadMod3$coefficients <- dyadMod3a$est[-length(dyadMod3a$est)]
dyadMod3$theta <- 1/exp(dyadMod3a$est[length(dyadMod3a$est)])
########################################################


vcov2.1 <- ClusterNegbin(dyadMod1, "iyear", data.full, vcov(dyadMod1))
print(coeftest(dyadMod1, vcov2.1$vcov))

vcov2.2 <- ClusterNegbin(dyadMod2, "iyear", data.full, vcov(dyadMod2a)[-length(dyadMod2a$est),
                                                                       -length(dyadMod2a$est)])
print(coeftest(dyadMod2, vcov2.2$vcov))

vcov2.3 <- ClusterNegbin(dyadMod3, "iyear", data.full, vcov(dyadMod3a)[-length(dyadMod3a$est),
                                                                       -length(dyadMod3a$est)])
print(coeftest(dyadMod3, vcov2.3$vcov))



## Figure 1


X <- dyadMod1$model
X[,1] <- 1
X <- as.data.frame(av.data(X, length.out=3))
X$left <- c(0,0,1)
X$center <- c(0,1,0)
BetaPar <- mvrnorm(1000, dyadMod1$coef, vcov2.1$vcov)
XB <- exp(as.matrix(X)%*% t(BetaPar))
fa <- data.frame(hi = apply(XB, 1, quantile, 0.975),
                 lo = apply(XB, 1, quantile, 0.025),
                 fit= apply(XB, 1, quantile, 0.5),
                 Party=factor(c("Right", "Center", "Left"),
                              levels=c("Right", "Center", "Left"))
)

fda <-rbind(XB[3,]-XB[1,],
            XB[3,]-XB[2,],
            XB[1,]-XB[2,])

fda <- data.frame(hi = apply(fda, 1, quantile, 0.975),
                  lo = apply(fda, 1, quantile, 0.025),
                  fit= apply(fda, 1, quantile, 0.5),
                  From=c("Right", "Center", "Center"),
                  To=c("Left", "Left", "Right"),
                  FromTo=factor(c("Right->Left", "Center->Left", "Center->Right"),
                                levels=c("Right->Left", "Center->Left", "Center->Right"))
)

grid.arrange(
     ggplot(fa)+
          geom_pointrange(aes(x=Party, y=fit, ymin=lo, ymax=hi), size=1)+
          theme_bw(22)+
          ylab("Number of Attacks")+
          ggtitle("Estimated Number of Attacks"),
     ggplot(fda)+
          geom_pointrange(aes(x=FromTo, y=fit, ymin=lo, ymax=hi), size=1)+
          theme_bw(22)+ 
          geom_hline(yintercept = 0, alpha=.3)+
          ylab("Difference")+
          xlab("Change in Party")+
          scale_x_discrete(labels=expression(Right%->%Left, 
                                             Center%->%Left, 
                                             Center%->%Right))+
          ggtitle("Estimated Change in Attacks"),          
     nrow=2
)


## Table 2 Models
interMod1 <- glm.nb(attacks~left
                    +center
                    +tpop
                    +milper
                    +rgdpl
                    +polity2
                    +PR
                    +govfrac  
                    +lag.attacks
                    +press.use
                    +growth_gdp
                    +ethfrac
                    +relfrac
                    +coldwar
                    +I(Size>2)
                    +Incomp
                    +intraD
                    +duration
                    +I(duration^2)
                    , data=data.full, subset=Left==0,
                    maxit=1000,
                    x=TRUE, 
                    y=TRUE)
##Convergence issues, so try it in maxLik
interMod1a <- maxLik(NBreg, NBgrad,
                     start=c(interMod1$coef, 1/interMod1$theta), 
                     X=interMod1$x, Y=interMod1$y,
                     control=list(tol=1e-20, gradtol=1e-20))
#MaxLik converged so let's use it
interMod1$coefficients <- interMod1a$est[-length(interMod1a$est)]
interMod1$theta <- 1/exp(interMod1a$est[length(interMod1a$est)])


interMod2 <- glm.nb(attacks~left*Left
                    +center
                    +tpop
                    +milper
                    +rgdpl
                    +polity2
                    +PR
                    +govfrac  
                    +lag.attacks
                    +press.use
                    +growth_gdp
                    +ethfrac
                    +relfrac
                    +coldwar
                    +I(Size>2)
                    +Incomp
                    +intraD
                    +duration
                    +I(duration^2)
                    , data=data.full,
                    maxit=500,
                    x=TRUE, 
                    y=TRUE)
##Convergence issues, so try it in maxLik
interMod2a <- maxLik(NBreg, NBgrad,
                     start=c(interMod2$coef, 1/interMod2$theta), 
                     X=interMod2$x, Y=interMod2$y,
                     control=list(tol=1e-20, gradtol=1e-20))
#MaxLik converged so let's use it
interMod2$coefficients <- interMod2a$est[-length(interMod2a$est)]
interMod2$theta <- 1/exp(interMod2a$est[length(interMod2a$est)])

######Size interactions models######     


interSizeMod1 <- glm.nb(attacks~left*I(Size>2)
                        +center
                        +tpop
                        +milper
                        +rgdpl
                        +polity2
                        +PR
                        +govfrac  
                        +lag.attacks
                        +press.use
                        +growth_gdp
                        +ethfrac
                        +relfrac
                        +coldwar
                        +Incomp
                        +intraD
                        +duration
                        +I(duration^2)
                        , data=data.full,
                        maxit=5000,
                        x=TRUE, 
                        y=TRUE)
# ##Convergence issues, so try it in maxLik
interSizeMod1a <- maxLik(NBreg, NBgrad,
                         start=c(interSizeMod1$coef, 1/interSizeMod1$theta),
                         X=interSizeMod1$x, Y=interSizeMod1$y,
                         control=list(tol=1e-20, gradtol=1e-20))
#MaxLik converged so let's use it
interSizeMod1$coefficients <- interSizeMod1a$est[-length(interSizeMod1a$est)]
interSizeMod1$theta <- 1/exp(interSizeMod1a$est[length(interSizeMod1a$est)])


vcov3.1 <- ClusterNegbin(interMod1, "iyear", data.full, vcov(interMod1))
print(coeftest(interMod1, vcov3.1$vcov))

vcov3.2 <- ClusterNegbin(interMod2, "iyear", data.full, vcov(interMod2a)[-length(interMod2a$est),-length(interMod2a$est)])
print(coeftest(interMod2, vcov3.2$vcov))

vcov3.3 <- ClusterNegbin(interSizeMod1, "iyear", data.full, vcov(interSizeMod1a)[-length(interSizeMod1a$est),
                                                                                 -length(interSizeMod1a$est)])
print(coeftest(interSizeMod1, vcov3.3$vcov))







## Table 3 Models
load("FullDataPlusNELDA.rdata")
m1 <- vglm(mnl~lag.attacks
           +tpop
           +rgdpl
           +govfrac  
           +intraD
           +press.use
           +growth_gdp
           +coldwar
           +duration
           +I(duration^2)
           , x=T, y=T, model=T,
           data=data.full, family=multinomial(refLevel = 1))
m2 <- vglm(mnl~lag.attacks
           +tpop
           +rgdpl
           +govfrac  
           +intraD
           +press.use
           +growth_gdp
           +coldwar
           +duration
           +I(duration^2)
           ,
           data=data.full, family=multinomial(refLevel = 1), subset=electionlastyear==1)



m3 <- vglm(mnl~lag.attacks + lag.attacks2
           +tpop
           +rgdpl
           +govfrac  
           +intraD
           +press.use
           +growth_gdp
           +coldwar
           +duration
           +I(duration^2)
           , x=T, y=T, model=T,
           data=data.full, family=multinomial(refLevel = 1))
m4 <- vglm(mnl~lag.attacks + lag.attacks2
           +tpop
           +rgdpl
           +govfrac  
           +intraD
           +press.use
           +growth_gdp
           +coldwar
           +duration
           +I(duration^2)
           , 
           data=data.full, family=multinomial(refLevel = 1), subset=electionlastyear==1)

print(summary(m1))
print(summary(m2))
print(summary(m3))
print(summary(m4))

## Table 4 Models
step1build1 <- brglm(
     agreement~left
     +center
     +log(cumAtt+1)
     +rgdpl
     +PR
     +tpop
     +intraD
     +milper
     +Incomp
     +govfrac
     +polity2
     +press.use
     +I(Size>2)
     +I(duration/10)
     , 
     data=data.full, 
     family=binomial,
     x=TRUE, y=TRUE, model=TRUE)

step1build2 <- brglm(
     agreement~left*Left
     +center
     +log(cumAtt+1)
     +rgdpl
     +PR
     +tpop
     +intraD
     +milper
     +Incomp
     +govfrac
     +polity2
     +press.use
     +I(Size>2)
     +I(duration/10)
     , 
     data=data.full, 
     family=binomial,
     x=TRUE, y=TRUE, model=TRUE)



v1 <- ClusterLogit(step1build1, "iyear", data.full)
v2 <- ClusterLogit(step1build2, "iyear", data.full)

print(coeftest(step1build1, v1$vcov))
print(coeftest(step1build2, v2$vcov))

#seed argument doesn't work, results will vary slightly from the printed version
step1GLMM <- brm(
     agreement~left
     +center
     +log(cumAtt+1)
     +rgdpl
     +PR
     +tpop
     +intraD
     +milper
     +govfrac
     +polity2
     +Incomp
     +press.use
     +I(Size>2)
     +I(duration/10)
     +(1|ccode:DyadID)
     +(1|ccode)
     , 
     data=data.full, 
     core=4, iter=4000, seed=12345, #seed doesn't work, results vary, but are close
     prior=set_prior("cauchy(0,2.5)", class='b'),
     family=bernoulli)


glmm.coef <- colMeans(posterior_samples(step1GLMM)[-(1:4000), 1:15])
glmm.se <- colSds(as.matrix(posterior_samples(step1GLMM)[-(1:4000), 1:15]))
names(glmm.coef) <- names(glmm.se) <- names(step1build1$coef)
print(cbind(glmm.coef, glmm.se))

## Table 5 Models

load("Israel.RData") #dugan and chenowith

filler <- c(161/6570,
            167/6690,
            167/6809)


Israel <- x
Israel$date <- str_c(x$year, x$month, sep="-")
Israel$date <- as.yearmon(Israel$date, "%Y-%m")


Idat<- subset(data.full, ccode==666, select=c('left', 
                                              'iyear', 
                                              'govfrac',
                                              'rgdpl',
                                              'coldwar',
                                              'milper',
                                              'growth_gdp',
                                              'gini'))

Idat <- unique(Idat)
Idat$milper[Idat$iyear>=2002 & Idat$iyear<=2004] <- log(filler+.1)
Israel <- merge(Israel, Idat, by.x="year", by.y="iyear")
Israel$lagAttack <- lag(as.zoo(Israel$att93miss), -1, na.pad=TRUE)

Israel$left[Israel$date< "July 1992"] <- 0
Israel$left[Israel$date>= "July 1992" & Israel$date<= "June 1996"] <- 1
Israel$left[Israel$date>  "June 1996" & Israel$date< "July 1999"] <- 0
Israel$left[Israel$date>= "July 1999" & Israel$date< "March 2001"] <- 1
Israel$left[Israel$date>= "March 2001"] <- 0

nb1 <- glm.nb(concil~left+lagAttack+log(rgdpl)+coldwar+govfrac, data=Israel, x=TRUE)
nb2 <- glm.nb(repress~left+lagAttack+log(rgdpl)+coldwar+govfrac, data=Israel, x=TRUE)
print(summary(nb1))
print(summary(nb2))
