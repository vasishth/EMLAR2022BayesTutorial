# Analytical example Bayes factor

plik1 <- function(theta) {
  dbinom(x = 0, size = 10, prob = theta) *
    dbeta(x = theta, shape1 = 1, shape2 = 1) }
# Then we integrate (compute the area under the curve):
(MargLik1 <- integrate(f = plik1, lower = 0, upper = 1)$value)
library(MASS)
fractions(MargLik1)

## load example data-set:
gw<-read.table("data/gibsonwu2012data.txt",
               header=TRUE)
## sum-contrast coding of predictor:
gw$so <- ifelse(
  gw$type%in%c("subj-ext"),-1,1)
## subset critical region
gw1<-subset(gw,region=="headnoun")

## load second data-set:
gw2<-read.table("data/gibsonwu2012datarepeat.txt",
                header=TRUE)
gw2$so <- ifelse(
  gw2$condition%in%c("subj-ext"),-1,1)

## frequentist analysis:
library(lme4)
m_lmer<-lmer(rt~so + (1+so|subj)+(1+so|item),gw1)
summary(m_lmer)

boxplot(rt~so,gw1)

m_lmerlog<-lmer(log(rt)~so + (1+so|subj)+(1+so|item),gw1)
summary(m_lmerlog)

boxplot(log(rt)~so,gw1)


library(brms)
priors <- c(set_prior("normal(6, 1.5)", class = "Intercept"),
            set_prior("normal(0, .02)", class = "b", 
                      coef = "so"),
            set_prior("normal(0, 1)", class = "sd"),
            set_prior("normal(0, 1)", class = "sigma"),
            set_prior("lkj(2)", class = "cor"))

m_gw<-brm(rt~so + (1+so|subj) + (1+so|item),gw1,family=lognormal(),
          prior=priors)
summary(m_gw)

pp_check(m_gw)

## graphical visualization:
library(bayesplot)
postgw<-posterior_samples(m_gw)
## extract posteriors:
alpha<-postgw$b_Intercept
beta<-postgw$b_so
cor<-posterior_samples(m_gw,"^cor")
sd<-posterior_samples(m_gw,"^sd")
sigma<-posterior_samples(m_gw,"sigma")

## subject level effects:
subj_re<-posterior_samples(m_gw,"^r_subj")
item_re<-posterior_samples(m_gw,"^r_item")


## mean effect on ms scale:
meandiff<- exp(alpha + beta) - exp(alpha - beta)
mean(meandiff)
round(quantile(meandiff,prob=c(0.025,0.975)),0)

## mean effect:
hist(meandiff,freq=FALSE,
     main="Mean OR vs SR processing cost",
     xlab=expression(exp(alpha + beta)- exp(alpha - beta)))

## individual level estimates:
nsubj<-37
subjdiff<-matrix(rep(NA,nsubj*4000),nrow=nsubj)
for(i in 1:nsubj){
  subjdiff[i,]<-exp(alpha + subj_re[,i]  + (beta+subj_re[,i+nsubj])) - 
    exp(alpha + subj_re[,i] - 
          (beta+subj_re[,i+nsubj]))
}

subjdiff<-t(subjdiff)

subjdiff<-as.data.frame(subjdiff)
colnames(subjdiff)<-paste("s",c(1:nsubj),sep="")
mns <- colMeans(subjdiff)
subjdiff<-subjdiff[,order(mns)]
mcmc_areas(subjdiff)

## lmer
m_lmer<-lmer(log(rt)~so + (1+so|subj),gw1)
summary(m_lmer)


alphalmer<-coef(m_lmer)$subj[,1]
betalmer<-coef(m_lmer)$subj[,2]
lmerest<-exp(alphalmer+betalmer)-exp(alphalmer-betalmer)


## lmList:
m_lmlist<-lmList(log(rt)~so|subj,gw1)
summary(m_lmlist)

alphalmlist<-coef(m_lmlist)[1]
betalmlist<-coef(m_lmlist)[2]
lmlistest<-exp(alphalmlist+betalmlist)-exp(alphalmlist-betalmlist)

comparison<-data.frame(subj=factor(1:37),brm=mns,lmlist=lmlistest,lmer=lmerest)
colnames(comparison)[c(3,4)]<-c("lmlist","lmer")

plot(lmlist~subj,comparison,ylab="estimate",main="lmlist vs Bayes\n shrinkage in action")
points(brm~subj,pch=3,comparison)
arrows(x0=1:37,y0=comparison$lmlist,x1=1:37,y1=comparison$brm,
       angle=45,code=2,length=0.1,col="red")

plot(lmlist~subj,comparison,ylab="estimate",main="lmlist vs lmer\n shrinkage in action")
points(lmer~subj,pch=3,comparison)
arrows(x0=1:37,y0=comparison$lmlist,x1=1:37,y1=comparison$lmer,
       angle=45,code=2,length=0.1,col="red")

plot(lmer~subj,comparison,ylab="estimate",main="lmer vs Bayes\n regularization in action")
points(brm~subj,pch=3,comparison)
arrows(x0=1:37,y0=comparison$lmer,x1=1:37,y1=comparison$brm,
       angle=45,code=2,length=0.1,col="red")



## Bayes factors
## Bayes factor analysis:
m_gw<-brm(rt~so + (1+so|subj) + (1+so|item),gw1,family=lognormal(),
          prior=priors,warmup=5000,iter=20000,
          ssave_pars=save_pars(all=TRUE))
summary(m_gw)

priors0 <- c(set_prior("normal(6, 1.5)", class = "Intercept"),
             set_prior("normal(0, 1)", class = "sd"),
             set_prior("normal(0, 1)", class = "sigma"),
             set_prior("lkj(2)", class = "cor"))

m_gw0<-brm(rt~1 + (1+so|subj),gw1,family=lognormal(),
           prior=priors0,,warmup=5000,iter=20000,
           save_pars=save_pars(all=TRUE))

bayes_factor(m_gw,m_gw0)
