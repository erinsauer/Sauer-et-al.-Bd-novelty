####################################################################
######## This code was developed by Erin L. Sauer     ##############
######## for: Sauer et al. "Are novel or locally      ##############
######## adapted pathogens more devastating and why?: ##############
######## Resolving opposing hypotheses"               ##############
####################################################################

###### Figure 1 ##########
library(tidyverse)
# make hypothetical databas
set.seed(125)
a <- rbeta(25, 2,5)
hist(a)
s <- as.numeric(rbeta(25,4,5))
hist(s)
id <- 1:25
al <- c("novel")
sy <- c("local")
x1 <-  c(1)
x2 <- c(2)
db1 <- cbind(a, id, al, x2)
db2 <- cbind(s, id, sy, x1)
db <- data.frame(rbind(db1, db2))
str(db)
View(db)
db$x2 <- as.numeric(db$x2)
db$id <- as.character(db$id)
db$a <- as.numeric(db$a)

# Kaltz & Skykoff 1998 style figure ---------------------------------------
#Figure 1B
base <- ggplot(db, aes(x = x2, y = a, colors=id)) +
  geom_point(colour="grey", size=5) +
  geom_line(colour="grey", size=2) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=25,face="bold")) +
  xlab("Local host                          Novel host") +
  ylab("Pathogen performance") 
base
summary(a) #0.31699
summary(s) #0.4875

pos <- data.frame(cbind(c(0.37246255, 0.59860060), c(50,50), c(1,2)))
colnames(pos) <- c("a", "id","x2")
posline <- base+
  geom_point(data = pos, colour="purple", size=5) +
  geom_line(data = pos, colour="purple", size=2) 

vir <- data.frame(cbind(c(0.58483373, 0.56933427), c(60,60), c(1,2)))
colnames(vir) <- c("a", "id","x2")
virline <- posline+
  geom_point(data = vir, colour="orange", size=5) +
  geom_line(data = vir, colour="orange", size=2) 

summary(a);summary(s)
mean <- data.frame(cbind(c(0.4875, 0.31699), c(57,57), c(1,2)))
colnames(mean) <- c("a", "id","x2")

meanline <- virline+
  geom_point(data = mean, colour="black", size =5) +
  geom_line(data = mean, colour="black", linetype = "dashed", size=2) 
meanline  

# Curve comparision version -----------------------------------------------
#Figure 1A

set.seed(125)
a2 <- rbeta(10000, 2,5)
s2 <- rbeta(10000, 4,5)
al2 <- c("zallopatric")
sy2 <- c("sympatric")
adb <- cbind(a2, al2)
sdb <- cbind(s2, sy2)
db.2 <- data.frame(rbind(adb, sdb))
db.2$a2 <- as.numeric(db.2$a2)
str(db.2)

fig1a <- ggplot(db.2, aes(a2, color = al2)) +
  geom_density(kernel = "gaussian", size =2, adjust=2) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=25,face="bold"),
        legend.position = "none") +
  annotate(geom="text", x=.1, y=2.5, label="Novel interactions",
           color="#01BFC4", size=7.5)+
  annotate(geom="text", x=0.7, y=2.4, label="Local interactions",
           color="#F8766D", size=7.5)+
  xlim(-.2,1)+ 
  labs(y="Density", x="Pathogen performance")
fig1a
library(ggpubr)
fig1 <- ggarrange(fig1a,meanline, ncol=2, labels=c("A", "B"),
                  font.label = list(size = 25))
fig1

###### analysis of the experimental data ##########
library(survival)
library(coxme)
library(car)
library(tidyverse)

surv70 <- read.csv("5x6 experimental data.csv")

#survival model
surv <-coxph(Surv(Date, mortality) ~   Host + Bd +mass+ log(Dist2+1), 
             data = surv70, na.action="na.fail")
summary(surv)
Anova(surv)

#Figure S2a&b
figS2a <- emmeans(surv, list(pairwise ~ Host), adjust = "tukey")
str(summary(figS2a$`emmeans of Host`))
figS2a <- summary(figS2a$`emmeans of Host`)[1:6]
str(figS2a)
log2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
figS2a$emmean <- log2prob(figS2a$emmean)
figS2a$asymp.LCL <- log2prob(figS2a$asymp.LCL)
figS2a$asymp.UCL <- log2prob(figS2a$asymp.UCL)
plot(emmean~Host, data=figS2a)
figS2a$tukey <- c("A","B", "A", "B", "A")
figS2a

figS2b <- emmeans(surv, list(pairwise ~ Bd), adjust = "tukey")
str(summary(figS2b$`emmeans of Bd`))
figS2b <- summary(figS2b$`emmeans of Bd`)[1:6]
figS2b$emmean <- log2prob(figS2b$emmean)
figS2b$asymp.LCL <- log2prob(figS2b$asymp.LCL)
figS2b$asymp.UCL <- log2prob(figS2b$asymp.UCL)
plot(emmean~Bd, data=figS2b)
figS2b$tukey <- c("A","B", "C", "A,C", "A", "A")
figS2b

#Figure 2A

points <- predict(surv, surv70)
View(points)
unname(points)
pps <- as.data.frame(cbind(points, surv70$Dist))

plot(points~V2, data=pps)
fig2a <- ggplot(pps, aes(x=log10(V2+1), y=points))+
  geom_smooth(method=glm, color="black")+
  geom_point()+
  theme(axis.text=element_text(size=30),
        panel.background=element_blank(),
        panel.grid.major.y=element_line(color="grey", size=.2),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=40),
        legend.position="none")+
  ylab("Relative mortality risk (log odds)") +
  xlab("Distance (log10(km)) between \nhost & Bd collection sites")
fig2a

###### zoospore load analysis #########
zglm <- glm(logZoospore ~ Host + Bd + mass + log(Dist2+1) + Date14, data=surv70, 
            na.action="na.fail")
summary(zglm)
Anova(zglm)


###### Bd prevalence analysis ########
surv70$inf.status <- ifelse(surv70$Zoospore>0, 1, ifelse(surv70$Zoospore==0, 0, NA)) 
prvl <- glm(inf.status ~ Host + Bd  + log(Dist2+1) + mass, data=surv70, family=binomial(link="logit"))
summary(prvl)
Anova(prvl)

#figure S1
pointszoo <- predict(zglm, surv70)
View(pointszoo)
unname(pointszoo)
pps.z <- as.data.frame(cbind(pointszoo, surv70$Dist))
plot(pointszoo~V2, data=pps.z)

pointsp <- predict(prvl, surv70)
View(pointsp)
unname(pointsp)
pps.p <- as.data.frame(cbind(pointsp, surv70$Dist))
plot(pointsp~V2, data=pps.p)

figS1A <- ggplot(pps.z, aes(x=log10(V2+1), y=pointszoo))+
  geom_smooth(method=glm, color="black")+
  geom_point()+
  theme(axis.text=element_text(size=30),
        panel.background=element_blank(),
        panel.grid.major.y=element_line(color="grey", size=.2),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=40),
        legend.position="none")+
  ylab("log10(Zoospore equivalent)") +
  xlab("Distance (log10(km)) between \nhost & Bd collection sites")
figS1B <- ggplot(pps.p, aes(x=log10(V2+1), y=pointsp))+
  geom_smooth(method=glm, color="black")+
  geom_point()+
  theme(axis.text=element_text(size=30),
        panel.background=element_blank(),
        panel.grid.major.y=element_line(color="grey", size=.2),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=40),
        legend.position="none")+
  ylab("Probability of infection") +
  xlab("Distance (log10(km)) between \nhost & Bd collection sites")
ggarrange(figS1A, figS1B, nrow=1, ncol=2)

#Figure S2
prev.fig <- emmeans(prvl, list(pairwise ~ Bd), adjust = "tukey")
str(summary(prev.fig$`emmeans of Bd`))
prev.fig <- summary(prev.fig$`emmeans of Bd`)[1:6]
prev.fig$emmean <- log2prob(prev.fig$emmean)
prev.fig$asymp.LCL <- log2prob(prev.fig$asymp.LCL)
prev.fig$asymp.UCL <- log2prob(prev.fig$asymp.UCL)
plot(emmean~Bd, data=prev.fig)
prev.fig$tukey <- c("A,B","C", "A", "A,B", "B", "A,B")
prev.fig
figS2f <- ggplot(prev.fig, aes(x=Bd, y=emmean))+
  geom_point(size=3)+ylim(0, 1.1)+
  geom_errorbar(ymin=prev.fig$asymp.LCL, ymax=prev.fig$asymp.UCL, size=.3,width=.2)+
  theme(axis.text=element_text(size=15),
        panel.background=element_blank(),
        panel.grid.major.y=element_line(color="grey", size=.2),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=20),
        plot.title = element_text(size=25, hjust=0.5))+
  geom_text(data=prev.fig, vjust=0,
            aes(x=Bd, y=0.02+asymp.UCL,label=tukey))+
  labs(y="", x="", title= )
figS2f 

prev.fig2 <- emmeans(prvl, list(pairwise ~ Host), adjust = "tukey")
str(summary(prev.fig2$`emmeans of Host`))
prev.fig2 <- summary(prev.fig2$`emmeans of Host`)[1:6]
prev.fig2$emmean <- log2prob(prev.fig2$emmean)
prev.fig2$asymp.LCL <- log2prob(prev.fig2$asymp.LCL)
prev.fig2$asymp.UCL <- log2prob(prev.fig2$asymp.UCL)
plot(emmean~Host, data=prev.fig2)
prev.fig2$tukey <- c("A,B", "A,B", "A", "A", "B")
prev.fig2
figS2e <- ggplot(prev.fig2, aes(x=Host, y=emmean))+
  geom_point(size=3)+ylim(0, 1.1)+
  geom_errorbar(ymin=prev.fig2$asymp.LCL, ymax=prev.fig2$asymp.UCL, size=.3,width=.2)+
  theme(axis.text=element_text(size=15),
        panel.background=element_blank(),
        panel.grid.major.y=element_line(color="grey", size=.2),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=20),
        plot.title = element_text(size=25, hjust=0.5))+
  geom_text(data=prev.fig2, vjust=0,
            aes(x=Host, y=0.02+asymp.UCL,label=tukey))+
  labs(y="Probability of infection", x="", title= )
figS2e  

figS2c <- emmeans(zglm, list(pairwise ~ Host), adjust = "tukey")
str(summary(figS2c$`emmeans of Host`))
figS2c <- summary(figS2c$`emmeans of Host`)[1:6]
figS2c$emmean <- log2prob(figS2c$emmean)
figS2c$asymp.LCL <- log2prob(figS2c$asymp.LCL)
figS2c$asymp.UCL <- log2prob(figS2c$asymp.UCL)
plot(emmean~Host, data=figS2c)
figS2c$tukey <- c("A","A,B", "A,B", "A,B", "B")
figS2c

figS2d <- emmeans(zglm, list(pairwise ~ Bd), adjust = "tukey")
str(summary(figS2d$`emmeans of Bd`))
figS2d <- summary(figS2d$`emmeans of Bd`)[1:6]
figS2d$emmean <- log2prob(figS2d$emmean)
figS2d$asymp.LCL <- log2prob(figS2d$asymp.LCL)
figS2d$asymp.UCL <- log2prob(figS2d$asymp.UCL)
plot(emmean~Bd, data=figS2d)
figS2d$tukey <- c("A","B", "C", "A", "A", "B")
figS2d

S2a<-ggplot(figS2a, aes(x=Host, y=emmean))+
  geom_point(size=3)+ylim(0, 0.4)+
  geom_errorbar(ymin=figS2a$asymp.LCL, ymax=figS2a$asymp.UCL, size=.3,width=.2)+
  theme(axis.text=element_text(size=15),
        panel.background=element_blank(),
        panel.grid.major.y=element_line(color="grey", size=.2),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=20),
        plot.title = element_text(size=25, hjust=0.5))+
  geom_text(data=figS2a, vjust=0,
            aes(x=Host, y=0.02+asymp.UCL,label=tukey))+
  labs(y="Probability of mortality", x="", title= "Host identity")
S2b<-ggplot(figS2b, aes(x=Bd, y=emmean))+
  geom_point(size=3)+ylim(0, 0.4)+
  geom_errorbar(ymin=figS2b$asymp.LCL, ymax=figS2b$asymp.UCL, size=.3,width=.2)+
  theme(axis.text=element_text(size=15),
        panel.background=element_blank(),
        panel.grid.major.y=element_line(color="grey", size=.2),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=20),
        plot.title = element_text(size=25, hjust=0.5))+
  geom_text(data=figS2b, vjust=0,
            aes(x=Bd, y=0.02+asymp.UCL,label=tukey))+
  labs(y="", x="", title= "Bd identity")
S2c<-ggplot(figS2c, aes(x=Host, y=emmean))+
  geom_point(size=3)+ylim(0.4, 0.8)+
  geom_errorbar(ymin=figS2c$asymp.LCL, ymax=figS2c$asymp.UCL, size=.3,width=.2)+
  theme(axis.text=element_text(size=15),
        panel.background=element_blank(),
        panel.grid.major.y=element_line(color="grey", size=.2),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=20),
        plot.title = element_text(size=25, hjust=0.5))+
  geom_text(data=figS2c, vjust=0,
            aes(x=Host, y=0.02+asymp.UCL,label=tukey))+
  labs(y="log10(Zoospore equivalent)", x="", title= "")
S2d<-ggplot(figS2d, aes(x=Bd, y=emmean))+
  geom_point(size=3)+ylim(0.4, 0.8)+
  geom_errorbar(ymin=figS2d$asymp.LCL, ymax=figS2d$asymp.UCL, size=.3,width=.2)+
  theme(axis.text=element_text(size=15),
        panel.background=element_blank(),
        panel.grid.major.y=element_line(color="grey", size=.2),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=20),
        plot.title = element_text(size=25, hjust=0.5))+
  geom_text(data=figS2d, vjust=0,
            aes(x=Bd, y=0.02+asymp.UCL,label=tukey))+
  labs(y="", x="", title= "")
S2a;S2b;S2c;S2d
library(ggpubr)
ggarrange(S2a,S2b,S2c,S2d,figS2e,figS2f, ncol=2, nrow=3)

###### novel metric model comparisons #############
#survival
surv <-coxph(Surv(Date, mortality) ~   Host + Bd +mass+ log(Dist2+1), 
             data = surv70, na.action="na.fail")
surv.null <-coxph(Surv(Date, mortality) ~   Host + Bd +mass, 
             data = surv70, na.action="na.fail")
surv.binary <-coxph(Surv(Date, mortality) ~   Host + Bd +mass + LoN, 
                  data = surv70, na.action="na.fail")
surv.null <-coxph(Surv(Date, mortality) ~   Host + Bd +mass, 
                  data = surv70, na.action="na.fail")

surv70.sub <- subset(surv70, Bd != "QC")
surv.sub <-coxph(Surv(Date, mortality) ~   Host + Bd +mass +log(Dist2+1), 
                  data = surv70.sub, na.action="na.fail")
surv.host <-coxph(Surv(Date, mortality) ~   Host + Bd +mass +HostGenDist, 
                  data = surv70.sub, na.action="na.fail")
surv.IB <-coxph(Surv(Date, mortality) ~   Host + Bd +mass +IB, 
                  data = surv70.sub, na.action="na.fail")
surv.LB <-coxph(Surv(Date, mortality) ~   Host + Bd +mass +LB, 
                  data = surv70.sub, na.action="na.fail")
surv.EB <-coxph(Surv(Date, mortality) ~   Host + Bd +mass +EB, 
                  data = surv70.sub, na.action="na.fail")
AIC(surv.null);AIC(surv);AIC(surv.binary);
AIC(surv.sub);AIC(surv.host);AIC(surv.IB);
AIC(surv.LB);AIC(surv.EB)

#prevalence
prvl <- glm(inf.status ~ Host + Bd  + log(Dist2+1) + mass, data=surv70, family=binomial(link="logit"))
prvl.null <- glm(inf.status ~ Host + Bd  + mass, data=surv70, family=binomial(link="logit"))
prvl.binary <- glm(inf.status ~ Host + Bd  + mass +LoN, data=surv70, family=binomial(link="logit"))
prvl.sub <- glm(inf.status ~ Host + Bd  + log(Dist2+1) + mass, data=surv70.sub, family=binomial(link="logit"))
prvl.host <- glm(inf.status ~ Host + Bd  + HostGenDist + mass, data=surv70.sub, family=binomial(link="logit"))
prvl.IB <- glm(inf.status ~ Host + Bd  + IB + mass, data=surv70.sub, family=binomial(link="logit"))
prvl.EB <- glm(inf.status ~ Host + Bd  + EB + mass, data=surv70.sub, family=binomial(link="logit"))
prvl.LB <- glm(inf.status ~ Host + Bd  + LB + mass, data=surv70.sub, family=binomial(link="logit"))
AIC(prvl.null);AIC(prvl);AIC(prvl.binary);
AIC(prvl.sub);AIC(prvl.host);AIC(prvl.IB);
AIC(prvl.LB);AIC(prvl.EB)

# zoospore load
zglm <- glm(logZoospore ~ Host + Bd + mass + log(Dist2+1) + Date14, data=surv70, 
            na.action="na.fail")
zglm.null <- glm(logZoospore ~ Host + Bd + mass + log(Dist2+1) + Date14, data=surv70, 
            na.action="na.fail")
zglm.binar <- glm(logZoospore ~ Host + Bd + mass + log(Dist2+1) + Date14, data=surv70, 
            na.action="na.fail")
zglm.sub <- glm(logZoospore ~ Host + Bd + mass + log(Dist2+1) + Date14, data=surv70.sub, 
            na.action="na.fail")
zglm.host <- glm(logZoospore ~ Host + Bd + mass + HostGenDist + Date14, data=surv70.sub, 
            na.action="na.fail")
zglm.IB <- glm(logZoospore ~ Host + Bd + mass + IB + Date14, data=surv70.sub, 
            na.action="na.fail")
zglm.LB <- glm(logZoospore ~ Host + Bd + mass + LB + Date14, data=surv70.sub, 
            na.action="na.fail")
zglm.EB <- glm(logZoospore ~ Host + Bd + mass + EB + Date14, data=surv70.sub, 
            na.action="na.fail")

AIC(zglm.null);AIC(zglm);AIC(zglm.binary);
AIC(zglm.sub);AIC(zglm.host);AIC(zglm.IB);
AIC(zglm.LB);AIC(zglm.EB)


######### meta analysis ##############

library(metafor)
library(blme)
library(car)

bdma <- read.csv("Sauer et al. 2020 data subset.csv")

blme.ma <- blmer(ln.or ~ logdose + Superfam + logdistkm +
                    (1|genus.and.species)+(1|ID)+(1|bd.ID), 
                  data=bdma, weights = 1/v.ln.or, 
                  resid.prior = point(1.0), 
                  cov.prior = NULL, na.action="na.fail")
summary(blme.ma)
Anova(blme.ma)

