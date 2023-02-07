library(readr)
library(lmerTest)
library(sjPlot)
df<- read.csv("D:/Online_Trust_Game/Data_files/age_srs_lme_all.csv")
library(nlme)
df[,'gender']<-factor(df[,'gender'])
df[,'group']<-factor(df[,'group'])
df[,'sage']<-scale(df[,'age'], center = TRUE, scale = TRUE)
df[,'SRS']<-scale(df[,'SRS'], center = TRUE, scale = TRUE)

y <- scale(df[["SRS"]])
df[,'SRS'] <- y

#TD uncertainty age
model1 <- lm(pars1 ~ age*SRS, data=df)
tab_model(model1)
summary(model1)

model2 <- lm(pars2 ~ age*SRS , data=df)
tab_model(model2)
summary(model2)

model3 <- lm(pars3 ~ age*SRS, data=df)
tab_model(model3)
summary(model3)

model4 <- lm(pars4 ~ age*SRS, data=df)
tab_model(model4)
summary(model4)

model1 <- lm(BIC ~ group*age*SRS, data=df)
tab_model(model1)
summary(model1)

model2 <- lm(pars1 ~ group*age*SRS, data=df)
tab_model(model2)
summary(model2)

model2 <- lm(pars2 ~ group*age*SRS, data=df)
tab_model(model2)
summary(model2)

model3 <- lm(pars3 ~ group*age*SRS, data=df)
tab_model(model3)
summary(model3)

model4 <- lm(pars4 ~ group*age*SRS, data=df)
tab_model(model4)
summary(model4)

library(tidyverse)
library(extrafont)
df<- read.csv("D:/Online_Trust_Game/Data_files/age_srs_lme_TD_plot.csv")

ggplot(df,aes(x=SRS,y=pars3,color=young))+
  geom_point(size=5)+
  geom_smooth(method="lm",aes(fill=young),size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  theme(axis.line.x = element_line(color = 'black',size = 1),
        axis.line.y = element_line(color = 'black',size = 1))+
  scale_fill_manual(values=c("#999999", "#E69F00")) + 
  scale_color_manual(values=c("#999999", "#E69F00")) + 
  labs(fill="Age group", x="SRS raw", y="Prior belief")+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(-2,2))+
  coord_cartesian(ylim=c(-0.05,1.05),xlim=c(0,160))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))

ggplot(df,aes(x=age,y=pars1))+
  geom_point(size=5)+
  geom_smooth(method="lm",size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  theme(axis.line.x = element_line(color = 'black',size = 1),
        axis.line.y = element_line(color = 'black',size = 1))+
  scale_fill_manual(values=c("#999999", "#E69F00")) + 
  scale_color_manual(values=c("#999999", "#E69F00")) + 
  labs(x="Age", y="Decision noise")+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(-500,1600))+
  coord_cartesian(ylim=c(-0.05,1600),xlim=c(80,160))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))


library(ggplot2)
# Basic barplot
p <- ggplot(data = df, aes(x=groupplot, y=pars4)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=pars4-sd, ymax=pars4+sd), width=.2,
                position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") + theme_minimal()

p<-ggplot(data=df, aes(x=groupplot, y=pars4)) +
  geom_bar(stat="identity")
p

df<- read.csv("D:/Online_Trust_Game/Data_files/age_srs_lme_ASD_plot.csv")

ggplot(df,aes(x=SRS,y=pars3,color=young))+
  geom_point(size=5)+
  geom_smooth(method="lm",aes(fill=young),size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  theme(axis.line.x = element_line(color = 'black',size = 1),
        axis.line.y = element_line(color = 'black',size = 1))+
  scale_fill_manual(values=c("#999999", "#E69F00")) + 
  scale_color_manual(values=c("#999999", "#E69F00")) + 
  labs(fill="Age group", x="SRS raw", y="Prior belief")+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(-2,2))+
  coord_cartesian(ylim=c(-0.05,1.05),xlim=c(0,160))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))